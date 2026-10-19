#include "interaction/fast-multipole.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <mutex>

#include "agents/charged_antiparticle.h"
#include "core/resource_manager.h"
#include "utils/constants.h"

namespace bdm {

// ---------------------------------------------------------------------------
// Math used throughout (all in rescaled coordinates, see Compute()).
//
// Multi-indices alpha = (a, b, c) label monomials x^a y^b z^c.
// With f(r) = 1/|r| and D_alpha(r) = d^alpha f(r):
//
//   multipole moments of a cell centered at c:
//       M_alpha = sum_i q_i (y_i - c)^alpha / alpha!
//   potential of the cell at x:
//       phi(x)  = sum_alpha (-1)^|alpha| M_alpha D_alpha(x - c)
//
//   M2M (child center c_k -> parent center c_p, d = c_k - c_p):
//       M^p_alpha = sum_{beta <= alpha} M^k_beta d^(alpha-beta) / (alpha-beta)!
//   M2L (source center c -> target center t, r = t - c):
//       L_gamma   = 1/gamma! sum_alpha (-1)^|alpha| M_alpha D_{alpha+gamma}(r)
//       phi(x)    = sum_gamma L_gamma (x - t)^gamma
//   L2L (parent center t_p -> child center t_c, s = t_c - t_p):
//       L^c_beta  = sum_{gamma >= beta} L^p_gamma C(gamma, beta) s^(gamma-beta)
//   L2P: the field is E = -grad phi, with
//       d_k phi   = sum_gamma gamma_k L_gamma (x - t)^(gamma - e_k)
//
// Derivatives of 1/r come from a recurrence derived from r^2 d_i f = -x_i f.
// ---------------------------------------------------------------------------

struct FmmTables {
  struct Triple {
    int a;
    int b;
    int c;
  };
  struct Recipe {
    int axis = 0;
    int delta[3] = {0, 0, 0};
    int minus1[3] = {-1, -1, -1};  // index of delta - e_j, or -1
    int minus2[3] = {-1, -1, -1};  // index of delta - 2 e_j, or -1
  };

  int order = 0;
  int n_m = 0;  // multi-indices with degree <= P
  int n_d = 0;  // multi-indices with degree <= 2P
  int size = 0;  // 2P + 1
  std::vector<std::array<int, 3>> multi;  // sorted by total degree
  std::vector<int> lookup;
  std::vector<double> fact, inv_fact, sign;
  std::vector<Recipe> recipes;
  std::vector<int> m2l_index;                  // [gamma * n_m + alpha]
  std::vector<Triple> m2m_pairs, l2l_pairs;
  std::vector<std::array<int, 3>> grad_index;  // index of gamma - e_k, or -1

  int Index(int a, int b, int c) const {
    if (a < 0 || b < 0 || c < 0 || a >= size || b >= size || c >= size) {
      return -1;
    }
    return lookup[(a * size + b) * size + c];
  }

  explicit FmmTables(int p) : order(p), size(2 * p + 1) {
    const int max_degree = 2 * p;
    lookup.assign(static_cast<size_t>(size) * size * size, -1);
    for (int degree = 0; degree <= max_degree; ++degree) {
      for (int a = degree; a >= 0; --a) {
        for (int b = degree - a; b >= 0; --b) {
          const int c = degree - a - b;
          lookup[(a * size + b) * size + c] = static_cast<int>(multi.size());
          multi.push_back({a, b, c});
        }
      }
    }
    n_d = static_cast<int>(multi.size());
    n_m = (p + 1) * (p + 2) * (p + 3) / 6;

    std::vector<double> factorial(max_degree + 1, 1.0);
    for (int i = 1; i <= max_degree; ++i) {
      factorial[i] = factorial[i - 1] * i;
    }
    fact.resize(n_d);
    inv_fact.resize(n_d);
    sign.resize(n_d);
    for (int k = 0; k < n_d; ++k) {
      const auto& m = multi[k];
      fact[k] = factorial[m[0]] * factorial[m[1]] * factorial[m[2]];
      inv_fact[k] = 1.0 / fact[k];
      sign[k] = ((m[0] + m[1] + m[2]) % 2 == 0) ? 1.0 : -1.0;
    }

    recipes.resize(n_d);
    for (int k = 1; k < n_d; ++k) 
    {
      const auto& m = multi[k];
      Recipe& r = recipes[k];
      r.axis = m[0] > 0 ? 0 : (m[1] > 0 ? 1 : 2);
      for (int j = 0; j < 3; ++j) {
        r.delta[j] = m[j];
        std::array<int, 3> t = m;
        if (m[j] >= 1) {
          t[j] -= 1;
          r.minus1[j] = Index(t[0], t[1], t[2]);
        }
        if (m[j] >= 2) {
          t[j] -= 1;
          r.minus2[j] = Index(t[0], t[1], t[2]);
        }
      }
    }

    m2l_index.resize(static_cast<size_t>(n_m) * n_m);
    for (int g = 0; g < n_m; ++g) 
    {
      for (int a = 0; a < n_m; ++a) 
      {
        m2l_index[static_cast<size_t>(g) * n_m + a] = Index(multi[a][0] + multi[g][0], multi[a][1] + multi[g][1], multi[a][2] + multi[g][2]);
      }
    }

    for (int a = 0; a < n_m; ++a) 
    {
      for (int b = 0; b < n_m; ++b) 
      {
        if (multi[b][0] <= multi[a][0] && multi[b][1] <= multi[a][1] && multi[b][2] <= multi[a][2]) 
        {
          m2m_pairs.push_back({a, b, Index(multi[a][0] - multi[b][0], multi[a][1] - multi[b][1], multi[a][2] - multi[b][2])});
        }
        if (multi[a][0] <= multi[b][0] && multi[a][1] <= multi[b][1] && multi[a][2] <= multi[b][2]) 
        {
          l2l_pairs.push_back({a, b, Index(multi[b][0] - multi[a][0], multi[b][1] - multi[a][1], multi[b][2] - multi[a][2])});
        }
      }
    }

    grad_index.resize(n_m);
    for (int g = 0; g < n_m; ++g) 
    {
      for (int k = 0; k < 3; ++k) 
      {
        std::array<int, 3> t = multi[g];
        if (t[k] >= 1) 
        {
          t[k] -= 1;
          grad_index[g][k] = Index(t[0], t[1], t[2]);
        } 
        else 
        {
          grad_index[g][k] = -1;
        }
      }
    }
  }
};

namespace 
{

constexpr int kMaxDepth = 32;
constexpr int kMinOrder = 1;
constexpr int kMaxOrder = 8;
// Pairs closer than this (m^2, physical units) are skipped.
constexpr double kMinSquaredDistance = 1e-30;

using Tables = FmmTables;

// out[k] = x^a y^b z^c for every multi-index of degree <= P.
void Monomials(const Tables& t, double x, double y, double z, double* out) 
{
  double px[kMaxOrder + 1], py[kMaxOrder + 1], pz[kMaxOrder + 1];
  px[0] = py[0] = pz[0] = 1.0;
  for (int i = 1; i <= t.order; ++i) 
  {
    px[i] = px[i - 1] * x;
    py[i] = py[i - 1] * y;
    pz[i] = pz[i - 1] * z;
  }
  for (int k = 0; k < t.n_m; ++k) 
  {
    out[k] = px[t.multi[k][0]] * py[t.multi[k][1]] * pz[t.multi[k][2]];
  }
}

// D[k] = d^alpha (1/|r|) for every multi-index of degree <= 2P.
void Derivatives(const Tables& t, double x, double y, double z, double* D) 
{
  const double inv_r2 = 1.0 / (x * x + y * y + z * z);
  D[0] = std::sqrt(inv_r2);
  const double xs[3] = {x, y, z};
  for (int k = 1; k < t.n_d; ++k) 
  {
    const auto& rc = t.recipes[k];
    const int i = rc.axis;
    double sum = (2 * rc.delta[i] - 1) * xs[i] * D[rc.minus1[i]];
    if (rc.minus2[i] >= 0) 
    {
      sum += (rc.delta[i] - 1) * (rc.delta[i] - 1) * D[rc.minus2[i]];
    }
    for (int j = 0; j < 3; ++j) 
    {
      if (j == i || rc.delta[j] == 0) continue;
      sum += 2.0 * rc.delta[j] * xs[j] * D[rc.minus1[j]];
      if (rc.minus2[j] >= 0) 
      {
        sum += rc.delta[j] * (rc.delta[j] - 1) * D[rc.minus2[j]];
      }
    }
    D[k] = -sum * inv_r2;
  }
}

inline int Octant(double x, double y, double z, double cx, double cy, double cz) 
{
  return (x >= cx ? 1 : 0) | (y >= cy ? 2 : 0) | (z >= cz ? 4 : 0);
}

struct Cell 
{
  double cx = 0, cy = 0, cz = 0;  // geometric center
  double half = 0;                // half edge length
  double radius = 0;              // max distance from center to a particle
  int children[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
  size_t begin = 0, end = 0;  // particle range in sorted order
  bool is_leaf = true;
};

// Everything built per Compute() call.
struct Workspace 
{
  std::vector<double> x, y, z, q;  // rescaled input, original order
  std::vector<size_t> order, scratch;
  std::vector<Cell> cells;
  size_t leaf_capacity = 32;

  int BuildCell(size_t begin, size_t end, double cx, double cy, double cz, double half, int depth) {
    const int index = static_cast<int>(cells.size());
    cells.emplace_back();
    cells[index].cx = cx;
    cells[index].cy = cy;
    cells[index].cz = cz;
    cells[index].half = half;
    cells[index].begin = begin;
    cells[index].end = end;

    if (end - begin <= leaf_capacity || depth >= kMaxDepth) 
    {
      cells[index].is_leaf = true;
      return index;
    }
    cells[index].is_leaf = false;

    std::array<size_t, 8> counts{};
    for (size_t k = begin; k < end; ++k) 
    {
      const size_t p = order[k];
      ++counts[Octant(x[p], y[p], z[p], cx, cy, cz)];
    }
    std::array<size_t, 9> offsets{};
    offsets[0] = begin;
    for (int o = 0; o < 8; ++o) {
      offsets[o + 1] = offsets[o] + counts[o];
    }
    std::array<size_t, 8> cursor{};
    for (int o = 0; o < 8; ++o) 
    {
      cursor[o] = offsets[o];
    }
    for (size_t k = begin; k < end; ++k) 
    {
      const size_t p = order[k];
      scratch[cursor[Octant(x[p], y[p], z[p], cx, cy, cz)]++] = p;
    }
    std::copy(scratch.begin() + begin, scratch.begin() + end,
              order.begin() + begin);

    const double ch = 0.5 * half;
    for (int o = 0; o < 8; ++o) 
    {
      if (counts[o] == 0)   continue;
      
      const int child = BuildCell(
          offsets[o], offsets[o + 1], cx + ((o & 1) ? ch : -ch),
          cy + ((o & 2) ? ch : -ch), cz + ((o & 4) ? ch : -ch), ch,
          depth + 1);
      cells[index].children[o] = child;
    }
    return index;
  }
};

std::shared_ptr<const FmmTables> GetTables(int order) 
{
  static std::mutex mutex;
  static std::shared_ptr<const FmmTables> cache[kMaxOrder + 1];
  std::lock_guard<std::mutex> lock(mutex);
  if (!cache[order])    cache[order] = std::make_shared<const FmmTables>(order);
  return cache[order];
}

}  // namespace

FastMultipoleSolver::FastMultipoleSolver(int order, real_t theta, size_t leaf_capacity)
                        : order_(std::min(std::max(order, kMinOrder), kMaxOrder)), theta_(theta),
                        leaf_capacity_(leaf_capacity < 1 ? 1 : leaf_capacity), tables_(GetTables(order_)) {}

void FastMultipoleSolver::Compute(const std::vector<Source>& sources, std::vector<Real3>* forces) const 
{
  const size_t n = sources.size();
  forces->assign(n, Real3{0.0, 0.0, 0.0});
  if (n == 0)   return;
  const FmmTables& t = *tables_;
  const int n_m = t.n_m;

  double low[3], high[3];
  for (int d = 0; d < 3; ++d) 
  {
    low[d] = high[d] = static_cast<double>(sources[0].position[d]);
  }
  for (const auto& s : sources) 
  {
    for (int d = 0; d < 3; ++d) 
    {
      const double v = static_cast<double>(s.position[d]);
      low[d] = std::min(low[d], v);
      high[d] = std::max(high[d], v);
    }
  }
  const double extent = std::max({high[0] - low[0], high[1] - low[1], high[2] - low[2]});
  const double length = std::max(extent, 1e-15);  // physical cube edge
  const double mid[3] = {0.5 * (low[0] + high[0]), 0.5 * (low[1] + high[1]), 0.5 * (low[2] + high[2])};
  
  const double force_scale = constants::kCoulombConstant / (length * length);
  const double min_sq = kMinSquaredDistance / (length * length);

  Workspace ws;
  ws.leaf_capacity = leaf_capacity_;
  ws.x.resize(n);
  ws.y.resize(n);
  ws.z.resize(n);
  ws.q.resize(n);
  for (size_t i = 0; i < n; ++i) {
    ws.x[i] = (static_cast<double>(sources[i].position[0]) - mid[0]) / length;
    ws.y[i] = (static_cast<double>(sources[i].position[1]) - mid[1]) / length;
    ws.z[i] = (static_cast<double>(sources[i].position[2]) - mid[2]) / length;
    ws.q[i] = static_cast<double>(sources[i].charge);
  }

  // --- Tree ----------------------------------------------------------------
  ws.order.resize(n);
  ws.scratch.resize(n);
  for (size_t i = 0; i < n; ++i) {
    ws.order[i] = i;
  }
  ws.cells.reserve(2 * n / std::max<size_t>(leaf_capacity_, 1) + 16);
  ws.BuildCell(0, n, 0.0, 0.0, 0.0, 0.5 * (1.0 + 1e-9), 0);
  std::vector<Cell>& cells = ws.cells;
  const int n_cells = static_cast<int>(cells.size());

  // Particle data in tree order, so each leaf is a contiguous range.
  std::vector<double> sx(n), sy(n), sz(n), sq(n);
  for (size_t k = 0; k < n; ++k) 
  {
    const size_t p = ws.order[k];
    sx[k] = ws.x[p];
    sy[k] = ws.y[p];
    sz[k] = ws.z[p];
    sq[k] = ws.q[p];
  }

  // Cell radii: leaves from their particles; parents bound their children.
  for (int c = n_cells - 1; c >= 0; --c) 
  {
    Cell& cell = cells[c];
    double r2 = 0.0;
    if (cell.is_leaf) 
    {
      for (size_t k = cell.begin; k < cell.end; ++k) 
      {
        const double dx = sx[k] - cell.cx, dy = sy[k] - cell.cy, dz = sz[k] - cell.cz;
        r2 = std::max(r2, dx * dx + dy * dy + dz * dz);
      }
      cell.radius = std::sqrt(r2);
    } 
    else 
    {
      double r = 0.0;
      for (int o = 0; o < 8; ++o) 
      {
        if (cell.children[o] < 0) continue;
        const Cell& ch = cells[cell.children[o]];
        const double dx = ch.cx - cell.cx, dy = ch.cy - cell.cy, dz = ch.cz - cell.cz;
        r = std::max(r, std::sqrt(dx * dx + dy * dy + dz * dz) + ch.radius);
      }
      cell.radius = r;
    }
  }

  std::vector<double> M(static_cast<size_t>(n_cells) * n_m, 0.0);
  std::vector<double> L(static_cast<size_t>(n_cells) * n_m, 0.0);

  // --- P2M (leaves, parallel) ------------------------------------------------
#pragma omp parallel
  {
    std::vector<double> mono(n_m);
#pragma omp for schedule(dynamic, 16)
    for (int c = 0; c < n_cells; ++c) {
      const Cell& cell = cells[c];
      if (!cell.is_leaf) {
        continue;
      }
      double* m = &M[static_cast<size_t>(c) * n_m];
      for (size_t k = cell.begin; k < cell.end; ++k) {
        if (sq[k] == 0.0) {
          continue;
        }
        Monomials(t, sx[k] - cell.cx, sy[k] - cell.cy, sz[k] - cell.cz,
                  mono.data());
        for (int a = 0; a < n_m; ++a) {
          m[a] += sq[k] * mono[a] * t.inv_fact[a];
        }
      }
    }
  }

  // --- M2M (children have larger indices than their parent) -----------------
  {
    std::vector<double> mono(n_m);
    for (int c = n_cells - 1; c >= 0; --c) 
    {
      const Cell& cell = cells[c];
      if (cell.is_leaf)  continue;
      double* mp = &M[static_cast<size_t>(c) * n_m];
      for (int o = 0; o < 8; ++o) 
      {
        const int child = cell.children[o];
        if (child < 0)  continue;
        const Cell& ch = cells[child];
        Monomials(t, ch.cx - cell.cx, ch.cy - cell.cy, ch.cz - cell.cz, mono.data());
        const double* mc = &M[static_cast<size_t>(child) * n_m];
        for (const auto& pr : t.m2m_pairs) 
        {
          mp[pr.a] += mc[pr.b] * mono[pr.c] * t.inv_fact[pr.c];
        }
      }
    }
  }

  // --- Dual-tree traversal: build interaction lists (serial) ----------------
  std::vector<std::vector<int>> m2l_list(n_cells), p2p_list(n_cells);
  const double theta = static_cast<double>(theta_);
  const double theta2 = theta * theta;
  {
    std::vector<std::pair<int, int>> stack;
    stack.emplace_back(0, 0);
    while (!stack.empty()) {
      const auto [ti, si] = stack.back();
      stack.pop_back();
      const Cell& T = cells[ti];
      const Cell& S = cells[si];
      const double dx = T.cx - S.cx, dy = T.cy - S.cy, dz = T.cz - S.cz;
      const double dist2 = dx * dx + dy * dy + dz * dz;
      const double rsum = T.radius + S.radius;
      if (rsum * rsum < theta2 * dist2) 
      {
        m2l_list[ti].push_back(si);
        continue;
      }
      if (T.is_leaf && S.is_leaf) 
      {
        p2p_list[ti].push_back(si);
        continue;
      }
      const bool split_target = S.is_leaf || (!T.is_leaf && T.radius >= S.radius);
      if (split_target) 
      {
        for (int o = 0; o < 8; ++o) 
        {
          if (T.children[o] >= 0)   stack.emplace_back(T.children[o], si);
        }
      } 
      else 
      {
        for (int o = 0; o < 8; ++o) 
        {
          if (S.children[o] >= 0)   stack.emplace_back(ti, S.children[o]);
        }
      }
    }
  }

  // --- M2L (parallel over target cells) -----------------------------------
#pragma omp parallel
  {
    std::vector<double> D(t.n_d);
#pragma omp for schedule(dynamic, 8)
    for (int ti = 0; ti < n_cells; ++ti) 
    {
      if (m2l_list[ti].empty())    continue;
      double* l = &L[static_cast<size_t>(ti) * n_m];
      for (const int si : m2l_list[ti]) 
      {
        const Cell& T = cells[ti];
        const Cell& S = cells[si];
        Derivatives(t, T.cx - S.cx, T.cy - S.cy, T.cz - S.cz, D.data());
        const double* m = &M[static_cast<size_t>(si) * n_m];
        for (int g = 0; g < n_m; ++g) 
        {
          const int* idx = &t.m2l_index[static_cast<size_t>(g) * n_m];
          double acc = 0.0;
          for (int a = 0; a < n_m; ++a) 
          {
            acc += t.sign[a] * m[a] * D[idx[a]];
          }
          l[g] += acc * t.inv_fact[g];
        }
      }
    }
  }

  // --- L2L (parents have smaller indices than their children) ---------------
  {
    std::vector<double> mono(n_m);
    for (int c = 0; c < n_cells; ++c) 
    {
      const Cell& cell = cells[c];
      if (cell.is_leaf)     continue;
      const double* lp = &L[static_cast<size_t>(c) * n_m];
      for (int o = 0; o < 8; ++o) 
      {
        const int child = cell.children[o];
        if (child < 0)  continue;
        const Cell& ch = cells[child];
        Monomials(t, ch.cx - cell.cx, ch.cy - cell.cy, ch.cz - cell.cz, mono.data());
        double* lc = &L[static_cast<size_t>(child) * n_m];
        for (const auto& pr : t.l2l_pairs) {
          lc[pr.a] += lp[pr.b] * t.fact[pr.b] * t.inv_fact[pr.a] *
                      t.inv_fact[pr.c] * mono[pr.c];
        }
      }
    }
  }

  // --- L2P + P2P (parallel over leaves) -------------------------------------
  // Field per unit target charge, rescaled units: E_i = sum_j q_j (x_i - x_j)/r^3
  std::vector<double> ex(n, 0.0), ey(n, 0.0), ez(n, 0.0);
#pragma omp parallel
  {
    std::vector<double> mono(n_m);
#pragma omp for schedule(dynamic, 8)
    for (int ti = 0; ti < n_cells; ++ti) 
    {
      const Cell& T = cells[ti];
      if (!T.is_leaf)     continue;
      const double* l = &L[static_cast<size_t>(ti) * n_m];
      for (size_t i = T.begin; i < T.end; ++i) 
      {
        Monomials(t, sx[i] - T.cx, sy[i] - T.cy, sz[i] - T.cz, mono.data());
        double gx = 0.0, gy = 0.0, gz = 0.0;
        for (int g = 1; g < n_m; ++g) 
        {
          const auto& gi = t.grad_index[g];
          if (gi[0] >= 0) gx += t.multi[g][0] * l[g] * mono[gi[0]];
          if (gi[1] >= 0) gy += t.multi[g][1] * l[g] * mono[gi[1]];
          if (gi[2] >= 0) gz += t.multi[g][2] * l[g] * mono[gi[2]];
        }
        double fx = -gx, fy = -gy, fz = -gz;

        // P2P
        for (const int si : p2p_list[ti]) 
        {
          const Cell& S = cells[si];
          for (size_t j = S.begin; j < S.end; ++j) 
          {
            if (j == i || sq[j] == 0.0)     continue;
            const double rx = sx[i] - sx[j], ry = sy[i] - sy[j],
                         rz = sz[i] - sz[j];
            const double r2 = rx * rx + ry * ry + rz * rz;
            if (r2 < min_sq)    continue;
            const double f = sq[j] / (r2 * std::sqrt(r2));
            fx += rx * f;
            fy += ry * f;
            fz += rz * f;
          }
        }
        ex[i] = fx;
        ey[i] = fy;
        ez[i] = fz;
      }
    }
  }

  // --- Back to physical forces, in the caller's order ---------------------
  for (size_t k = 0; k < n; ++k) 
  {
    const size_t p = ws.order[k];
    const double s = force_scale * sq[k];
    (*forces)[p] = Real3{static_cast<real_t>(s * ex[k]),  static_cast<real_t>(s * ey[k]), static_cast<real_t>(s * ez[k])};
  }
}

void ComputeFastMultipoleCoulombForces(ResourceManager* rm, int order, real_t theta) 
{
  std::vector<ChargedAntiparticle*> charged;
  rm->ForEachAgent([&](Agent* agent) 
  {
    auto* particle = FastCast(agent);
    if (particle == nullptr || particle->IsAnnihilated())  return;
    if (particle->GetCharge() == 0.0) 
    {
      particle->SetCoulombForce(Real3{0.0, 0.0, 0.0});  // neutral: no force
    } 
    else 
    {
      charged.push_back(particle);
    }
  });
  if (charged.empty())  return;

  std::vector<FastMultipoleSolver::Source> sources(charged.size());
  for (size_t i = 0; i < charged.size(); ++i) 
  {
    sources[i].position = charged[i]->GetPosition();
    sources[i].charge = charged[i]->GetCharge();
  }

  FastMultipoleSolver solver(order, theta);
  std::vector<Real3> forces;
  solver.Compute(sources, &forces);

  for (size_t i = 0; i < charged.size(); ++i) 
  {
    charged[i]->SetCoulombForce(forces[i]);
  }
}
}