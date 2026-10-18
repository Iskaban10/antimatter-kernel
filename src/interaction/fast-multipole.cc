#include "interaction/fast-multipole.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

#include "agents/charged_antiparticle.h"
#include "core/resource_manager.h"
#include "utils/constants.h"

namespace bdm
{
    namespace
    {
        constexpr int kMaxDepth = 32;
        constexpr double kMinSquaredDistance = 1e-30;
        constexpr double kSqrt3 = 1.7320508075688772;

        inline int Octant(const std::array<double, 3>& p, const std::array<double, 3>& c)
        {
            return (p[0] >= c[0] ? 1 : 0) | (p[1] >= c[1] ? 2 : 0) | (p[2] >= c[2] ? 4 : 0);
        }
    }

    FmmTree::FmmTree(int order, real_t theta, size_t leaf_capacity)
        : p_(order < 1 ? 1 : order),
          squared_theta_(static_cast<double>(theta) * static_cast<double>(theta)),
          leaf_capacity_(leaf_capacity < 1 ? 1 : leaf_capacity)
    {
        BuildTables();
    }

    // Multi-indices are ordered by total degree, so the first n_p_ entries are
    // exactly those with |alpha| <= p, and any recurrence on lower degree only
    // reads entries that are already computed.
    void FmmTree::BuildTables()
    {
        const int P2 = 2 * p_;
        s_ = P2 + 1;
        index_.assign(static_cast<size_t>(s_) * s_ * s_, -1);
        mi_.clear();
        for (int n = 0; n <= P2; ++n)
        {
            for (int a = n; a >= 0; --a)
            {
                for (int b = n - a; b >= 0; --b)
                {
                    const int c = n - a - b;
                    index_[(a * s_ + b) * s_ + c] = static_cast<int>(mi_.size());
                    mi_.push_back({a, b, c});
                }
            }
        }
        n_2p_ = static_cast<int>(mi_.size());
        n_p_ = (p_ + 1) * (p_ + 2) * (p_ + 3) / 6;

        std::vector<double> fact(P2 + 1, 1.0);
        for (int i = 1; i <= P2; ++i) fact[i] = fact[i - 1] * i;

        inv_fact_.resize(n_2p_);
        for (int k = 0; k < n_2p_; ++k)
        {
            inv_fact_[k] = 1.0 / (fact[mi_[k][0]] * fact[mi_[k][1]] * fact[mi_[k][2]]);
        }

        sign_.resize(n_p_);
        for (int k = 0; k < n_p_; ++k)
        {
            sign_[k] = ((mi_[k][0] + mi_[k][1] + mi_[k][2]) & 1) ? -1.0 : 1.0;
        }

        sum_idx_.resize(static_cast<size_t>(n_p_) * n_p_);
        for (int beta = 0; beta < n_p_; ++beta)
        {
            for (int alpha = 0; alpha < n_p_; ++alpha)
            {
                sum_idx_[static_cast<size_t>(beta) * n_p_ + alpha] =
                    Idx(mi_[alpha][0] + mi_[beta][0], mi_[alpha][1] + mi_[beta][1], mi_[alpha][2] + mi_[beta][2]);
            }
        }
    }

    // out[k] = d0^a * d1^b * d2^c for the first `count` multi-indices.
    void FmmTree::Powers(const Vec3& d, int count, double* out) const
    {
        out[0] = 1.0;
        for (int k = 1; k < count; ++k)
        {
            std::array<int, 3> m = mi_[k];
            const int dim = m[0] > 0 ? 0 : (m[1] > 0 ? 1 : 2);
            --m[dim];
            out[k] = out[Idx(m)] * d[dim];
        }
    }

    // D[k] = partial derivative of 1/|x| with multi-index mi_[k], at x = R.
    //
    // Recurrence (d = first nonzero component of n, m = n - e_d), derived from
    // r^2 * d_d(1/r) + x_d * (1/r) = 0 by differentiating m times:
    //   r^2 D^n = - x_d D^{n-e_d}
    //             - sum_k 2 m_k x_k D^{n-e_k}
    //             - sum_k m_k (m_k - 1) D^{n-2e_k}
    //             - m_d D^{n-2e_d}
    void FmmTree::Derivatives(const Vec3& R, double* D) const
    {
        const double r2 = R[0] * R[0] + R[1] * R[1] + R[2] * R[2];
        const double inv_r2 = 1.0 / r2;
        D[0] = std::sqrt(inv_r2);
        for (int k = 1; k < n_2p_; ++k)
        {
            const std::array<int, 3>& n = mi_[k];
            const int d = n[0] > 0 ? 0 : (n[1] > 0 ? 1 : 2);
            std::array<int, 3> m = n;
            --m[d];

            double acc = R[d] * D[Idx(m)];
            for (int kk = 0; kk < 3; ++kk)
            {
                if (m[kk] <= 0) continue;
                std::array<int, 3> t = n;
                --t[kk];
                acc += 2.0 * m[kk] * R[kk] * D[Idx(t)];
                if (m[kk] > 1)
                {
                    --t[kk];
                    acc += static_cast<double>(m[kk]) * (m[kk] - 1) * D[Idx(t)];
                }
            }
            if (m[d] > 0)
            {
                std::array<int, 3> t = n;
                t[d] -= 2;
                acc += static_cast<double>(m[d]) * D[Idx(t)];
            }
            D[k] = -acc * inv_r2;
        }
    }

    // parent_alpha += sum_{gamma <= alpha} child_gamma * d^{alpha-gamma} / (alpha-gamma)!
    // d = child_center - parent_center.
    void FmmTree::M2M(double* parent, const double* child, const Vec3& d, double* w) const
    {
        Powers(d, n_p_, w);
        for (int k = 0; k < n_p_; ++k) w[k] *= inv_fact_[k];

        for (int k = 0; k < n_p_; ++k)
        {
            const auto& al = mi_[k];
            double s = 0.0;
            for (int g0 = 0; g0 <= al[0]; ++g0)
                for (int g1 = 0; g1 <= al[1]; ++g1)
                    for (int g2 = 0; g2 <= al[2]; ++g2)
                    {
                        s += child[Idx(g0, g1, g2)] * w[Idx(al[0] - g0, al[1] - g1, al[2] - g2)];
                    }
            parent[k] += s;
        }
    }

    // child_delta += sum_{beta >= delta} parent_beta * e^{beta-delta} / (beta-delta)!
    // e = child_center - parent_center.
    void FmmTree::L2L(const double* parent, double* child, const Vec3& e, double* w) const
    {
        Powers(e, n_p_, w);
        for (int k = 0; k < n_p_; ++k) w[k] *= inv_fact_[k];

        for (int k = 0; k < n_p_; ++k)
        {
            const auto& be = mi_[k];
            const double lb = parent[k];
            for (int d0 = 0; d0 <= be[0]; ++d0)
                for (int d1 = 0; d1 <= be[1]; ++d1)
                    for (int d2 = 0; d2 <= be[2]; ++d2)
                    {
                        child[Idx(d0, d1, d2)] += lb * w[Idx(be[0] - d0, be[1] - d1, be[2] - d2)];
                    }
        }
    }

    void FmmTree::Build(const std::vector<Source>& sources)
    {
        const size_t n = sources.size();
        pos_.resize(n);
        q_.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            pos_[i] = {static_cast<double>(sources[i].position[0]),
                       static_cast<double>(sources[i].position[1]),
                       static_cast<double>(sources[i].position[2])};
            q_[i] = static_cast<double>(sources[i].charge);
        }

        nodes_.clear();
        perm_.resize(n);
        scratch_.resize(n);
        std::iota(perm_.begin(), perm_.end(), size_t{0});
        if (n == 0) return;

        Vec3 low = pos_[0];
        Vec3 high = pos_[0];
        for (const auto& x : pos_)
        {
            for (int d = 0; d < 3; ++d)
            {
                low[d] = std::min(low[d], x[d]);
                high[d] = std::max(high[d], x[d]);
            }
        }
        const Vec3 center{0.5 * (low[0] + high[0]), 0.5 * (low[1] + high[1]), 0.5 * (low[2] + high[2])};
        const double extent = std::max({high[0] - low[0], high[1] - low[1], high[2] - low[2]});
        double half = 0.5 * extent * (1.0 + 1e-9);
        if (half < 1e-30) half = 1e-30;

        nodes_.reserve(2 * n / leaf_capacity_ + 16);
        BuildNode(0, n, center, half, 0);
    }

    int FmmTree::BuildNode(size_t begin, size_t end, const Vec3& center, double half, int depth)
    {
        const int index = static_cast<int>(nodes_.size());
        nodes_.emplace_back();
        nodes_[index].center = center;
        nodes_[index].half_size = half;
        nodes_[index].begin = begin;
        nodes_[index].end = end;

        if (end - begin <= leaf_capacity_ || depth >= kMaxDepth)
        {
            nodes_[index].is_leaf = true;
            return index;
        }
        nodes_[index].is_leaf = false;

        std::array<size_t, 8> counts{};
        for (size_t k = begin; k < end; ++k) ++counts[Octant(pos_[perm_[k]], center)];

        std::array<size_t, 9> offsets{};
        offsets[0] = begin;
        for (int o = 0; o < 8; ++o) offsets[o + 1] = offsets[o] + counts[o];

        std::array<size_t, 8> cursor{};
        for (int o = 0; o < 8; ++o) cursor[o] = offsets[o];
        for (size_t k = begin; k < end; ++k)
        {
            const int o = Octant(pos_[perm_[k]], center);
            scratch_[cursor[o]++] = perm_[k];
        }
        std::copy(scratch_.begin() + begin, scratch_.begin() + end, perm_.begin() + begin);

        const double child_half = 0.5 * half;
        for (int o = 0; o < 8; ++o)
        {
            if (counts[o] == 0) continue;
            const Vec3 child_center{center[0] + ((o & 1) ? child_half : -child_half),
                                    center[1] + ((o & 2) ? child_half : -child_half),
                                    center[2] + ((o & 4) ? child_half : -child_half)};
            const int child = BuildNode(offsets[o], offsets[o + 1], child_center, child_half, depth + 1);
            nodes_[index].children[o] = child;
        }
        return index;
    }

    // "Targets in a receive from sources in b."
    void FmmTree::Traverse(int a, int b)
    {
        const Node& A = nodes_[a];
        const Node& B = nodes_[b];

        const double dx = A.center[0] - B.center[0];
        const double dy = A.center[1] - B.center[1];
        const double dz = A.center[2] - B.center[2];
        const double dist2 = dx * dx + dy * dy + dz * dz;
        const double rsum = kSqrt3 * (A.half_size + B.half_size);

        if (a != b && rsum * rsum < squared_theta_ * dist2)
        {
            m2l_[a].push_back(b);
            return;
        }
        if (A.is_leaf && B.is_leaf)
        {
            p2p_[a].push_back(b);
            return;
        }
        if (a == b)
        {
            for (int i = 0; i < 8; ++i)
            {
                if (A.children[i] < 0) continue;
                for (int j = 0; j < 8; ++j)
                {
                    if (A.children[j] >= 0) Traverse(A.children[i], A.children[j]);
                }
            }
            return;
        }
        if (B.is_leaf || (!A.is_leaf && A.half_size >= B.half_size))
        {
            for (int o = 0; o < 8; ++o)
            {
                if (A.children[o] >= 0) Traverse(A.children[o], b);
            }
        }
        else
        {
            for (int o = 0; o < 8; ++o)
            {
                if (B.children[o] >= 0) Traverse(a, B.children[o]);
            }
        }
    }

    void FmmTree::Solve()
    {
        const size_t nn = nodes_.size();
        const size_t n = pos_.size();
        field_.assign(n, Vec3{0.0, 0.0, 0.0});
        if (nn == 0) return;

        mult_.assign(nn * n_p_, 0.0);
        loc_.assign(nn * n_p_, 0.0);
        m2l_.assign(nn, std::vector<int>());
        p2p_.assign(nn, std::vector<int>());
        Traverse(0, 0);

        const int64_t NN = static_cast<int64_t>(nn);

        // P2M: multipole of each leaf about its center
        #pragma omp parallel
        {
            std::vector<double> pw(n_p_);
            #pragma omp for schedule(dynamic, 16)
            for (int64_t a = 0; a < NN; ++a)
            {
                const Node& nd = nodes_[a];
                if (!nd.is_leaf) continue;
                double* M = &mult_[a * n_p_];
                for (size_t k = nd.begin; k < nd.end; ++k)
                {
                    const size_t j = perm_[k];
                    const Vec3 b{pos_[j][0] - nd.center[0], pos_[j][1] - nd.center[1], pos_[j][2] - nd.center[2]};
                    Powers(b, n_p_, pw.data());
                    for (int m = 0; m < n_p_; ++m) M[m] += q_[j] * pw[m] * inv_fact_[m];
                }
            }
        }

        // M2M: children always have a larger index than their parent
        {
            std::vector<double> w(n_p_);
            for (int64_t a = NN - 1; a >= 0; --a)
            {
                const Node& nd = nodes_[a];
                if (nd.is_leaf) continue;
                for (int o = 0; o < 8; ++o)
                {
                    const int c = nd.children[o];
                    if (c < 0) continue;
                    const Vec3 d{nodes_[c].center[0] - nd.center[0], nodes_[c].center[1] - nd.center[1],
                                 nodes_[c].center[2] - nd.center[2]};
                    M2M(&mult_[a * n_p_], &mult_[static_cast<size_t>(c) * n_p_], d, w.data());
                }
            }
        }

        // M2L: each target node writes only its own local expansion
        #pragma omp parallel
        {
            std::vector<double> D(n_2p_);
            #pragma omp for schedule(dynamic, 16)
            for (int64_t a = 0; a < NN; ++a)
            {
                double* L = &loc_[a * n_p_];
                for (int b : m2l_[a])
                {
                    const Vec3 R{nodes_[a].center[0] - nodes_[b].center[0], nodes_[a].center[1] - nodes_[b].center[1],
                                 nodes_[a].center[2] - nodes_[b].center[2]};
                    Derivatives(R, D.data());
                    const double* M = &mult_[static_cast<size_t>(b) * n_p_];
                    for (int beta = 0; beta < n_p_; ++beta)
                    {
                        double s = 0.0;
                        const int* si = &sum_idx_[static_cast<size_t>(beta) * n_p_];
                        for (int alpha = 0; alpha < n_p_; ++alpha) s += sign_[alpha] * M[alpha] * D[si[alpha]];
                        L[beta] += s;
                    }
                }
            }
        }

        // L2L: parents come before children, so a parent's local is final when pushed
        {
            std::vector<double> w(n_p_);
            for (int64_t a = 0; a < NN; ++a)
            {
                const Node& nd = nodes_[a];
                if (nd.is_leaf) continue;
                for (int o = 0; o < 8; ++o)
                {
                    const int c = nd.children[o];
                    if (c < 0) continue;
                    const Vec3 e{nodes_[c].center[0] - nd.center[0], nodes_[c].center[1] - nd.center[1],
                                 nodes_[c].center[2] - nd.center[2]};
                    L2L(&loc_[a * n_p_], &loc_[static_cast<size_t>(c) * n_p_], e, w.data());
                }
            }
        }

        // L2P and P2P: each leaf writes only the field of its own particles
        #pragma omp parallel
        {
            std::vector<double> w(n_p_);
            #pragma omp for schedule(dynamic, 16)
            for (int64_t a = 0; a < NN; ++a)
            {
                const Node& nd = nodes_[a];
                if (!nd.is_leaf) continue;
                const double* L = &loc_[a * n_p_];

                for (size_t k = nd.begin; k < nd.end; ++k)
                {
                    const size_t i = perm_[k];
                    const Vec3 x{pos_[i][0] - nd.center[0], pos_[i][1] - nd.center[1], pos_[i][2] - nd.center[2]};
                    Powers(x, n_p_, w.data());
                    for (int m = 0; m < n_p_; ++m) w[m] *= inv_fact_[m];

                    double g[3] = {0.0, 0.0, 0.0};  // gradient of the potential
                    for (int beta = 0; beta < n_p_; ++beta)
                    {
                        for (int d = 0; d < 3; ++d)
                        {
                            if (mi_[beta][d] == 0) continue;
                            std::array<int, 3> t = mi_[beta];
                            --t[d];
                            g[d] += L[beta] * w[Idx(t)];
                        }
                    }
                    field_[i][0] -= g[0];
                    field_[i][1] -= g[1];
                    field_[i][2] -= g[2];
                }

                for (int b : p2p_[a])
                {
                    const Node& nb = nodes_[b];
                    for (size_t k = nd.begin; k < nd.end; ++k)
                    {
                        const size_t i = perm_[k];
                        double ex = 0.0, ey = 0.0, ez = 0.0;
                        for (size_t m = nb.begin; m < nb.end; ++m)
                        {
                            const size_t j = perm_[m];
                            if (j == i) continue;
                            const double qj = q_[j];
                            if (qj == 0.0) continue;
                            const double rx = pos_[i][0] - pos_[j][0];
                            const double ry = pos_[i][1] - pos_[j][1];
                            const double rz = pos_[i][2] - pos_[j][2];
                            const double r2 = rx * rx + ry * ry + rz * rz;
                            if (r2 < kMinSquaredDistance) continue;
                            const double f = qj / (r2 * std::sqrt(r2));
                            ex += f * rx;
                            ey += f * ry;
                            ez += f * rz;
                        }
                        field_[i][0] += ex;
                        field_[i][1] += ey;
                        field_[i][2] += ez;
                    }
                }
            }
        }
    }

    Real3 FmmTree::GetForce(size_t i) const
    {
        const double s = static_cast<double>(constants::kCoulombConstant) * q_[i];
        return Real3{static_cast<real_t>(s * field_[i][0]), static_cast<real_t>(s * field_[i][1]),
                     static_cast<real_t>(s * field_[i][2])};
    }

    void ComputeFmmCoulombForces(ResourceManager* rm, int order, real_t theta)
    {
        std::vector<ChargedAntiparticle*> particles;
        rm->ForEachAgent([&](Agent* agent)
        {
            auto* particle = FastCast(agent);
            if (particle != nullptr && !particle->IsAnnihilated()) particles.push_back(particle);
        });

        if (particles.empty()) return;

        std::vector<FmmTree::Source> sources(particles.size());
        for (size_t i = 0; i < particles.size(); ++i)
        {
            sources[i].position = particles[i]->GetPosition();
            sources[i].charge = particles[i]->GetCharge();
        }

        FmmTree tree(order, theta);
        tree.Build(sources);
        tree.Solve();

        const int64_t n = static_cast<int64_t>(particles.size());
        #pragma omp parallel for schedule(static)
        for (int64_t i = 0; i < n; ++i)
        {
            particles[static_cast<size_t>(i)]->SetCoulombForce(tree.GetForce(static_cast<size_t>(i)));
        }
    }
}