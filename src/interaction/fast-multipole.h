#ifndef INTERACTION_FMM_H_
#define INTERACTION_FMM_H_

#include <array>
#include <cstddef>
#include <vector>

#include "core/container/math_array.h"
#include "core/real_t.h"

namespace bdm
{
    class ResourceManager;

    // Fast Multipole Method for the Coulomb force between point charges.
    //
    //                          HOW IT WORKS
    // Particles are sorted into an octree. Every cell stores a multipole
    // expansion (Cartesian Taylor moments, order p) of the charges inside it,
    // and a local expansion of the field produced by all well-separated cells.
    //   P2M : particles -> multipole of a leaf
    //   M2M : child multipole -> parent multipole        (upward pass)
    //   M2L : multipole of a far cell -> local of a target cell
    //   L2L : parent local -> child local                 (downward pass)
    //   L2P : local of a leaf -> field at its particles
    //   P2P : direct sum between neighbouring leaves
    // Two cells A, B interact through M2L when
    //   sqrt(3) * (halfA + halfB) < theta * |centerA - centerB|
    // otherwise the larger one is opened (or P2P if both are leaves).
    //   - Cost is O(N) for fixed order and theta (BH is O(N log N)).
    //   - Error falls roughly like theta^(order+1). Because moments are signed,
    //     error is measured against sum|q|, not the net charge of a cell.
    //   - Larger order or smaller theta = more accurate and slower.

    class FmmTree
    {
    public:
        struct Source
        {
            Real3 position;
            real_t charge = 0.0;
        };

        explicit FmmTree(int order = 4, real_t theta = 0.5, size_t leaf_capacity = 32);

        void Build(const std::vector<Source>& sources);
        void Solve();

        // Coulomb force (including the Coulomb constant) on source i.
        // Valid after Solve().
        Real3 GetForce(size_t i) const;

        size_t NumNodes() const { return nodes_.size(); }

    // private:
        using Vec3 = std::array<double, 3>;

        struct Node
        {
            Vec3 center = {0.0, 0.0, 0.0};
            double half_size = 0.0;
            int children[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
            size_t begin = 0;  // particle range [begin, end) within perm_
            size_t end = 0;
            bool is_leaf = true;
        };

        void BuildTables();
        int BuildNode(size_t begin, size_t end, const Vec3& center, double half, int depth);
        void Traverse(int a, int b);

        int Idx(int a, int b, int c) const { return index_[(a * s_ + b) * s_ + c]; }
        int Idx(const std::array<int, 3>& m) const { return Idx(m[0], m[1], m[2]); }
        void Powers(const Vec3& d, int count, double* out) const;
        void Derivatives(const Vec3& R, double* D) const;
        void M2M(double* parent, const double* child, const Vec3& d, double* w) const;
        void L2L(const double* parent, double* child, const Vec3& e, double* w) const;

        // settings
        int p_ = 4;
        double squared_theta_ = 0.25;
        size_t leaf_capacity_ = 32;

        // multi-index tables (|alpha| <= 2p)
        int s_ = 0;       // 2p + 1
        int n_p_ = 0;     // number of multi-indices with |alpha| <= p
        int n_2p_ = 0;    // number of multi-indices with |alpha| <= 2p
        std::vector<int> index_;
        std::vector<std::array<int, 3>> mi_;
        std::vector<double> inv_fact_;  // 1 / (a! b! c!)
        std::vector<double> sign_;      // (-1)^|alpha|, size n_p_
        std::vector<int> sum_idx_;      // [beta * n_p_ + alpha] -> index of alpha + beta

        // particles
        std::vector<Vec3> pos_;
        std::vector<double> q_;
        std::vector<size_t> perm_;
        std::vector<size_t> scratch_;
        std::vector<Vec3> field_;  // electric field with k = 1

        // tree and expansions
        std::vector<Node> nodes_;
        std::vector<double> mult_;  // nodes * n_p_
        std::vector<double> loc_;   // nodes * n_p_
        std::vector<std::vector<int>> m2l_;
        std::vector<std::vector<int>> p2p_;
    };

    void ComputeFmmCoulombForces(ResourceManager* rm, int order, real_t theta);
}

#endif  // INTERACTION_FMM_H_