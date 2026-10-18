#ifndef INTERACTION_BARNES_HUT_H_
#define INTERACTION_BARNES_HUT_H_

#include <cstddef>
#include <vector>

#include "core/container/math_array.h"
#include "core/real_t.h"

namespace bdm 
{
    class ResourceManager;

    // Barnes-Hut octree for the Coulomb force between point charges.
    //
    //                          HOW IT WORKS
    // Particles are sorted into an octree. When computing the force on one
    // particle, a tree node that is "far enough away" is replaced by point
    // charges at its charge centers (a monopole approximation); otherwise the
    // node is opened and its children are examined. A node of width s at
    // distance d from the target is accepted when s / d < theta.
    //   1> theta = 0 never accepts a node, so every pair is summed directly:
    //     the exact all-pairs Coulomb force, with NO cutoff radius. Useful as
    //     a reference for validating larger theta values.
    //   2> Larger theta is faster and less accurate. 0.3 to 0.7 is the usual
    //     range for monopole-only Barnes-Hut.

    class BarnesHutTree 
    {
        public:

        struct Source 
        {
            Real3 position;
            real_t charge = 0.0;
        };

        // Passed as `self_index` when the target is not one of the sources.
        static constexpr size_t kNoSelf = static_cast<size_t>(-1);

        explicit BarnesHutTree(real_t theta, size_t leaf_capacity = 8);
        void Build(const std::vector<Source>& sources);

        Real3 ComputeForce(const Real3& position, real_t charge, size_t self_index = kNoSelf) const;

        size_t NumNodes() const { return nodes_.size(); }

        // private:
        struct Node 
        {
            Real3 center = {0.0, 0.0, 0.0};  // geometric center of the cube
            real_t half_size = 0.0;          // half of the cube's edge length
            real_t positive_charge = 0.0;    // sum of positive charges
            real_t negative_charge = 0.0;    // sum of |negative charges| (>= 0)
            Real3 positive_center = {0.0, 0.0, 0.0};
            Real3 negative_center = {0.0, 0.0, 0.0};
            int children[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
            size_t begin = 0;  // particle range [begin, end) within order_
            size_t end = 0;
            bool is_leaf = true;
        };

        int BuildNode(size_t begin, size_t end, const Real3& center, real_t half, int depth);
        void Accumulate(int node_index, const Real3& position, real_t charge, size_t self_index, Real3* force) const;

        real_t theta_ = 0.5;
        real_t squared_theta_ = 0.25;
        size_t leaf_capacity_ = 8;
        std::vector<Source> sources_;
        std::vector<size_t> order_;
        std::vector<size_t> scratch_;
        std::vector<Node> nodes_;
    };

    void ComputeBarnesHutCoulombForces(ResourceManager* rm, real_t theta);
}
#endif