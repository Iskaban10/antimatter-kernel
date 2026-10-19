#ifndef INTERACTION_FAST_MULTIPOLE_H_
#define INTERACTION_FAST_MULTIPOLE_H_

#include <cstddef>
#include <memory>
#include <vector>

#include "core/container/math_array.h"
#include "core/real_t.h"

namespace bdm 
{

    class ResourceManager;

    // Precomputed index tables for the expansions (defined in the .cc).
    struct FmmTables;

    // Fast multipole method (FMM) for the Coulomb force between point charges.

    //                         HOW IT WORKS
    // Particles are sorted into an adaptive octree. Every cell carries two
    // Cartesian Taylor expansions of the potential 1/|x - y|:
    //   - a MULTIPOLE expansion of the charges inside the cell (valid far
    //     away from the cell), and
    //   - a LOCAL expansion of the potential produced by all distant cells
    //     (valid close to the cell's center).
    // A dual-tree traversal pairs a target cell with a source cell:
    //   - if the two are well separated, the source's multipole expansion is
    //     converted into the target's local expansion (M2L);
    //   - otherwise, if both are leaves, their particles interact directly
    //     (P2P); otherwise the larger cell is split and the test is repeated.
    // Local expansions are then pushed down to the leaves (L2L) and evaluated
    // at each particle (L2P). The cost is O(N) for fixed accuracy, versus
    // O(N^2) for the direct sum and O(N log N) for Barnes-Hut.
    //
    // Two settings control accuracy:
    //   - `order` P: expansions keep all terms up to total degree P. The
    //     error falls roughly geometrically as P increases.
    //   - `theta`: two cells are well separated when
    //     (radius_target + radius_source) < theta * distance_between_centers.
    //     Smaller theta is more accurate and slower. theta = 0 never accepts
    //     an expansion, so everything is summed directly and the result is
    //     the exact all-pairs force (a validation reference).



    class FastMultipoleSolver 
    {
        public:

            struct Source 
            {
                Real3 position;
                real_t charge = 0.0;
            };

            FastMultipoleSolver(int order, real_t theta, size_t leaf_capacity = 32);

            // Computes the Coulomb force (Newtons) on every source from all the others. `forces` is resized to sources.size()
            void Compute(const std::vector<Source>& sources, std::vector<Real3>* forces) const;

        private:
            int order_ = 4;
            real_t theta_ = 0.5;
            size_t leaf_capacity_ = 32;
            std::shared_ptr<const FmmTables> tables_;
    };
    void ComputeFastMultipoleCoulombForces(ResourceManager* rm, int order, real_t theta);
}
#endif