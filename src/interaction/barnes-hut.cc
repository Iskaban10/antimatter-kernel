#include "interaction/barnes-hut.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

#include "agents/charged_antiparticle.h"
#include "core/resource_manager.h"
#include "utils/constants.h"

namespace bdm 
{
    namespace 
    {
        constexpr int kMaxDepth = 32;   // deepest allowed subdivision
        constexpr real_t kMinSquaredDistance = 1e-30;

        // adds Coulomb force on a charge `target_q` at `target_pos` due to point charge `source_q` at `source_pos`.
        inline void AddPointCharge(const Real3& target_pos, real_t target_q, const Real3& source_pos, real_t source_q, Real3* force) 
        {
            const Real3 r = target_pos - source_pos;
            const real_t squared_distance = r * r;
            if (squared_distance < kMinSquaredDistance) return;
            const real_t inverse_distance = 1.0 / std::sqrt(squared_distance);
            const real_t magnitude = constants::kCoulombConstant * target_q * source_q / squared_distance;
            *force = *force + r * (magnitude * inverse_distance);
        }

        // Index (0-7) of the octant of `position` relative to `center`.
        inline int Octant(const Real3& position, const Real3& center) 
        {
            return (position[0] >= center[0] ? 1 : 0) | (position[1] >= center[1] ? 2 : 0) | (position[2] >= center[2] ? 4 : 0);
        }
    }

    BarnesHutTree::BarnesHutTree(real_t theta, size_t leaf_capacity) : theta_(theta), squared_theta_(theta * theta), leaf_capacity_(leaf_capacity < 1 ? 1 : leaf_capacity) {}

    void BarnesHutTree::Build(const std::vector<Source>& sources) 
    {
        sources_ = sources;
        nodes_.clear();
        order_.resize(sources_.size());
        scratch_.resize(sources_.size());
        for (size_t i = 0; i < order_.size(); ++i) order_[i] = i;
        if (sources_.empty()) return;

        // Smallest cube containing every particle.
        Real3 low = sources_[0].position;
        Real3 high = sources_[0].position;
        for (const auto& s : sources_) 
        {
            for (int d = 0; d < 3; ++d) 
            {
            low[d] = std::min(low[d], s.position[d]);
            high[d] = std::max(high[d], s.position[d]);
            }
        }
        const Real3 center{0.5 * (low[0] + high[0]), 0.5 * (low[1] + high[1]), 0.5 * (low[2] + high[2])};
        const real_t extent = std::max( {high[0] - low[0], high[1] - low[1], high[2] - low[2]});
        real_t half = 0.5 * extent * (1.0 + 1e-9);

        if (half < 1e-30) half = 1e-30;  // all particles at one point

        nodes_.reserve(sources_.size());
        BuildNode(0, sources_.size(), center, half, 0);
    }

    int BarnesHutTree::BuildNode(size_t begin, size_t end, const Real3& center,real_t half, int depth) 
    {
        const int index = static_cast<int>(nodes_.size());
        nodes_.emplace_back();
        nodes_[index].center = center;
        nodes_[index].half_size = half;
        nodes_[index].begin = begin;
        nodes_[index].end = end;

        real_t positive_charge = 0.0;
        real_t negative_charge = 0.0;
        Real3 positive_moment{0.0, 0.0, 0.0};  // sum of q * position
        Real3 negative_moment{0.0, 0.0, 0.0};  // sum of |q| * position

        const size_t count = end - begin;
        if (count <= leaf_capacity_ || depth >= kMaxDepth) 
        {
            nodes_[index].is_leaf = true;
            for (size_t k = begin; k < end; ++k) 
            {
            const Source& s = sources_[order_[k]];
            if (s.charge > 0.0) 
            {
                positive_charge += s.charge;
                positive_moment = positive_moment + s.position * s.charge;
            } 
            else if (s.charge < 0.0) 
            {
                negative_charge += -s.charge;
                negative_moment = negative_moment + s.position * (-s.charge);
            }
            }
        } 
        else 
        {
            nodes_[index].is_leaf = false;

            std::array<size_t, 8> counts{};
            for (size_t k = begin; k < end; ++k) ++counts[Octant(sources_[order_[k]].position, center)];

            std::array<size_t, 9> offsets{};
            offsets[0] = begin;

            for (int o = 0; o < 8; ++o) offsets[o + 1] = offsets[o] + counts[o];

            std::array<size_t, 8> cursor{};

            for (int o = 0; o < 8; ++o) cursor[o] = offsets[o];

            for (size_t k = begin; k < end; ++k) 
            {
                const int o = Octant(sources_[order_[k]].position, center);
                scratch_[cursor[o]++] = order_[k];
            }
            std::copy(scratch_.begin() + begin, scratch_.begin() + end, order_.begin() + begin);

            const real_t child_half = 0.5 * half;
            for (int o = 0; o < 8; ++o) 
            {
                if (counts[o] == 0) continue;
                const Real3 child_center
                {
                    center[0] + ((o & 1) ? child_half : -child_half),
                    center[1] + ((o & 2) ? child_half : -child_half),
                    center[2] + ((o & 4) ? child_half : -child_half)
                };
                const int child = BuildNode(offsets[o], offsets[o + 1], child_center, child_half, depth + 1);
                nodes_[index].children[o] = child;
            }

            // Combines the children's charges and charge centers.
            for (int o = 0; o < 8; ++o) 
            {
                const int c = nodes_[index].children[o];
                if (c < 0) continue;
                
                const Node& child = nodes_[c];
                if (child.positive_charge > 0.0) 
                {
                    positive_charge += child.positive_charge;
                    positive_moment = positive_moment + child.positive_center * child.positive_charge;
                }
                if (child.negative_charge > 0.0) 
                {
                    negative_charge += child.negative_charge;
                    negative_moment = negative_moment + child.negative_center * child.negative_charge;
                }
            }
        }

        Node& node = nodes_[index];
        node.positive_charge = positive_charge;
        node.negative_charge = negative_charge;
        node.positive_center = positive_charge > 0.0 ? positive_moment * (1.0 / positive_charge) : center;
        node.negative_center = negative_charge > 0.0 ? negative_moment * (1.0 / negative_charge) : center;
        return index;
    }

    Real3 BarnesHutTree::ComputeForce(const Real3& position, real_t charge, size_t self_index) const 
    {
        Real3 force{0.0, 0.0, 0.0};
        if (nodes_.empty() || charge == 0.0) return force;

        Accumulate(0, position, charge, self_index, &force);
        return force;
    }

    void BarnesHutTree::Accumulate(int node_index, const Real3& position, real_t charge, size_t self_index, Real3* force) const 
    {
        const Node& node = nodes_[node_index];
        if (node.positive_charge == 0.0 && node.negative_charge == 0.0) return;

        if (node.is_leaf) 
        {
            for (size_t k = node.begin; k < node.end; ++k) 
            {
            const size_t j = order_[k];
            if (j == self_index) continue;
            const Source& s = sources_[j];
            if (s.charge == 0.0) continue;
            AddPointCharge(position, charge, s.position, s.charge, force);
            }
            return;
        }

        const Real3 delta = position - node.center;
        const real_t squared_distance = delta * delta;
        const real_t width = 2.0 * node.half_size;
        if (width * width < squared_theta_ * squared_distance) 
        {
            if (node.positive_charge > 0.0) 
            {
            AddPointCharge(position, charge, node.positive_center, node.positive_charge, force);
            }
            if (node.negative_charge > 0.0) 
            {
            AddPointCharge(position, charge, node.negative_center, -node.negative_charge, force);
            }
            return;
        }

        for (int o = 0; o < 8; ++o) 
        {
            if (node.children[o] >= 0) { Accumulate(node.children[o], position, charge, self_index, force); }
        }
    }

    void ComputeBarnesHutCoulombForces(ResourceManager* rm, real_t theta) 
    {
        std::vector<ChargedAntiparticle*> particles;
        rm->ForEachAgent([&](Agent* agent) 
        {
            auto* particle = FastCast(agent);
            if (particle != nullptr && !particle->IsAnnihilated()) particles.push_back(particle);
        });

        if (particles.empty()) return;

        std::vector<BarnesHutTree::Source> sources(particles.size());
        for (size_t i = 0; i < particles.size(); ++i) 
        {
            sources[i].position = particles[i]->GetPosition();
            sources[i].charge = particles[i]->GetCharge();
        }

        BarnesHutTree tree(theta);
        tree.Build(sources);

        const int64_t n = static_cast<int64_t>(particles.size());
        #pragma omp parallel for schedule(static)
        for (int64_t i = 0; i < n; ++i) 
        {
            const size_t idx = static_cast<size_t>(i);
            particles[idx]->SetCoulombForce(tree.ComputeForce(sources[idx].position, sources[idx].charge, idx));
        }
    }
}