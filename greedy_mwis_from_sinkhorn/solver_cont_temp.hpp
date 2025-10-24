#ifndef LIBMPOPT_MWIS_SOLVER_CONT_TEMP_HPP
#define LIBMPOPT_MWIS_SOLVER_CONT_TEMP_HPP

#include <random>
#include <iostream>
#include <vector>
#include <iterator>

namespace mwis
{

    // Can be switch to float for possibly faster dual updates. Loss of precision
    // will not affect the solution quality, because the exp-costs are used to
    // compute a equivalance transformation (reparametrization) at full precsision.
    using cost_exp = float;
    using index = unsigned int;
    using cost = double;

    constexpr const cost epsilon = 1e-8;
    constexpr const cost infinity = std::numeric_limits<cost>::infinity();
    static_assert(std::numeric_limits<cost>::has_infinity);
    static_assert(std::numeric_limits<cost>::has_signaling_NaN);
    static_assert(std::numeric_limits<cost>::has_quiet_NaN);

    struct range
    {
        index begin;
        index end;
        index size;
    };

    constexpr cost_exp default_stabilization_threshold = 10e30;
    constexpr cost_exp default_sparsity_threshold = 1e-8;
    constexpr bool initial_reparametrization = true;

    template <typename T>
    bool feasibility_check(const T sum) { return std::abs(sum - 1.0) < 1e-6; }
    template <>
    bool feasibility_check(const int sum) { return sum == 1; }

    class solver
    {
    public:
        solver()
            : finalized_graph_(false), finalized_costs_(false), constant_(0.0), scaling_(1.0), gen_(42)
#ifdef ENABLE_QPBO
              ,
              qpbo_(0, 0)
#endif
              ,
              temperature_drop_factor_(0.5)
        {
#ifndef ENABLE_QPBO
            // std::cerr << "!!!!!!!!!!\n"
            //           << "ENABLE_QPBO was not activated during configuration of libmpopt.\n"
            //           << "No fusion moves are performed and the the quality of the computed assignment is degraded.\n"
            //           << "!!!!!!!!!!\n" << std::endl;
#endif
        }

        index add_node(cost cost)
        {
            assert(!finalized_graph_);
            assert(no_cliques() == 0);
            costs_.push_back(cost / scaling_);
            orig_.push_back(cost);
            return costs_.size() - 1;
        }

        index add_clique(const std::vector<index> &indices)
        {
            assert(!finalized_graph_);

            range cl;
            cl.begin = clique_index_data_.size();
            cl.size = indices.size() + 1;
            cl.end = cl.begin + cl.size;
            clique_indices_.push_back(cl);

            for (auto index : indices)
                clique_index_data_.push_back(index);
            clique_index_data_.push_back(costs_.size());
            costs_.push_back(0.0);

            assert(cl.end == clique_index_data_.size());
            return clique_indices_.size() - 1;
        }

        void finalize()
        {
            temperature_ = 1;
            finalize_graph();
            finalize_costs();
        }

        bool finalized() const { return finalized_graph_ && finalized_costs_; }

        cost constant() const { return constant_ * scaling_; }
        void constant(cost c) { constant_ = c / scaling_; }

        template <bool reduced = false>
        cost node_cost(index i) const
        {
            assert(finalized_graph_);
            assert(i < no_nodes() && i < no_orig());
            return reduced ? costs_[i] * scaling_ : orig_[i];
        }

        void node_cost(index i, cost c)
        {
            assert(finalized_graph_);
            assert(i < no_nodes() && i < no_orig());
            const auto shift = c - orig_[i];
            orig_[i] += shift;
            costs_[i] += shift / scaling_;
            finalized_costs_ = false;
        }

        template <bool reduced = false>
        cost clique_cost(index i) const
        {
            assert(finalized());
            assert(i < no_cliques());
            const auto j = no_orig() + i;
            assert(j < no_nodes());
            return reduced ? costs_[j] * scaling_ : 0.0;
        }

        cost dual_relaxed() const { return dual_relaxed_scaled() * scaling_; }

        cost primal() const { return value_best_ * scaling_; }

        cost primal_relaxed() const { return value_relaxed_best_ * scaling_; }
        cost primal_relaxed_projected() const { return value_relaxed_ * scaling_; }

        template <typename OUTPUT_ITERATOR>
        void assignment(OUTPUT_ITERATOR begin, OUTPUT_ITERATOR end) const
        {
            assert(finalized_graph_);
            assert(end - begin == orig_.end() - orig_.begin());
            assert(assignment_best_.size() >= orig_.size());
            std::copy(assignment_best_.begin(), assignment_best_.begin() + orig_.size(), begin);
        }

        bool assignment(index node_idx) const
        {
            assert(finalized_graph_);
            assert(node_idx >= 0 && node_idx < orig_.size());
            return assignment_best_[node_idx];
        }

        /* Outputs the best found relaxed assignment for the slack-padded MWIS formulation */
        template <typename OUTPUT_ITERATOR>
        void assignment_relaxed(OUTPUT_ITERATOR begin, OUTPUT_ITERATOR end) const
        {
            assert(finalized_graph_);
            assert(end - begin == no_nodes());
            assert(assignment_relaxed_best_.size() == no_nodes());
            std::copy(assignment_relaxed_best_.cbegin(), assignment_relaxed_best_.cend(), begin);
        }

        /* Outputs the latest relaxed assignment found by the truncation projection for the slack-padded MWIS formulation */
        template <typename OUTPUT_ITERATOR>
        void assignment_relaxed_projected(OUTPUT_ITERATOR begin, OUTPUT_ITERATOR end) const
        {
            assert(finalized_graph_);
            assert(end - begin == no_nodes());
            assert(assignment_relaxed_best_.size() == no_nodes());
            std::copy(assignment_relaxed_.cbegin(), assignment_relaxed_.cend(), begin);
        }

        /* Outputs the costs for the slack-padded MWIS formulation */
        template <typename OUTPUT_ITERATOR>
        void reduced_costs(OUTPUT_ITERATOR begin, OUTPUT_ITERATOR end) const
        {
            assert(finalized_graph_);
            assert(end - begin == no_nodes());
            assert(costs_.size() == no_nodes());
            OUTPUT_ITERATOR begin_copy = begin;
            std::copy(costs_.cbegin(), costs_.cend(), begin);
            std::transform(begin_copy, end, begin_copy,
                           [this](cost x)
                           { return x * scaling_; });
        }

        int iterations() const { return iterations_; }

        void dualRun(const int batch_size,
                     const int max_batches,
                     const double min_duality_gap)
        {
            const std::string red = "\033[31m";
            const std::string blue = "\033[34m";
            const std::string green = "\033[32m";
            const std::string reset = "\033[0m"; // Reset to default color
            assert(finalized());

            // std::cout << green
            //           << "Dual Solver for Max-Weight Independent Set Problem based on the Bregman-Sinkhorn smoothed dual coordinate descent." << std::endl
            //           << red << "No primal heuristic" << green << " is computed, only dual updates and truncation projection to the feasible set of the LP relaxation." << std::endl
            //           << "Randomized greedy solutions can be generated additionally."
            //           << reset << std::endl;

            // std::cout << std::endl
            //           << "# nodes = " << no_nodes() << std::endl
            //           << "# cliques = " << no_cliques() << std::endl
            //           << "# dual iterations per batch = " << batch_size << std::endl
            //           << "max. # batches = " << max_batches << std::endl
            //           << "integer primal interations starts only after reaching " << min_duality_gap << "% relaxed duality gap" << std::endl
            //           << std::endl;

            for (int i = 0; i < max_batches; ++i)
            {
                init_exponential_domain();

                reduced_init_cliques();
                stabilization_counter_ = 0;
                int stabilization_iter_counter = 0;
                int max_stabilization_counter = 0;
                for (int j = 0; j < batch_size; ++j)
                {
                    reduced_single_pass();
                    if (stabilization_counter_ > 0 && j < (batch_size - 1)) // stabilization needed and not the last iteration
                    {
                        reparametrize();
                        init_exponential_domain();
                        reduced_init_cliques<false>(); // non-verbose
                        max_stabilization_counter = std::max(max_stabilization_counter, stabilization_counter_);
                        stabilization_counter_ = 0;
                        ++stabilization_iter_counter;
                    }
                }

                // if (stabilization_iter_counter > 0)
                // {
                //     std::cout
                //         << std::endl
                //         << red
                //         << "Stabilization was called "
                //         << stabilization_iter_counter
                //         << " times for at most "
                //         << max_stabilization_counter
                //         << " variables."
                //         << reset << std::endl;
                // }

                // Reparametrize alpha_ into costs_.
                reparametrize();

                // Second, we compute the primal projection by using costs_.
                compute_relaxed_truncated_projection();

                // Last, we can run greedy primal generator and fusion if the relaxed duality gap is reached.
                const auto d = dual_relaxed();
                auto p_relaxed = primal_relaxed();
                auto gap_relaxed = (d - p_relaxed) / d * 100.0;
                // t_dual.stop();

                if (gap_relaxed < min_duality_gap)
                {
                    // std::cout << "Relaxed duality gap = " << gap_relaxed << "< 0.1%, dual optimization stopped." << std::endl;
                    break;
                    // OUTPUT greedy solutions in JSON format: end ===============================
                }
                // else
                // {
                //     std::cout << "Relaxed duality gap limit of " << std::fixed << min_duality_gap << "% not reached, primal heuristic skipped." << std::endl;
                // }

                iterations_ += batch_size;

                // std::cout.precision(std::numeric_limits<cost>::max_digits10);
                p_relaxed = primal_relaxed(); // it can be eventually updated by the primal heuristic
                gap_relaxed = (d - p_relaxed) / d * 100.0;
                const auto p = primal();
                const auto gap = (d - p) / d * 100.0;

                // std::cout << std::endl
                //           << "it=" << iterations_ << " "
                //           << "d=" << d << " "
                //           << "p=" << p << " "
                //           << "gap=" << gap << "% "
                //           << "p_relaxed=" << p_relaxed << " "
                //           << "gap_relaxed=" << gap_relaxed << "% "
                //           << "T=" << temperature_ << " "
                //           << std::endl;

                // std::cout << "Leading constant for thresholding =" << gap_relaxed / temperature_ / no_nodes() << std::endl;

                // Now, we can update temperature (uses costs_ and primal projection result).
                update_temperature();
            }

            const auto d = dual_relaxed();
            const auto p = primal_relaxed();
            // std::cout.precision(std::numeric_limits<cost>::max_digits10);
            // std::cout << std::endl
            //           << std::endl
            //           << "Final optimization results: "
            //           << "primal relaxed=" << p << " "
            //           << "dual bound=" << d << " "
            //           << "duality gap=" << (d - p) / d * 100.0 << "% "
            //           << std::endl;
        }

        index no_nodes() const { return costs_.size(); }
        index no_orig() const { return orig_.size(); }
        index no_cliques() const { return clique_indices_.size(); }

        double temperature_drop_factor() const { return temperature_drop_factor_; }
        void temperature_drop_factor(const double v) { temperature_drop_factor_ = v; }

        double temperature() const { return temperature_; }
        void temperature(const double v) { temperature_ = v; }

        template <typename OUTPUT_ITERATOR>
        void generate_greedy_solution(OUTPUT_ITERATOR begin, OUTPUT_ITERATOR end)
        {
            assert(finalized_graph_);
            assert(end - begin == orig_.end() - orig_.begin());
            assert(assignment_latest_.size() >= orig_.size());
            greedy(); // generate greedy integer solution
            std::copy(assignment_latest_.begin(), assignment_latest_.begin() + orig_.size(), begin);
        }

    protected:
        template <typename F>
        void foreach_clique(F f) const
        {
            for (index clique_idx = 0; clique_idx < no_cliques(); ++clique_idx)
                f(clique_idx);
        }

        template <bool only_orig = false, typename F>
        void foreach_node(F f) const
        {
            for (index node_idx = 0; node_idx < (only_orig ? no_orig() : no_nodes()); ++node_idx)
                f(node_idx);
        }

        template <typename F>
        void foreach_node_in_clique(const index clique_idx, F f) const
        {
            const auto &r = clique_indices_[clique_idx];
            for (index idx = r.begin; idx < r.end; ++idx)
            {
                const index node_idx = clique_index_data_[idx];
                f(node_idx);
            }
        }

        template <typename F>
        void foreach_clique_of_node(const index node_idx, F f) const
        {
            const auto &r = node_cliques_[node_idx];
            for (index idx = r.begin; idx < r.end; ++idx)
            {
                const index clique_idx = node_cliques_data_[idx];
                f(clique_idx);
            }
        }

        template <typename F>
        void foreach_node_neigh(const index node_idx, F f) const
        {
            const auto &r = node_neighs_[node_idx];
            for (index idx = r.begin; idx < r.end; ++idx)
            {
                const index other_node_idx = node_neigh_data_[idx];
                f(other_node_idx);
            }
        }

        cost dual_relaxed_scaled() const
        {
            // If alphas are non-zero this computation fails. Reparametrization would be required before calling this
            // function then.
            assert_unset_alphas();

            // Compute $D(\lambda) = \sum_i \lambda_i + \sum_i \max_{x \in [0, 1]^N} <c^\lambda, x>$.
            // Note that the sum of all lambdas = constant_.
            if constexpr (initial_reparametrization)
            {
                // We know that the second sum is always zero, because this is ensured after the first
                // reparametrization run.
                return constant_;
            }
            else
            {
                // Compute $\sum_i \max_{x \in [0, 1]^N} <c^\lambda, x>$.
                const auto f = [](const auto a, const auto b)
                {
                    return a + std::max(b, cost{0});
                };
                const auto sum_max_i = std::accumulate(costs_.cbegin(), costs_.cend(), cost{0}, f);
                return constant_ + sum_max_i;
            }
        }

        cost dual_smoothed_scaled() const
        {
            // Compute $D^T(\lambda) = \sum_i \lambda_i + \max_{x \in [0, 1]^N} [ <c^\lambda, x> + T H(x) ]$,
            // where $H(x) = - \sum_i (x_i \log x_i - x_i)$, i.e., entropy shifted by linear function.
            // The max of $c_i^\lambda x_i - T (x_i log x_i - x_i)$ is obtained at $x_i = exp(c_i^\lambda / T)$.
            // Hence $D^T(\lambda) = sum_i \lambda_i + T \sum_i exp(c_i / T)$.
            assert_negative_node_costs();

            // sum of all lambdas = constant_
            auto f = [this](const auto a, const auto c)
            { return a + std::exp(c / temperature_); };
            return constant_ + temperature_ * std::accumulate(costs_.cbegin(), costs_.cend(), 0.0, f);
        }

        cost entropy_of_relaxed_assignment_scaled() const
        {
            // Estimate entropy $- \sum_i x_i log x_i - x_i$. We use `assignment_relaxed_` for x.
            auto f = [this](const auto a, const auto x_i)
            {
                if (x_i == 0)
                    return a;
                else
                    return a + x_i * std::log(x_i) - x_i;
            };

            return -std::accumulate(assignment_relaxed_.cbegin(), assignment_relaxed_.cend(), 0.0, f);
        }

        void update_temperature()
        {
            // Note: This is the implementation for section "7.2 Method 2: Duality gap" of the paper.
            const auto d = dual_smoothed_scaled();
            const auto p = value_relaxed_;

            auto new_temp = (d - p) / (entropy_of_relaxed_assignment_scaled() / temperature_drop_factor_);

            assert(std::isnormal(new_temp) && new_temp >= 0.0);

            temperature_ = std::max(std::min(temperature_, new_temp), 1e-10);
        }

        template <typename T>
        bool compute_feasibility(const std::vector<T> &assignment) const
        {
            for (const auto &cl : clique_indices_)
            {
                T sum = 0;
                for (index idx = cl.begin; idx < cl.end; ++idx)
                    sum += assignment[clique_index_data_[idx]];
                if (!feasibility_check(sum))
                    return false;
            }
            return true;
        }

        cost compute_primal(const std::vector<int> &assignment) const
        {
            // Same as relaxed objective, we just check that $x_i \in {0, 1}$.
#ifndef NDEBUG
            for (auto a : assignment)
                assert(a == 0 || a == 1);
#endif

            return compute_primal_relaxed(assignment);
        }

        template <typename T>
        cost compute_primal_relaxed(const std::vector<T> &assignment) const
        {
            // Compute $<c, x>$ s.t. uniqueness constraints.
            assert(assignment.size() == no_nodes());
            cost result = constant_;

            for (index node_idx = 0; node_idx < no_nodes(); ++node_idx)
            {
                const auto x = assignment[node_idx];
                assert(x >= 0 && x <= 1);
                result += costs_[node_idx] * x;
            }

#ifndef NDEBUG
            assert(compute_feasibility(assignment));
#endif

            return result;
        }

        void finalize_graph()
        {
            if (finalized_graph_)
                return;

            // Look for nodes that are not part of any clique. If a node is not part of
            // any clique, there is no lambda that can affect the reduced cost of the
            // node. However, the code assumes that we can shift the reduced cost of
            // any node below zero. To simplify the code we add a fake clique.
            {
                std::vector<bool> tmp(no_nodes(), false);
                for (index clique_idx = 0; clique_idx < no_cliques(); ++clique_idx)
                {
                    const auto &cl = clique_indices_[clique_idx];
                    for (index idx = cl.begin; idx < cl.end; ++idx)
                    {
                        const auto node_idx = clique_index_data_[idx];
                        tmp[node_idx] = true;
                    }
                }

                std::vector<index> clique;
                for (index node_idx = 0; node_idx < no_orig(); ++node_idx)
                {
                    if (!tmp[node_idx])
                    {
                        clique.resize(1);
                        clique[0] = node_idx;
                        add_clique(clique);
                    }
                }
            }

            //
            // Construct node to clique mapping.
            //

            std::vector<std::vector<index>> tmp(no_nodes());
            for (index clique_idx = 0; clique_idx < no_cliques(); ++clique_idx)
            {
                const auto &cl = clique_indices_[clique_idx];
                for (index idx = cl.begin; idx < cl.end; ++idx)
                {
                    const auto node_idx = clique_index_data_[idx];
                    tmp[node_idx].push_back(clique_idx);
                }
            }

            node_cliques_.resize(no_nodes());
            for (index nidx = 0; nidx < no_nodes(); ++nidx)
            {
                auto &current = tmp[nidx];
                auto &nc = node_cliques_[nidx];
                nc.begin = node_cliques_data_.size();
                node_cliques_data_.insert(node_cliques_data_.end(), current.begin(), current.end());
                nc.end = node_cliques_data_.size();
                nc.size = nc.end - nc.begin;
            }

            //
            // Construct node neighborhood.
            //

            for (auto &vec : tmp)
                vec.clear();

            for (const auto &cl : clique_indices_)
            {
                for (index idx0 = cl.begin; idx0 < cl.end; ++idx0)
                    for (index idx1 = cl.begin; idx1 < cl.end; ++idx1)
                        if (idx0 != idx1)
                            tmp[clique_index_data_[idx0]].push_back(clique_index_data_[idx1]);
            }

            node_neighs_.resize(no_nodes());
            for (index nidx = 0; nidx < no_nodes(); ++nidx)
            {
                auto &current = tmp[nidx];
                std::sort(current.begin(), current.end());
                auto end = std::unique(current.begin(), current.end());

                auto &neighbors = node_neighs_[nidx];
                neighbors.begin = node_neigh_data_.size();
                node_neigh_data_.insert(node_neigh_data_.end(), current.begin(), end);
                neighbors.end = node_neigh_data_.size();
                neighbors.size = neighbors.end - neighbors.begin;
            }

            //
            // Construct trivial labeling of zero cost.
            //

            // Set assignment_latest_ to zero labeling.
            assignment_latest_.resize(no_nodes());
            for (index nidx = 0; nidx < no_nodes(); ++nidx)
                assignment_latest_[nidx] = nidx < no_orig() ? 0 : 1;
            value_latest_ = compute_primal(assignment_latest_);
            assert(std::abs(value_latest_) < 1e-8);

            // Set assignment_best_ to zero labeling.
            assignment_best_ = assignment_latest_;
            value_best_ = 0.0;

            // Set assignment_relaxed_ to zero labeling.
            assignment_relaxed_.assign(assignment_latest_.cbegin(), assignment_latest_.cend());
            value_relaxed_ = 0.0;

            // Set assignment_relaxed_best_ to zero labeling.
            assignment_relaxed_best_ = assignment_relaxed_;
            value_relaxed_best_ = 0.0;

            //
            // Initialize remaining things.
            //

            alphas_.resize(no_cliques());

            scratch_greedy_indices_.resize(no_cliques());
            std::iota(scratch_greedy_indices_.begin(), scratch_greedy_indices_.end(), 0);

            scratch_qpbo_indices_.resize(no_orig());

            finalized_graph_ = true;
        }

        void finalize_costs()
        {
            if (finalized_costs_)
                return;

            if constexpr (initial_reparametrization)
            {
                // Update all lambdas (without smoothing, invariants do not hold) to ensure
                // that invariants (negative node costs) hold.
                single_pass<false>();
                auto it = std::min_element(costs_.cbegin(), costs_.cend());
                scaling_ = std::abs(*it);
            }
            else
            {
                auto it = std::max_element(costs_.cbegin(), costs_.cend());
                scaling_ = std::abs(*it);
            }

            constant_ /= scaling_;
            for (auto &c : costs_)
                c /= scaling_;

            // We update the cached values for the corresponding assignments (costs
            // have most likely changed the assignment between calls to
            // finalize_costs).
            value_relaxed_ = compute_primal_relaxed(assignment_relaxed_);
            update_relaxed_best(assignment_relaxed_, value_relaxed_);
            value_best_ = compute_primal(assignment_best_);
            iterations_ = 0;

            // This will set up the corresponding costs in the exponential domain. We
            // do this so that all invariants hold, e.g., for all functions that call
            // `assert_unset_alphas`.
            init_exponential_domain();

            // std::cout << "initial reparametrization: dual=" << dual_relaxed() << " primal=" << primal() << std::endl;

            // Try to improve naive assignment by greedily sampling an assignment.
            // This will be used for inital temperature selection.
            greedy();
            if (value_latest_ > value_best_)
            {
                value_best_ = value_latest_;
                assignment_best_ = assignment_latest_;
            }

            compute_relaxed_truncated_projection();
            update_temperature();

            finalized_costs_ = true;
        }

        template <bool verbose = false>
        cost_exp sparsity_estimation_pass() // BSD
        {
            cost_exp min_max_clique_cost = infinity;

            for (auto cl : clique_indices_)
            {
                copy_clique_in(cl); // copy clique to scratch_

                if (scratch_.cbegin() != scratch_.cend())
                    min_max_clique_cost = std::min(min_max_clique_cost, *std::max_element(scratch_.cbegin(), scratch_.cend()));
            }

            cost_exp threshold = min_max_clique_cost * default_sparsity_threshold;
            std::transform(
                assignment_relaxed_.begin(),
                assignment_relaxed_.end(),
                assignment_relaxed_.begin(),
                [threshold](auto value)
                { return value < threshold ? 0.0 : value; });

            // if (verbose)
            // {
            //     cost_exp num_too_small = std::count(assignment_relaxed_.cbegin(), assignment_relaxed_.cend(), 0.0);

            //     std::cout
            //         << "min-max clique cost="
            //         << min_max_clique_cost
            //         << ", threshold="
            //         << threshold
            //         << "; x-coordinates to be ignored= "
            //         << num_too_small / assignment_relaxed_.size() * 100 << "%"
            //         << std::endl;
            // }

            return threshold;
        }

        template <bool smoothing = true>
        void single_pass()
        {
            for (index clique_idx = 0; clique_idx < no_cliques(); ++clique_idx)
                update_lambda<smoothing>(clique_idx);
        }

        void reduced_single_pass()
        {
            for (index clique_idx = 0; clique_idx < no_cliques(); ++clique_idx)
                reduced_update_lambda(clique_idx);
        }

        void init_exponential_domain()
        {
            // Reset alphas to 1.
            std::fill(alphas_.begin(), alphas_.end(), 1.0);

            // Recompute exponentiated costs.
            assert_negative_node_costs();
            std::transform(costs_.cbegin(), costs_.cend(), assignment_relaxed_.begin(),
                           [this](const auto c)
                           { return std::exp(c / temperature_); });

#ifndef NDEBUG
            for (const auto v : assignment_relaxed_)
                assert(std::isfinite(v));
#endif
        }

        template <bool verbose = true>
        void reduced_init_cliques()
        {
            sparsity_estimation_pass<verbose>();

            reduced_clique_index_data_.resize(0);
            reduced_clique_indices_.resize(0);

            int ave_clique_size = 0;
            int ave_reduced_clique_size = 0;

            for (auto &cl : clique_indices_)
            {
                range reduced_cl;
                reduced_cl.begin = reduced_clique_index_data_.size();
                reduced_cl.size = 0;

                for (index idx = 0; idx < cl.size; ++idx)
                {
                    index current_index = clique_index_data_[cl.begin + idx];
                    if (assignment_relaxed_[current_index] > 0.0)
                    {
                        reduced_clique_index_data_.push_back(current_index);
                        ++reduced_cl.size;
                    }
                }
                reduced_cl.end = reduced_cl.begin + reduced_cl.size;
                reduced_clique_indices_.push_back(reduced_cl);

                ave_clique_size += cl.size;
                ave_reduced_clique_size += reduced_cl.size;

#ifndef NDEBUG
                assert(reduced_cl.size > 0);
#endif
            }

            // computing the average from the sum
            ave_clique_size /= no_cliques();
            ave_reduced_clique_size /= no_cliques();

            // if (verbose)
            // {
            //     std::cout << "average clique size=" << ave_clique_size
            //               << "; average reduced clique size=" << ave_reduced_clique_size << " (" << (float)ave_reduced_clique_size / ave_clique_size * 100 << " %)" << std::endl;
            //     std::cout << std::endl;
            // }
        }

        template <bool smoothing = true>
        void copy_clique_in(const range &cl)
        {
            scratch_.resize(cl.size, 0.0);
            for (index idx = 0; idx < cl.size; ++idx)
                if constexpr (smoothing)
                    scratch_[idx] = assignment_relaxed_[clique_index_data_[cl.begin + idx]];
                else
                    scratch_[idx] = costs_[clique_index_data_[cl.begin + idx]];
        }

        template <bool smoothing = true>
        void copy_clique_out(const range &cl)
        {
            for (index idx = 0; idx < cl.size; ++idx)
                if constexpr (smoothing)
                    assignment_relaxed_[clique_index_data_[cl.begin + idx]] = scratch_[idx];
                else
                    costs_[clique_index_data_[cl.begin + idx]] = scratch_[idx];
        }

        template <bool smoothing = true>
        void update_lambda(const index clique_idx, const std::vector<range> &clique_indices, const std::vector<index> &clique_index_data)
        {
            const auto &cl = clique_indices[clique_idx];
            if (smoothing)
            {
                cost_exp sum = 0.0;

                // sum up x within a clique
                for (index idx = 0; idx < cl.size; ++idx)
                    sum += assignment_relaxed_[clique_index_data[cl.begin + idx]];

                // normalize x within a clique
                for (index idx = 0; idx < cl.size; ++idx)
                    assignment_relaxed_[clique_index_data[cl.begin + idx]] /= sum;

                // update exponentiated duals
                auto &alpha = alphas_[clique_idx];
                alpha /= sum;

                if (alpha + 1 / alpha > default_stabilization_threshold)
                    ++stabilization_counter_;
            }
            else
            {
                cost c_max = -infinity;

                // find max cost within a clique
                for (index idx = 0; idx < cl.size; ++idx)
                    c_max = std::max(c_max, costs_[clique_index_data[cl.begin + idx]]);

                assert(std::isfinite(c_max));
                constant_ += c_max;

                // subtract x_max from all x within a clique
                for (index idx = 0; idx < cl.size; ++idx)
                    costs_[clique_index_data[cl.begin + idx]] -= c_max;
            }
        }

        void reduced_update_lambda(const index clique_idx)
        {
            update_lambda<true>(clique_idx, reduced_clique_indices_, reduced_clique_index_data_);
        };

        template <bool smoothing = true>
        void update_lambda(const index clique_idx)
        {
            update_lambda<smoothing>(clique_idx, clique_indices_, clique_index_data_);
        }

        void reparametrize(const index clique_idx, const cost v)
        {
            assert(clique_idx >= 0 && clique_idx < no_cliques());
            assert(std::isfinite(v));
            assert(std::isfinite(constant_));

            foreach_node_in_clique(clique_idx, [&](const auto node_idx)
                                   {
      auto& cost = costs_[node_idx]; assert(std::isfinite(cost));
      cost -= v; assert(std::isfinite(cost)); });

            constant_ += v;
            assert(std::isfinite(constant_));
        }

        void reparametrize(const index clique_idx)
        {
            assert(clique_idx >= 0 && clique_idx < no_cliques());
            auto &alpha = alphas_[clique_idx];
            assert(std::isfinite(alpha));
            reparametrize(clique_idx, -temperature_ * std::log(alpha));
            alpha = 1.0;
        }

        void reparametrize()
        {
            foreach_clique([this](const auto clique_idx)
                           { reparametrize(clique_idx); });
        }

        template <typename T>
        void update_relaxed_best(const std::vector<T> &assignment, const cost &value)
        {
            assert(assignment_relaxed_best_.size() == assignment.size());
            if (value > value_relaxed_best_)
            {
                assignment_relaxed_best_.assign(assignment.cbegin(), assignment.cend());
                value_relaxed_best_ = value;
            }
        }

        void compute_relaxed_truncated_projection()
        {
            auto node_cost = [this](const index node_idx) -> cost_exp *
            {
                assert(node_idx < no_orig());
                return &assignment_relaxed_[node_idx];
            };

            auto slack = [this](const index clique_idx) -> cost_exp *
            {
                assert(no_orig() + clique_idx < assignment_relaxed_.size());
                return &assignment_relaxed_[no_orig() + clique_idx];
            };

            auto max_allowed = [this, slack](const index node_idx) -> cost_exp
            {
                assert(node_idx < no_orig());
                const auto &nc = node_cliques_[node_idx];
                cost_exp result = 1.0;
                for (index idx = nc.begin; idx < nc.end; ++idx)
                {
                    const auto clique_idx = node_cliques_data_[idx];
                    result = std::min(result, *slack(clique_idx));
                }
                assert(!std::isinf(result));
                return result;
            };

            auto reduce_max_allowed = [this, &slack](const index node_idx, const cost value)
            {
                assert(node_idx < no_orig());
                const auto &nc = node_cliques_[node_idx];
                for (index idx = nc.begin; idx < nc.end; ++idx)
                {
                    const auto clique_idx = node_cliques_data_[idx];
                    assert(value <= *slack(clique_idx));
                    *slack(clique_idx) -= value;
                }
            };

            for (index node_idx = 0; node_idx < no_orig(); ++node_idx)
                *node_cost(node_idx) = std::exp(costs_[node_idx] / temperature_);

            for (index clique_idx = 0; clique_idx < no_cliques(); ++clique_idx)
                *slack(clique_idx) = 1.0;

            foreach_node<true>([&](const auto node_idx)
                               {
      cost_exp* x = node_cost(node_idx);
      *x = std::min(*x, max_allowed(node_idx));
      reduce_max_allowed(node_idx, *x); });

            value_relaxed_ = compute_primal_relaxed(assignment_relaxed_);
            update_relaxed_best(assignment_relaxed_, value_relaxed_);
        }

        bool update_integer_assignment(int greedy_generations)
        {
            bool has_improved = false;
            for (int i = 0; i < greedy_generations; ++i)
                has_improved |= update_integer_assignment();
            // if (greedy_generations > 0)
            //     std::cout << std::endl;
            return has_improved;
        }

        bool update_integer_assignment()
        {
            greedy();
#ifdef ENABLE_QPBO
            const auto old_value_best_ = value_best_;
            fusion_move();
            return value_best_ > old_value_best_;
#else
            if (value_latest_ > value_best_)
            {
                value_best_ = value_latest_;
                assignment_best_ = assignment_latest_;
                return true;
            }
            return false;
#endif
        }

        void greedy()
        {
            // std::cout << "g " << std::flush;
            // std::cout << std::flush;
            std::shuffle(scratch_greedy_indices_.begin(), scratch_greedy_indices_.end(), gen_);

            std::fill(assignment_latest_.begin(), assignment_latest_.end(), -1);
            for (const auto clique_idx : scratch_greedy_indices_)
                greedy_clique(clique_idx);

            value_latest_ = compute_primal(assignment_latest_);
            update_relaxed_best(assignment_latest_, value_latest_);
        }

        void greedy_clique(const index clique_idx)
        {
            assert(clique_idx >= 0 && clique_idx < no_cliques());

            int count = 0;
            foreach_node_in_clique(clique_idx, [&](const auto node_idx)
                                   { count += assignment_latest_[node_idx] == 1 ? 1 : 0; });
            assert(count == 0 || count == 1);

            if (count > 0)
                return;

            cost max = -infinity;
            index argmax = -1;
            foreach_node_in_clique(clique_idx, [&](const auto node_idx)
                                   {
      if (assignment_latest_[node_idx] != 0 && costs_[node_idx] > max) {
        max = costs_[node_idx];
        argmax = node_idx;
      } });

            assignment_latest_[argmax] = 1;
            foreach_node_neigh(argmax, [&](const auto node_idx)
                               {
      assert(node_idx != argmax);
      assert(assignment_latest_[node_idx] != 1);
      assignment_latest_[node_idx] = 0; });
        }

#ifdef ENABLE_QPBO
        bool fuse_two_assignments(std::vector<int> &a0, const std::vector<int> &a1)
        {
            constexpr cost QPBO_INF = 1e20;
            assert(a0.size() == no_nodes());
            assert(a1.size() == no_nodes());

            auto reset_qpbo_indices = [&]()
            {
                std::fill(scratch_qpbo_indices_.begin(), scratch_qpbo_indices_.end(), -1);
            };

            auto is_node_present = [&](auto nidx)
            {
                return nidx < no_orig() && scratch_qpbo_indices_[nidx] != -1;
            };

            auto for_each_present_node_tuple = [&](auto func)
            {
                for (index nidx1 = 0; nidx1 < no_orig(); ++nidx1)
                {
                    if (is_node_present(nidx1))
                    {
                        for (index idx = node_neighs_[nidx1].begin; idx < node_neighs_[nidx1].end; ++idx)
                        {
                            const auto nidx2 = node_neigh_data_[idx];
                            if (nidx1 < nidx2 && is_node_present(nidx2))
                            {
                                func(nidx1, nidx2);
                            }
                        }
                    }
                }
            };

            auto enable_all_dumies = [&]()
            {
                std::fill(a0.begin() + no_orig(), a0.end(), 1);
            };

            auto disable_dummy_of_clique = [&](const auto clique_idx)
            {
                assert(clique_idx >= 0 && clique_idx < clique_indices_.size());
                assert(clique_idx >= 0 && clique_idx < clique_indices_.size());
                const auto &cl = clique_indices_[clique_idx];
                assert(cl.size >= 2);
                assert(cl.end - 1 >= 0 && cl.end - 1 < clique_index_data_.size());
                const auto nidx = clique_index_data_[cl.end - 1];
                assert(nidx >= no_orig());
                assert(nidx < no_nodes());
                assert(a0[nidx] == 1);
                a0[nidx] = 0;
            };

            auto disable_dummy_for_node = [&](const auto nidx)
            {
                assert(nidx >= 0 && nidx < no_orig());
                assert(a0[nidx] == 1);
                const auto &nc = node_cliques_[nidx];
                for (auto idx = nc.begin; idx < nc.end; ++idx)
                {
                    assert(idx >= 0 && idx < node_cliques_data_.size());
                    const auto clique_idx = node_cliques_data_[idx];
                    disable_dummy_of_clique(clique_idx);
                }
            };

            qpbo_.Reset();
            reset_qpbo_indices();

            for (index nidx = 0; nidx < no_orig(); ++nidx)
            {
                assert(a0[nidx] == 0 || a0[nidx] == 1);
                assert(a1[nidx] == 0 || a1[nidx] == 1);

                if (a0[nidx] != a1[nidx])
                {
                    const bool l0 = a0[nidx] == 1, l1 = a1[nidx] == 1;
                    cost c0 = l0 ? -orig_[nidx] : 0.0;
                    cost c1 = l1 ? -orig_[nidx] : 0.0;
                    const auto qpbo_index = qpbo_.AddNode();
                    qpbo_.AddUnaryTerm(qpbo_index, c0, c1);
                    scratch_qpbo_indices_[nidx] = qpbo_index;
                }
            }

#ifndef NDEBUG
            const auto qpbo_size = qpbo_.GetNodeNum();
            // std::cout << "[DBG] qpbo_size = " << qpbo_size << " / " << no_orig() << " ("
            //           << (100.0f * qpbo_size / no_orig()) << "%)" << std::endl;
#endif

            for_each_present_node_tuple([&](const auto nidx1, const auto nidx2)
                                        {
      const cost c01 = a0[nidx1] == 1 && a1[nidx2] == 1 ? QPBO_INF : 0.0;
      const cost c10 = a1[nidx1] == 1 && a0[nidx2] == 1 ? QPBO_INF : 0.0;

#ifndef NDEBUG
      const cost c00 = a0[nidx1] == 1 && a0[nidx2] == 1 ? QPBO_INF : 0.0;
      const cost c11 = a1[nidx1] == 1 && a1[nidx2] == 1 ? QPBO_INF : 0.0;
      assert(std::abs(c00) < 1e-8);
      assert(std::abs(c11) < 1e-8);
#else
      const cost c00 = 0.0, c11 = 0.0;
#endif

      if (c00 || c01 || c10 || c11) {
        const auto qpbo_idx1 = scratch_qpbo_indices_[nidx1];
        const auto qpbo_idx2 = scratch_qpbo_indices_[nidx2];
        qpbo_.AddPairwiseTerm(qpbo_idx1, qpbo_idx2, c00, c01, c10, c11);
      } });

            qpbo_.Solve();

            bool changed = false;
            enable_all_dumies();
            for (index nidx = 0; nidx < no_orig(); ++nidx)
            {
                const auto qpbo_idx = scratch_qpbo_indices_[nidx];
                const auto l = qpbo_idx != -1 ? qpbo_.GetLabel(qpbo_idx) : 0;
                changed = changed || (l != 0);
                a0[nidx] = l == 0 ? a0[nidx] : a1[nidx];

                if (a0[nidx] == 1)
                    disable_dummy_for_node(nidx);
            }

            assert(compute_feasibility(a0));

            return changed;
        }

        void fusion_move()
        {
            // std::cout << "f " << std::flush;
#ifndef NDEBUG
            const auto value_best_old = value_best_;
#endif
            if (fuse_two_assignments(assignment_best_, assignment_latest_))
                value_best_ = compute_primal(assignment_best_);
            assert(dbg::are_identical(value_best_, compute_primal(assignment_best_)));
            assert(value_best_ >= value_best_old - 1e-8);
            update_relaxed_best(assignment_best_, value_best_);
        }
#endif

        void assert_unset_alphas() const
        {
#ifndef NDEBUG
            for (const auto v : alphas_)
                assert(std::abs(v - 1) < epsilon);
#endif
        }

        void assert_negative_node_costs() const
        {
#ifndef NDEBUG
            for (const auto c : costs_)
                assert(c <= 0);
#endif
        }

        bool finalized_graph_, finalized_costs_;
        std::vector<cost> costs_;
        std::vector<cost> orig_;
        cost constant_;
        double scaling_;

        int iterations_;
        double temperature_;

        std::vector<range> clique_indices_;
        std::vector<index> clique_index_data_;

        std::vector<range> reduced_clique_indices_;
        std::vector<index> reduced_clique_index_data_;
        int stabilization_counter_;

        std::vector<range> node_cliques_;
        std::vector<index> node_cliques_data_;

        std::vector<range> node_neighs_;
        std::vector<index> node_neigh_data_;

        mutable std::vector<cost_exp> scratch_;
        mutable std::vector<index> scratch_greedy_indices_;
        mutable std::vector<index> scratch_qpbo_indices_;

        cost value_latest_, value_best_;
        std::vector<int> assignment_latest_, assignment_best_;

        cost value_relaxed_, value_relaxed_best_;
        std::vector<cost_exp> assignment_relaxed_, assignment_relaxed_best_;
        std::vector<cost_exp> alphas_;

        std::default_random_engine gen_;
#ifdef ENABLE_QPBO
        qpbo::QPBO<cost> qpbo_;
#endif

        cost temperature_drop_factor_;
    };

} // namespace mpopt::mwis

#endif

/* vim: set ts=8 sts=2 sw=2 et ft=cpp: */
