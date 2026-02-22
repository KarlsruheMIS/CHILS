#pragma once

/*
    API for CHILS

    You can either construct the graph using the add_vertex
    and add_edge functions or provide a graph directly on
    the CSR format. After construction, there are two
    functions to run either CHILS with N solutions or just
    the local search (recommended for low time limits). To
    get the independent set, use the functions prefixed by
    chils_solution_.

    Note that you need to call chils_free to free the
    memory used by the heuristic. Also, note that only one
    form of graph construction should be used. That is,
    either add edges and vertices or set the CSR graph
    directly. You can not give a CSR graph and then add
    more vertices or edges.
*/

#ifdef __cplusplus
extern "C"
{
#endif

    /**
     * @brief Initialize the heuristic.
     *
     * @returns A void pointer that should be used as the input
     * parameter 'solver' for the rest of the API.
     */
    void *chils_initialize();

    /**
     * @brief Free all allocated memory.
     *
     * @param solver pointer to the solver object.
     */
    void chils_release(void *solver);

    /**
     * @brief Add a new vertex with the given weight.
     *
     * @param solver pointer to the solver object.
     * @param weight the weight of the new vertex.
     *
     * @returns The ID of the new vertex. This will always be
     * one larger than the previous time the function was
     * called, starting with 0.
     */
    int chils_add_vertex(void *solver, long long weight);

    /**
     * @brief Add a new edge to the graph. The edge is considered
     * as undirected, meaning you do not need to add both {u,v}
     * and {v,u}.
     *
     * @attention It is not required that both u and v have
     * already been added to the graph at this point. That
     * requirement only applies when the heuristic starts.
     *
     * @param solver pointer to the solver object.
     * @param u the first endpoint of the edge.
     * @param v the second endpoint of the edge.
     */
    void chils_add_edge(void *solver, int u, int v);

    /**
     * @brief Set the input graph for the heuristic on the compressed
     * sparse row (CSR) format. The number of vertices can be at most
     * 2^31, but the number of edges may exceed this limit.
     *
     * @attention This function makes a copy of the graph.
     *
     * @example A cycle with four vertices, the input could look like
     * this: n = 4, xadj = {0,2,4,6,8}, adjncy = {1,2,0,3,0,3,1,2}.
     *
     * @param solver pointer to the solver object.
     * @param n the number of vertices in the graph.
     * @param xadj array of size n + 1 that stores the interval where the
     * neighborhood of each vertex is stored in adjncy.
     * @param adjncy compressed array of neighborhoods, stored one after the
     * other in one continuous array. The neighborhood of vertex i is placed
     * from adjncy[xadj[i]] until (not including) adjncy[xadj[i + 1]].
     * @param weights array of size n holding the vertex weights.
     */
    void chils_set_graph(void *solver, int n, const long long *xadj, const int *adjncy, const long long *weights);

    /**
     * @brief Run the heuristic for a certain number of seconds.
     *
     * @attention The alternating time between full graph and the
     * d-core is 10 seconds. For very short time limits, use the
     * local_search_only function instead. Also, repeated calls
     * to this function ignore the new values for n_solutions and
     * seed. It simply continues for time_limit seconds starting
     * from where the last call ended.
     *
     * @param solver pointer to the solver object.
     * @param time_limit time limit in seconds.
     * @param n_solutions number of solutions to use for CHILS.
     * @param seed seed for the random number generator
     */
    void chils_run_full(void *solver, double time_limit, int n_solutions, unsigned int seed);

    /**
     * @brief Run the heuristic using only the baseline local search.
     *
     * @attention Repeated calls to this function ignore the new seed.
     * It simply continues for time_limit seconds starting from where
     * the last call ended.
     *
     * @param solver pointer to the solver object.
     * @param time_limit time limit in seconds.
     * @param seed seed for the random number generator
     */
    void chils_run_local_search_only(void *solver, double time_limit, unsigned int seed);

    /**
     * @brief Extract the number of vertices in the best independent
     * set found by the heuristic.
     *
     * @param solver pointer to the solver object.
     *
     * @return The size of the best solution found.
     */
    int chils_solution_get_size(void *solver);

    /**
     * @brief Extract the weight of the best independent set found
     * by the heuristic.
     *
     * @param solver pointer to the solver object.
     *
     * @return The weight of the best solution found.
     */
    long long chils_solution_get_weight(void *solver);

    /**
     * @brief Extract the time when the heuristic found the best
     * independent set. Measured from the last run call.
     *
     * @param solver pointer to the solver object.
     *
     * @return The time in seconds it took to find the best solution.
     */
    double chils_solution_get_time(void *solver);

    /**
     * @brief Extract the configuration of a particular vertex
     * in the best independent set found.
     *
     * @param solver pointer to the solver object.
     * @param u the ID of the vertex in question.
     * @return 1 if the vertex is in the independent set, 0 otherwise.
     */
    int chils_solution_get_vertex_configuration(void *solver, int u);

    /**
     * @brief Extract a list of the vertices in the best
     * independent set found.
     *
     * @attention To know the number of vertices in the list, use the
     * function chils_solution_get_size. Note that this array must
     * be freed by the programmer.
     *
     * @param solver pointer to the solver object.
     * @return A pointer to an integer array with the IDs of the
     * vertices making up the best independent set.
     */
    int *chils_solution_get_independent_set(void *solver);

    /**
     * @brief Requests the heuristic to stop gracefully.
     * This function can be called from a signal handler or another thread to interrupt
     * the process.
     */
    void chils_request_stop();

#ifdef __cplusplus
}
#endif