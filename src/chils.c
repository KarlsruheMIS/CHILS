#include "chils.h"

#include "graph.h"
#include "chils_internal.h"
#include "local_search.h"

#include <stdlib.h>
#include <limits.h>

volatile sig_atomic_t keep_running = 1;

typedef struct
{
    int solver_mode;

    int n;
    long long _an, m, _am;
    long long *W;
    int *X, *Y;

    graph *g;
    chils *c;
    local_search *ls;

    long long cost;
    int size;
    double time;
    int *IS;
} api_data;

void *chils_initialize()
{
    api_data *d = malloc(sizeof(api_data));

    d->solver_mode = 0;

    d->n = 0;
    d->m = 0;
    d->_an = (1 << 10);
    d->_am = (1 << 10);

    d->W = malloc(sizeof(long long) * d->_an);
    d->X = malloc(sizeof(int) * d->_am);
    d->Y = malloc(sizeof(int) * d->_am);

    d->g = NULL;
    d->c = NULL;
    d->ls = NULL;

    d->cost = 0;
    d->size = 0;
    d->time = 0.0;
    d->IS = NULL;

    return d;
}

void chils_release(void *solver)
{
    api_data *d = (api_data *)solver;

    free(d->W);
    free(d->X);
    free(d->Y);

    graph_free(d->g);
    chils_free(d->c);
    local_search_free(d->ls);

    free(d);
}

int chils_add_vertex(void *solver, long long weight)
{
    api_data *d = (api_data *)solver;

    if (d->n >= INT_MAX - 1)
    {
        fprintf(stderr, "Number of vertices must be less than %d\n", INT_MAX);
        exit(1);
    }

    if (d->n == d->_an)
    {
        d->_an *= 2;
        d->W = realloc(d->W, sizeof(long long) * d->_an);
    }

    d->W[d->n] = weight;
    d->n++;

    return d->n - 1;
}

void chils_add_edge(void *solver, int u, int v)
{
    api_data *d = (api_data *)solver;

    if (d->m == d->_am)
    {
        d->_am *= 2;
        d->X = realloc(d->X, sizeof(int) * d->_am);
        d->Y = realloc(d->Y, sizeof(int) * d->_am);
    }

    d->X[d->m] = u;
    d->Y[d->m] = v;
    d->m++;
}

void chils_set_graph(void *solver, int n, const long long *xadj, const int *adjncy, const long long *weights)
{
    api_data *d = (api_data *)solver;

    d->g = malloc(sizeof(graph));
    d->g->n = n;
    d->g->m = xadj[n];
    d->g->V = malloc(sizeof(long long) * (n + 1));
    d->g->E = malloc(sizeof(int) * xadj[n]);
    d->g->W = malloc(sizeof(long long) * n);

    for (int u = 0; u <= n; u++)
        d->g->V[u] = xadj[u];

    for (long long i = 0; i < xadj[n]; i++)
        d->g->E[i] = adjncy[i];

    for (int u = 0; u < n; u++)
        d->g->W[u] = weights[u];
}

static inline int chils_api_compare(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

void chils_construct_graph(api_data *d)
{
    long long *V = malloc(sizeof(long long) * (d->n + 1));
    int *E = malloc(sizeof(int) * d->m * 2);
    long long *W = malloc(sizeof(long long) * d->n);

    for (int i = 0; i < d->n; i++)
        W[i] = d->W[i];

    for (int i = 0; i <= d->n; i++)
        V[i] = 0;

    for (int i = 0; i < d->m; i++)
    {
        if (d->X[i] < 0 || d->X[i] >= d->n || d->Y[i] < 0 || d->Y[i] >= d->n)
        {
            fprintf(stderr, "Invalid edge {%d,%d} in graph with %d vertices\n", d->X[i], d->Y[i], d->n);
            continue;
        }
        V[d->X[i] + 1]++;
        V[d->Y[i] + 1]++;
    }

    int ps = 0;
    for (int i = 1; i <= d->n; i++)
    {
        ps += V[i];
        V[i] = ps - V[i];
    }

    for (long long i = 0; i < d->m; i++)
    {
        if (d->X[i] < 0 || d->X[i] >= d->n || d->Y[i] < 0 || d->Y[i] >= d->n)
        {
            continue;
        }
        E[V[d->X[i] + 1]++] = d->Y[i];
        E[V[d->Y[i] + 1]++] = d->X[i];
    }

    for (int i = 0; i < d->n; i++)
        qsort(E + V[i], V[i + 1] - V[i], sizeof(int), chils_api_compare);

    long long p = 0;
    for (int u = 0; u < d->n; u++)
    {
        int s = p;
        for (long long i = V[u]; i < V[u + 1]; i++)
        {
            int v = E[i];
            if (u == v || (i > V[u] && v <= E[i - 1]))
                continue;

            E[p++] = v;
        }
        V[u] = s;
    }
    V[d->n] = p;

    E = realloc(E, sizeof(int) * p);

    graph *g = malloc(sizeof(graph));
    *g = (graph){.n = d->n, .m = p, .V = V, .E = E, .W = W};

    d->g = g;
}

void chils_run_full(void *solver, double time_limit, int n_solutions, unsigned int seed)
{
    api_data *d = (api_data *)solver;

    if (d->g == NULL)
        chils_construct_graph(d);

    if (!graph_validate(d->g))
    {
        fprintf(stderr, "Detected errors in the graph, make sure "
                        "there are no self edges or missing endpoints\n");
        return;
    }

    if (d->c == NULL)
        d->c = chils_init(d->g, n_solutions, seed);

    chils_run(d->g, d->c, time_limit, LLONG_MAX, 0);

    if (d->cost < d->c->cost || (d->cost == d->c->cost && d->time > d->c->time))
    {
        d->cost = d->c->cost;
        d->size = d->c->size;
        d->time = d->c->time;
        d->IS = chils_get_best_independent_set(d->c);
    }
}

void chils_run_local_search_only(void *solver, double time_limit, unsigned int seed)
{
    api_data *d = (api_data *)solver;

    if (d->g == NULL)
        chils_construct_graph(d);

    if (!graph_validate(d->g))
    {
        fprintf(stderr, "Detected errors in the graph, make sure "
                        "there are no self edges or missing endpoints\n");
        return;
    }

    if (d->ls == NULL)
        d->ls = local_search_init(d->g, seed);

    local_search_explore(d->g, d->ls, time_limit, LLONG_MAX, 0);

    if (d->cost < d->ls->cost || (d->cost == d->ls->cost && d->time > d->ls->time))
    {
        d->cost = d->ls->cost;
        d->size = d->ls->size;
        d->time = d->ls->time;
        d->IS = d->ls->independent_set;
    }
}

int chils_solution_get_size(void *solver)
{
    api_data *d = (api_data *)solver;

    return d->size;
}

long long chils_solution_get_weight(void *solver)
{
    api_data *d = (api_data *)solver;

    return d->cost;
}

double chils_solution_get_time(void *solver)
{
    api_data *d = (api_data *)solver;

    return d->time;
}

int chils_solution_get_vertex_configuration(void *solver, int u)
{
    api_data *d = (api_data *)solver;

    if (d->IS == NULL || u < 0 || u >= d->g->n)
        return 0;

    return d->IS[u];
}

int *chils_solution_get_independent_set(void *solver)
{
    api_data *d = (api_data *)solver;

    if (d->IS == NULL)
        return NULL;

    int *res = malloc(sizeof(int) * d->size);

    int p = 0;
    for (int i = 0; i < d->g->n; i++)
        if (d->IS[i])
            res[p++] = i;

    return res;
}

void chils_request_stop()
{
    keep_running = 0;
}