#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>

typedef struct
{
    int n;        // Number of vertices
    long long m;  // Number of edges
    long long *V; // Neighborhood pointers
    int *E;       // Edgelist
    long long *W; // Vertex weights
} graph;

graph *graph_parse(FILE *f);

void graph_free(graph *g);

int intersection(const int *A, int n, const int *B, int m, int *C)
{
    int i = 0, j = 0, k = 0;
    while (i < n && j < m)
    {
        if (A[i] < B[j])
            i++;
        else if (A[i] > B[j])
            j++;
        else
        {
            C[k] = A[i];
            i++;
            j++;
            k++;
        }
    }
    return k;
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: ./METIS_TO_CLIQUE [input graph] [output cliques]\n");
        return 1;
    }

    FILE *f = fopen(argv[1], "r");
    graph *g = graph_parse(f);
    fclose(f);

    int *E_mark = malloc(sizeof(int) * g->m / 2);
    for (long long i = 0; i < g->m / 2; i++)
    {
        E_mark[i] = 1;
    }

    long long *E_index = malloc(sizeof(long long) * g->m);
    long long *E_start = malloc(sizeof(long long) * g->n);
    long long ei = 0;

    int *Temp1 = malloc(sizeof(int) * g->n);
    int *Temp2 = malloc(sizeof(int) * g->n);

    for (int u = 0; u < g->n; u++)
    {
        int first = 1;
        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            long long p;
            if (u > v)
            {
                p = E_index[E_start[v]++];
            }
            else
            {
                p = ei++;
                if (first)
                {
                    E_start[u] = i;
                    first = 0;
                }
            }
            E_index[i] = p;
        }
    }

    int *T = malloc(sizeof(int) * g->m);
    long long *TV = malloc(sizeof(long long) * g->m);

    int *M = malloc(sizeof(int) * g->n);
    for (int i = 0; i < g->n; i++)
        M[i] = 0;

    long long nc = 0;
    long long ncc = 0;

    double t0 = omp_get_wtime();

    TV[0] = 0;

    for (int u = 0; u < g->n; u++)
    {
        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
        {
            if (!E_mark[E_index[i]])
                continue;

            int v = g->E[i];

            int cn = ncc;
            T[cn++] = u;
            T[cn++] = v;

            M[u] = 1;
            M[v] = 1;

            E_mark[E_index[i]] = 0;

            int k = intersection(g->E + g->V[u], g->V[u + 1] - g->V[u], g->E + g->V[v], g->V[v + 1] - g->V[v], Temp1);
            while (k > 0)
            {
                int w = Temp1[0];
                for (long long j = g->V[w]; j < g->V[w + 1]; j++)
                {
                    int z = g->E[j];
                    if (M[z])
                        E_mark[E_index[j]] = 0;
                }
                T[cn++] = w;
                M[w] = 1;

                k = intersection(g->E + g->V[w], g->V[w + 1] - g->V[w], Temp1, k, Temp2);
                int *t = Temp1;
                Temp1 = Temp2;
                Temp2 = t;
            }

            for (int j = ncc; j < cn; j++)
            {
                M[T[j]] = 0;
            }

            TV[++nc] = cn;
            ncc = cn;
        }
    }

    double elapsed = omp_get_wtime() - t0;

    printf("%lld,%lld,%.4lf\n", nc, ncc, elapsed);

    f = fopen(argv[2], "w");

    fprintf(f, "{\n\t\"nodes\": [\n");
    for (int u = 0; u < g->n; u++)
    {
        if (u == g->n - 1)
            fprintf(f, "\t\t%lld\n", g->W[u]);
        else
            fprintf(f, "\t\t%lld,\n", g->W[u]);
    }
    fprintf(f, "\t],\n\t\"cliques\": [\n");

    int first_print = 1;

    for (int i = 0; i < nc; i++)
    {
        if (!first_print)
            fprintf(f, "\t\t],\n");
        first_print = 0;

        fprintf(f, "\t\t[\n");

        for (long long j = TV[i]; j < TV[i + 1]; j++)
        {
            if (j == TV[i + 1] - 1)
                fprintf(f, "\t\t\t%d\n", T[j]);
            else
                fprintf(f, "\t\t\t%d,\n", T[j]);
        }
    }

    fprintf(f, "\t\t]\n\t]\n}\n");
    fclose(f);

    graph_free(g);

    free(E_mark);
    free(E_index);
    free(E_start);

    free(Temp1);
    free(Temp2);

    free(T);
    free(TV);
    free(M);

    return 0;
}

static inline void parse_id(char *Data, size_t *p, long long *v)
{
    while (Data[*p] < '0' || Data[*p] > '9')
        (*p)++;

    *v = 0;
    while (Data[*p] >= '0' && Data[*p] <= '9')
        *v = (*v) * 10 + Data[(*p)++] - '0';
}

graph *graph_parse(FILE *f)
{
    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *Data = malloc(size);
    size_t red = fread(Data, 1, size, f);
    size_t p = 0;

    long long n, m, t;
    parse_id(Data, &p, &n);
    parse_id(Data, &p, &m);
    parse_id(Data, &p, &t);

    if (n >= INT_MAX)
    {
        fprintf(stderr, "Number of vertices must be less than %d, got %lld\n", INT_MAX, n);
        exit(1);
    }

    long long *V = malloc(sizeof(long long) * (n + 1));
    int *E = malloc(sizeof(int) * (m * 2));

    long long *W = malloc(sizeof(long long) * n);

    long long ei = 0;
    for (int u = 0; u < n; u++)
    {
        parse_id(Data, &p, W + u);
        V[u] = ei;
        while (ei < m * 2)
        {
            while (Data[p] == ' ')
                p++;
            if (Data[p] == '\n')
                break;

            long long e;
            parse_id(Data, &p, &e);

            if (e > n)
            {
                fprintf(stderr, "Edge endpoint out of bounds, {%lld, %lld}\n", u + 1ll, e);
                exit(1);
            }

            E[ei++] = e - 1;
        }
        p++;
    }
    V[n] = ei;

    free(Data);

    graph *g = malloc(sizeof(graph));
    *g = (graph){.n = n, .m = ei, .V = V, .E = E, .W = W};

    return g;
}

void graph_free(graph *g)
{
    if (g == NULL)
        return;

    free(g->V);
    free(g->E);
    free(g->W);

    free(g);
}