#pragma once
#include <stdio.h>
#include <signal.h>

extern volatile sig_atomic_t keep_running;

typedef struct
{
    int n;        // Number of vertices
    long long m;  // Number of edges
    long long *V; // Neighborhood pointers
    int *E;       // Edgelist
    long long *W; // Vertex weights
} graph;

graph *graph_parse(FILE *f);

void graph_store(FILE *f, graph *g);

void graph_free(graph *g);

int graph_validate(graph *g);

graph *graph_subgraph(graph *g, int *Mask, int *RM);

// Should be called inside parallel region
void graph_subgraph_par(graph *g, graph *sg, int *Mask, int *RM, int *FM, long long *S1, long long *S2);