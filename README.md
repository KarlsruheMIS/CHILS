# CHILS &mdash; Concurrent Hybrid Iterated Local Search
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/829941140.svg)](https://doi.org/10.5281/zenodo.15173363)
[![Build and Test](https://github.com/KennethLangedal/CHILS/actions/workflows/github-actions.yml/badge.svg)](https://github.com/KennethLangedal/CHILS/actions/workflows/github-actions.yml)

![CHILS illustration](CHILS.png)

Part of the [KaMIS &mdash; Karlsruhe Maximum Independent Sets](https://github.com/KarlsruheMIS) open source framework to compute (weighted) independent sets and vertex covers.

## Installation

On linux with GCC, compile the project by running
```
make
```
This will produce an executable called CHILS and a library libCHILS.a.

There are no dependencies besides an OpenMP-compatible compiler. Adjust the CC variable in the Makefile if you prefer another compiler than GCC. For other common compilers and operating systems, see the later [installation](#installation-for-macos-and-windows) section. For details on how to use CHILS as a library in your project, see the [API](#api) section.

## Program Options

Note that the `-p N` option is necessary to run CHILS (using **N** concurrent solutions). Otherwise, the program defaults to baseline local search.

| Option | Description | Default | Mandatory
|-|-|-|-
| `-h` | Display help information | | 
| `-v` | Verbose mode, shows continuous updates to STDOUT | |
| `-b` | Blocked mode, output additional results after 10%, 50%, and 100% of the  max time/iterations | |
| `-g path` | Path to the input graph, see input forma | | &check;
| `-i path` | Path to the initial solution, see output format | |
| `-f path` | Path to folder containing initial solutions | |
| `-o path` | Path to store the best solution found, see output format | Not stored |
| `-p N` | Run CHILS with **N** concurrent solutions, use 1 for baseline local seach | 16 |
| `-t sec` | Timeout in seconds | 3600 (1h) |
| `-s sec` | Alternating interval for CHILS in seconds | 10 |
| `-q N` | Max queue size **N** after perturbe | 32 |
| `-c T` | Set number of threads to **T** | OMP_NUM_THREADS |
| `-r s` | Set random seed to **s** | time(NULL) |
| `-n it` | Set max CHILS iterations to **it** | inf |
| `-m it` | Set max local search iterations to **it** | inf |

The output of the program without the `-v` option is a single line on the format
```
instance_name,#vertices,#edges,WIS,time
```
For example, after running
```
./CHILS -g CR-T-D-4.graph
```
The output could look something like this
```
CR-T-D-4,651861,220480534,4922752,3557.2678
```
Where 4922752 is the solution size after 1 hour. 3557.2678 is the time when the best solution was found.

## How to Use

Examples of typical use cases are listed below. Naturally, change `-t` (time limit), `-q` (max queue size), and `-s` (alternating CHILS interval) as necessary.

### Baseline Local Search

```
./CHILS -g [path] -p 1
```

### Sequential CHILS

```
./CHILS -g [path] -c 1
```

### Parallel CHILS

For best parallel performance, use
```
export OMP_NUM_THREADS=N
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
```

And then run

```
./CHILS -g [path]
```

## Python Interface

You can also interact directly with CHILS through its Python API. Install CHILS for Python using

```
pip install CHILS
```

The following example illustrates how to run CHILS from Python. It creates a graph with three vertices and then runs the full CHILS algorithm using four concurrent solutions for 20 seconds.

```python
import chils

solver = chils.CHILS()

solver.add_vertex(10)
solver.add_vertex(10)
solver.add_vertex(20)

solver.add_edge(0, 1)
solver.add_edge(0, 2)

solver.run_full(20.0, 4, 0)

print(solver.get_solution_weight())
```

For more information, see the Python interface [repository](https://github.com/KennethLangedal/CHILS_Python), or the pypi [website](https://pypi.org/project/CHILS/).

## Input Format

CHILS expects graphs on the METIS graph format. A graph with **N** vertices is stored using **N + 1** lines. The first line lists the number of vertices, the number of edges, and the weight type. For CHILS, the first line should use 10 as the weight type to indicate integer vertex weights. Each subsequent line first gives the weight and then lists the neighbors of that node.

Here is an example of a graph with 3 vertices of weight 15, 15, and 20, where the weight 20 vertex is connected to the two 15-weight vertices.

```
3 2 10
15 3
15 3
20 1 2
```

Notice that vertices are 1-indexed, and edges appear in the neighborhoods of both endpoints.

## Output Format

The output format, also used as input format for the `-i` option, is simply a list of 1-indexed vertices, one vertex per line. For example, the solution to the graph used above would look like this

```
1
2
```

## API

The API for the library version of the project is defined in the [chils.h](include/chils.h) header. You only need the header file and the libCHILS.a to use CHILS as a library. The header file contains detailed information about the library's functionality. For example, consider the following C++ program constructing the same graph as mentioned above.

```c++
#include <iostream>

#include "chils.h"

int main(int argc, char **argv)
{
    void *solver = chils_initialize();

    // Creating a graph: 15---20---15

    chils_add_vertex(solver, 15); // Vertex 0
    chils_add_vertex(solver, 15); // Vertex 1
    chils_add_vertex(solver, 20); // Vertex 2

    chils_add_edge(solver, 0, 2); // Edge {0, 2}
    chils_add_edge(solver, 2, 1); // Edge {2, 1}

    // Local search is recommended for small time limits (< 5min)
    chils_run_local_search_only(solver, 1.0, 0);

    // CHILS is recommended for larger time limits (> 5min)
    chils_run_full(solver, 1.0, 8, 0);

    std::cout << chils_solution_get_weight(solver) << std::endl;

    chils_release(solver);

    return 0;
}
```

Assuming chils.h and libCHILS.a resides in the same directory as the main.cpp, you can compile this program using the following command.

```
g++ -fopenmp main.cpp -o prog -L. -lCHILS
```

## Installation for macOS and Windows

Here are some guidlines for use with other operating systems and compilers than GCC on linux. If you have issues compiling the code, it might help to check the [.github/workflows/github-actions.yml](.github/workflows/github-actions.yml) file for details on how we tested the code for these platforms. 

### Linux clang

```
sudo apt-get update
sudo apt-get install -y libomp-dev
```
```
make CC=clang
```

### macOS GCC

```
brew install gcc
```
```
make CC=$(find $(brew --prefix gcc)/bin -name "gcc-*" | head -n 1)
```
(You can also select a specific version if you prefer)

### macOS clang

```
brew install llvm libomp
```
```
make CC=$(brew --prefix llvm)/bin/clang CFLAGS="-Xpreprocessor -I$(brew --prefix libomp)/include"
```

### Windows GCC

```
choco install make
choco install mingw
```
```
make
```

## Reproducing Results from the Paper

An archive version of the paper is available on arXiv combined with the [Learn and Reduce](https://github.com/ernestine-grossmann/MWIS_learn_and_reduce) reduction framework. The CHILS part of this combined paper was accepted at SEA 2025 and will be published in the conference proceedings later in 2025.
```
@article{grossmann2024accelerating,
  title    = {Accelerating Reductions Using Graph Neural Networks and a New Concurrent Local Search for the Maximum Weight Independent Set Problem},
  author   = {Gro{\ss}mann, Ernestine and Langedal, Kenneth and Schulz, Christian},
  journal  = {arXiv preprint arXiv:2412.14198},
  year     = {2024}
}
```

Using the `-b` option gives additional information about solution size after 10%, 50%, and 100% of the time/iterations. The single output line then looks like
```
instance_name,#vertices,#edges,WIS_after_10%,WIS_after_50%,WIS_after_100%,time
```
For example, after running
```
./CHILS -g CR-T-D-4.graph -b
```
The output could look something like this
```
CR-T-D-4,651861,220480534,4900379,4919890,4922752,3557.2678
```
Where 4900379 is the solution size after 10 minutes, 4919890 after 30, and 4922752 after 1 hour. 3557.2678 is the time when the best solution was found.

### Parallel Experiments

The options `-n` and `-m` are used to set an exact number of CHILS/local search iterations. Together with `-r` to set a specific random seed, the execution can be made deterministic for parallel scalability experiments.
