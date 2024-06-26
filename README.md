# Implementation of the Reverse Shortest Path Problem in Unit-Disk Graps

Given a set $P$ of $n$ points in the plane, the unit disk graph $G_r(P)$ for a given $r$ is
an undirected graph, where the vertices are points of $P$, and there exists an edge
connecting two points $p, q ∈ P$, if the Euclidean distance between $p$ and $q$ is at most
$r$. The reverse shortest path problem is, for a given $λ$ and two points $s, t ∈ P$, to
compute the smallest value of $r$ such that the length of the shortest path between $s$
and $t$ in $G_r(P)$ is at most $λ$. In this repository we present three algorithms: one that runs
in $O(n log(n) + n log(r))$, the second one in $O(⌊λ⌋n log(n))$, and the third one: A brute-force baseline which runs in $O(n^2log(n))$

## Thesis

The repository contains the paper describing the problem and the implementation.

## Requirements

To build and run the project, following are required:
- CMake 3.1 or newer,
- C++17,
- CGAL

## Building The Project

To build the project run:

```
cmake .
cmake --build .
```

In the projects directory.

## Executables

There are four executable files created upon building:
- test_generator - which on input takes two values: the number of points, and the lambda value, and outputs a random test case in the format described below,
- alg_1_solver - which can be fed a test scenario, and outputs the result, solving the problem using the first algorithm,
- alg_2_solver - which can be fed a test scenario, and outputs the result,
solving the problem using the second algorithm,
- brut_solver - which can be fed a test scenario, and outputs the result,
solving the problem using the brut-force algorithm,


## Test scenario format

Files generated by the test generator, that are to be fed into the problem solving executables have the following format:
```
n λ s_i t_i

x_1 y_1
x_2 y_2
.
.
.
x_n y_n
```
In the first line four numbers separated by a white-space: $n$, $λ$, $s_i$, $t_i$, where $n$ is the number of
points, and $λ$, $s_i$ and $t_i$ is defined as in the RSP. Both indexes start at 0. Than in next $n$ lines, there are pairs of floats.
The i-th line describes the coordinates of the i-th line

## Contributions

This is a bachelor thesis of Rafał Kajca 