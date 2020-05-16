# ClutteredEnvPathOpt
[![Build Status](https://travis-ci.org/mpolson64/ClutteredEnvPathOpt.jl.svg?branch=master)](https://travis-ci.org/mpolson64/ClutteredEnvPathOpt.jl)
[![codecov](https://codecov.io/gh/mpolson64/ClutteredEnvPathOpt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mpolson64/ClutteredEnvPathOpt.jl)

## Design Rationale
### types.jl
I've created a `LabeledGraph` wrapper for the `LightGraphs.AbstractGraph`, which includes an abstract graph and a `Dict` from type `T` to `Int`.
Because LightGraphs graphs' vertices must all be consecutive positive integers starting from 1, this makes it imensely easier to insert, delete, etc. while keeping things in order.
I've wrapped a number of baisc LightGraphs functions for the `LabeledGraph`, all of which are named the same as they are in LightGraphs and prefixed with `lg_`.
There is also a set of functions to convert a set of LightGraphs vertices (i.e. integers) to their corresponding label of type `T` and vice versa.

### separator.jl
This file contains a separator finding algorithm for finite element graphs and and two separator finding algorithms for planar graphs (Lipton-Tarjan and fundamental cycle separator).
Each returns a triple containing the separator first and the two non-adjacent components second and third.
The finite element graph separator and fundamental cycle separator algorithms also each have additional functions which will run the algorithim with each vertex of the source graph as the root and return the best balanced separator (these algorithms are slower and run in O(n^2) time).
The algorithms are reconstructed as close to as presented in their respective papers as possible using the LightGraphs contructs, and should behave identically in terms of asymtotic time and space complexity.
Also included is a prostprocessing algorithm to shrink the size of a separator called `pp_expell`.
It does this through a process described by Engineering Planar Separator Algorithms called node expulsion, by which a vertex is removed from the separator if it is not in connected to both A and B.

#### Note
The finite element graph separator algorithm (and its related helper functions) take in face data as a set of sets of pairs, where each set of pairs corresponds to the edges of a given face in no particular order.
In functions written later, specifically those in `biclique.jl`, I realized I needed the ordering of the vertices of every face, hence why they appear there as a set of vectors of vertices.
For simplicity, these functions could be modified to take in the face data as a set of vectors of vertices in the future.

### biclique.jl
This file is effectively just the `find_biclique_cover` function with a few helpers and testing functions also included.
The find biclique cover algorithm take in the skeleton of the graph and its face data represented as a set of clockwise-ordered vertex vectors, and first constructs the full finite element graph and checks to see if it is complete, in which case it terminates returing the two empty sets.
If it is not complete, it runs a modified version of the finite element graph separator algorithm from `separator.jl`, the modification being if either A or B is empty after postprocessing the separator algorithm is run again with a new root vertex until a suitable separator is found.
Once a separator is found we reconstruct two skeletons and face datas for subgraphs of the input embedding containing only the vertices A union C and B union C respectively.
This step is slightly more complex than a naive approach may imagine.
This is because if a face loses a vertex a face edge (i.e. an edge that was part of the input finite element graph but not part of its skeleton) will need to become a skeleton edge.
The helper function `_find_skeleton_faces` does these recalculations and returns a tuple (skeleton, faces).
With these found, we simply return the union of A => B and a recursive call of the algorithm on these two new skeletons and face datas.
A possible future optimization I expect would yield significant speedups would be to wrap the two recurive calls in futures and offloading their calculation to another thread.