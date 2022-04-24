# Core Decomposition Algorithms

This repository implements the algorithms presented in our monograph “Cohesive Subgraph Computation over Large Sparse Graphs”. If you are using the code, please cite our monograph.
<pre>
Lijun Chang and Lu Qin.
<a href="https://www.springer.com/us/book/9783030035983">Cohesive Subgraph Computation over Large Sparse Graphs.</a>
Springer Series in the Data Sciences, 2018
</pre>

## Compile the code
```sh
$ make clean
$ make
```
It generates an executable "core_decompose"

## Run the code
```sh
$ ./core_decompose ../datasets/as-skitter/ list-heap output
$ ./core_decompose ../datasets/as-skitter/ array-heap output
$ ./core_decompose ../datasets/as-skitter/ h-index output
$ ./core_decompose ../datasets/as-skitter/ output
$ ./core_decompose ../datasets/as-skitter/ hierarchy output
```
Note that, parameter "output" is optional
* list-heap uses linked list-based linear heap for core decomposition
* array-heap uses array-based linear heap for core decomposition
* h-index uses h-index-based core decomposition
* The third one is similar to array-heap, but directly implements the logic of the data structure within core decomposition without explicitly constructing a new class
* hierarchy constructs the core spanning tree
