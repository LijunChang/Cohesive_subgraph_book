# Core Decomposition Algorithms

## compile

```
make
```
It generates an executable "core_decompose"

## run

```
./core_decompose ../datasets/as-skitter/ list-heap output
./core_decompose ../datasets/as-skitter/ array-heap output
./core_decompose ../datasets/as-skitter/ h-index output
./core_decompose ../datasets/as-skitter/ output
./core_decompose ../datasets/as-skitter/ hierarchy output
```

*list-heap uses linked list-based linear heap for core decomposition
*array-heap uses array-based linear heap for core decomposition
*h-index uses h-index-based core decomposition
*The third one is similar to array-heap, but directly implements the logic of the data structure within core decomposition without explicitly constructing a new class
*hierarchy constructs the core spanning tree
