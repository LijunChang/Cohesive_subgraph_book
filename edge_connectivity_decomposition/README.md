# Edge Connectivity-based Decomposition

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
It generates an executable "eco_decompose"

## 1. Run kecc or kecc-space

```
./eco_decompose ../datasets/as-skitter/ kecc 10 output
./eco_decompose ../datasets/as-skitter/ kecc-space 10 output
```
Note that, the fourth parameter is an integer that specifies the value of k. kecc-space is more space effective than kecc; that is, kecc-space consumes smaller main memory space.
