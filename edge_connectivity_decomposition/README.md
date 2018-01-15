# Edge Connectivity-based Graph Decomposition Algorithms

## Compile

```
make
```
It generates an executable "eco_decompose"

## 1. Run kecc or kecc-space

```
./eco_decompose ../datasets/as-skitter/ kecc 10 output
./eco_decompose ../datasets/as-skitter/ kecc-space 10 output
```
Note that, the fourth parameter is an integer that specifies the value of k. kecc-space is more space effective than kecc; that is, kecc-space consumes smaller main memory space.

This implements the algorithm proposed in the following SIGMOD'13 paper, which computes all k-edge connected components of a graph for a given k.

Lijun Chang, Jeffrey Xu Yu, Lu Qin, Xuemin Lin, Chengfei Liu, and Weifa Liang <br/>
**Efficiently Computing k-Edge Connected Components via Graph Decomposition** <br/>
*Proceedings of the ACM SIGMOD International Conference on Management of Data* (SIGMOD’13), 2013
