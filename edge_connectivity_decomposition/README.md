# Edge Connectivity-based Graph Decomposition Algorithms

## Compile

```
make
```
It generates an executable "eco_decompose"

## Run kecc

```
./eco_decompose ../datasets/as-skitter/ kecc 10 output
```
Note that, the fourth parameter is an integer that specifies the value of k.

This implements the algorithm proposed in the following paper.

Lijun Chang, Jeffrey Xu Yu, Lu Qin, Xuemin Lin, Chengfei Liu, and Weifa Liang <br/>
**Efficiently Computing k-Edge Connected Components via Graph Decomposition** <br />
*Proceedings of the ACM SIGMOD International Conference on Management of Data* (SIGMOD’13), 2013
