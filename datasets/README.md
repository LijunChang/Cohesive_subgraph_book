# Graph datasets

Each graph is presented by two binary files, [b_adj.bin](example/b_adj.bin) and [b_degree.bin](example/b_degree.bin), on disk.

The contents of [b_adj.bin](example/b_adj.bin) and [b_degree.bin](example/b_degree.bin) are the same as [adj.txt](example/adj.txt) and [degree.txt](example/degree.txt), respecitvely. But the binary files only store the numbers in binary form (i.e., without the spaces and the comments).

# main.cpp

It transforms a graph in the format of edge list into our binary form.

## compible

~~
make
~~
It generates an executable "edgelist2binary"
