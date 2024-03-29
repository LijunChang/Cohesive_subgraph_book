# Graph Dataset Format

Each graph is presented by two binary files, [b_adj.bin](example/b_adj.bin) and [b_degree.bin](example/b_degree.bin), on disk.

The contents of [b_adj.bin](example/b_adj.bin) and [b_degree.bin](example/b_degree.bin) are the same as [adj.txt](example/adj.txt) and [degree.txt](example/degree.txt), respecitvely. But the binary files only store the numbers in binary form (i.e., without the spaces and the comments).

# Converting from Edge List to the Binary Format

main.cpp converts a graph from the edge list format into our binary format.

## Compile the code
```sh
$ make clean
$ make
```
It generates an executable "edgelist2binary"

## Run the code
```sh
$ ./edgelist2binary example edges.txt
```
