# ********     Author: Lijun Chang    ******
# ******** Email: ljchang@outlook.com ******
#
CC=g++ -flto -O3
CFLAGS=-c -I. -std=c++11

all: edgelist2binary

edgelist2binary: .obj/main.o
	${CC} .obj/main.o -o edgelist2binary
	rm .obj/*.o

.obj/main.o: main.cpp
	${CC} ${CFLAGS} -o .obj/main.o main.cpp

clean:
	rm -rf *o .obj/
	mkdir .obj
