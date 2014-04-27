SHELL = /bin/sh
SYSTEM = $(shell uname)
C++ = g++
CC = gcc
DFLAGS = 
OFLAGS = -g -O3
LFLAGS = -llpsolve55
CFLAGS = $(OFLAGS) $(DFLAGS) -I. -I/home/wizardus/Desktop/lpsolve -Wno-deprecated

OBJS = bit_vector.o branchPath.o brutePath.o Contraction.o Dijkstra.o Graph.o graphlink.o GraphUtil.o HD.o lpPath.o oneTaxiPath.o path.o PerformanceTimer.o Query.o shortestPath.o simulator.o taxiPath.o taxitree.o treeTaxiPath.o treeSlackTaxiPath.o treeClusterTaxiPath.o Util.o
PROGS = ./simulator

all: $(OBJS) $(COBJS) $(PROGS)

./simulator: $(OBJS) $(COBJS)
	$(C++) -o ./simulator $(OBJS) $(COBJS) $(LFLAGS)

clean:
	rm -f $(OBJS) $(COBJS) $(PROGS)

$(OBJS): %.o: %.cpp
	$(C++) -o $@ $(CFLAGS) -c $<

$(COBJS): %.o: %.c
	$(CC) -o $@ $(CFLAGS) -c $<

./simulator: $(OBJS) $(COBJS)

all: $(PROGS)

bit_vector.o: bit_vector.h
branchPath.o: path.h taxiPath.h branchPath.h
brutePath.o: path.h taxiPath.h brutePath.h
Contraction.o: GraphUtil.h Contraction.h
Dijkstra.o: Graph.h Dijkstra.h
Graph.o: Definition.h Util.h Graph.h
graphlink.o: Contraction.h Query.h Dijkstra.h PerformanceTimer.h GraphUtil.h HD.h graphlink.h
GraphUtil.o: Graph.h GraphUtil.h
HD.o: Graph.h HD.h
lpPath.o: path.h taxiPath.h lpPath.h
oneTaxiPath.o: includes.h path.h shortestPath.h vertex.h taxiPath.h branchPath.h brutePath.h lpPath.h oneTaxiPath.h
path.o: includes.h taxiPath.h path.h
PerformanceTimer.o: PerformanceTimer.h
Query.o: Graph.h Query.h
shortestPath.o: includes.h Contraction.h Query.h Dijkstra.h PerformanceTimer.h GraphUtil.h HD.h graphlink.h PerformanceTimer.h vertex.h shortestPath.h
simulator.o: includes.h taxitree.h customerqueue.h vertex.h shortestPath.h taxilinklist.h
taxiPath.o: includes.h path.h shortestPath.h vertex.h taxiPath.h
taxitree.o: includes.h taxiPath.h location.h vertex.h shortestPath.h taxitree.h
treeTaxiPath.o: includes.h path.h shortestPath.h vertex.h taxiPath.h treeTaxiPath.h
treeSlackTaxiPath.o: includes.h path.h shortestPath.h vertex.h taxiPath.h treeSlackTaxiPath.h
treeClusterTaxiPath.o: includes.h path.h shortestPath.h vertex.h taxiPath.h treeClusterTaxiPath.h
Util.o: Config.cpp bit_vector.h PerformanceTimer.h Util.h
