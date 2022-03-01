// Centrality Measures API implementation
// Adam Chen (z5363476)
// This program contains the functions which facilitate the calculations of
// the betweenness and closeness centralities. These centralities are stored 
// in the returned NodeValues struct.
// COMP2521 Assignment 2

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "PQ.h"

#define INFINITY INT_MAX

static double closenessCalculation(Graph g, ShortestPaths sps);
static double betweennessCalculation(Graph g, Vertex v, int nodeCount);
static double vertexBetweenness(Graph g, Vertex v, Vertex start, int nodeCount);
static double pathCount(Vertex dest, Vertex start, ShortestPaths sps);
static void recursiveCount(Vertex start, Vertex dest, int *count, int *visited, PQ currentPath, ShortestPaths sps);
static double vPathCount(ShortestPaths sps, Vertex v, Vertex start, Vertex dest);

// This function returns a NodeValues struct which contains the closeness
// centrality of all the nodes in the graph
NodeValues closenessCentrality(Graph g) {
	NodeValues nvs = {0};
	nvs.numNodes = GraphNumVertices(g);
	nvs.values = malloc(sizeof(double) * nvs.numNodes);
	for (int v = 0; v < nvs.numNodes; v++) {
		ShortestPaths sps = dijkstra(g, v);
		nvs.values[v] = closenessCalculation(g, sps);
		freeShortestPaths(sps);
	}
	return nvs;
}

// This function returns the closeness centrality of the given vertex in the 
// given graph
static double closenessCalculation(Graph g, ShortestPaths sps) {
	double totalDist = 0;
	double reachableVertices = 0;
	for (int i = 0; i < sps.numNodes; i++) {
		if (sps.dist[i] < INFINITY) {
			// If the vertex is reachable from v
			totalDist += sps.dist[i];
			reachableVertices++;
		}
	}
	if (reachableVertices == 1) {
		// If the vertex is only reachable from itself
		return 0;
	}
	return (reachableVertices - 1) / (sps.numNodes - 1) * (reachableVertices - 1) / totalDist;
}

// This function returns a NodeValues struct which contains the betweenness
// centrality of all the nodes in the graph
NodeValues betweennessCentrality(Graph g) {
	NodeValues nvs = {0};
	nvs.numNodes = GraphNumVertices(g);
	nvs.values = malloc(sizeof(double) * nvs.numNodes);
	for (int v = 0; v < nvs.numNodes; v++) {
		nvs.values[v] = betweennessCalculation(g, v, nvs.numNodes);
	}
	return nvs;
}

// This function returns the betweenness centrality of the vertex 'v' in the 
// given graph
static double betweennessCalculation(Graph g, Vertex v, int nodeCount) {
	double totalBetweenness = 0;
	for (int start = 0; start < nodeCount; start++) {
		if (start != v) {
			totalBetweenness += vertexBetweenness(g, v, start, nodeCount);
		}
	}
	return totalBetweenness;
}

// This is a helper function for betweenessCalculation and completes the
// betweenness centrality from the start vertex to all other vertices in
// the graph
static double vertexBetweenness(Graph g, Vertex v, Vertex start, int nodeCount) {
	double betweenness = 0;
	ShortestPaths sps = dijkstra(g, start);
	for (int dest = 0; dest < nodeCount; dest++) {
		if (dest != v && sps.dist[dest] < INFINITY) {
			// If the dest vertex is not v and is reachable 
			// from the start vertex
			double totalPaths = pathCount(dest, start, sps);
			double vPaths = vPathCount(sps, v, start, dest);
			if (totalPaths != 0) {
				betweenness += (vPaths / totalPaths);
			}
		}
	}
	freeShortestPaths(sps);
	return betweenness;
}

// This function returns the number of paths from vertex v to w by using a dfs
static double pathCount(Vertex dest, Vertex start, ShortestPaths sps) {
	int count = 0;
	// Following lines mallocs a visited array and sets all vertices to 'unvisited'
	int *visited = malloc(sizeof(int) * sps.numNodes);
	for (int i = 0; i < sps.numNodes; i++) {
		visited[i] = false;
	}
	// currentPath will be the stack
	PQ currentPath = PQNew();
	recursiveCount(start, dest, &count, visited, currentPath, sps);
	free(visited);
	PQFree(currentPath);
	return count;
}

// This is a helper function for pathCount.
// This function determines the number of total paths from start to dest using recursion
static void recursiveCount(Vertex start, Vertex dest, int *count, int *visited, PQ currentPath, ShortestPaths sps) {
	if (!visited[dest]) {
		visited[dest] = true;
		PQInsert(currentPath, dest, *count);
		if (dest == start) {
			*count += 1;
			visited[dest] = false;
			PQDequeue(currentPath);
		}
		for (PredNode *temp = sps.pred[dest]; temp != NULL; temp = temp->next) {
			recursiveCount(start, temp->v, count, visited, currentPath, sps);
		}
		visited[dest] = false;
	}
}

// This function returns the number of paths from i to j which contain v
static double vPathCount(ShortestPaths sps, Vertex v, Vertex start, Vertex dest) {
	return pathCount(dest, v, sps) * pathCount(v, start, sps);
}

// This function returns a NodeValues struct which contains the normalised 
// betweenness centrality of all the nodes in the graph
NodeValues betweennessCentralityNormalised(Graph g) {
	NodeValues nvs = {0};
	nvs.numNodes = GraphNumVertices(g);
	nvs.values = malloc(sizeof(double) * nvs.numNodes);
	for (int v = 0; v < nvs.numNodes; v++) {
		nvs.values[v] = betweennessCalculation(g, v, nvs.numNodes) / (nvs.numNodes - 1) / (nvs.numNodes - 2);
	}
	return nvs;
}

void showNodeValues(NodeValues nvs) {
}

void freeNodeValues(NodeValues nvs) {
	free(nvs.values);
}

