// Dijkstra API implementation
// Adam Chen (z5363476)
// This program contains the functions nessecary to return a ShortestPaths 
// struct containing the predecessor and distance arrays which map out shortest 
// paths from the given vertex to all other vertices in the graph. This is 
// achieved through a modified version of Dijkstra's algorithm which is able to 
// store multiple paths.
// COMP2521 Assignment 2

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"

#define INFINITY INT_MAX

static ShortestPaths *newShortestPath(Vertex src, int vCount);
static void findShortestPaths(Graph g, PQ vSet, ShortestPaths *sps);
static void updateShortestPath(ShortestPaths *sps, int minEdge, Vertex dest, Vertex prev);
static PredNode *newPredNode(Vertex w);
static void freePredNode(PredNode *pred);

// This function returns a ShortestPaths struct which contains
// the shortest paths from the 'src' vertex to all reachable nodes
// in the given graph
ShortestPaths dijkstra(Graph g, Vertex src) {
	int vCount = GraphNumVertices(g);
	ShortestPaths *spsPointer = newShortestPath(src, vCount);
	// vSet is a priotity queue which stores vertices based on adjacency
	PQ vSet = PQNew();
	// Populating the predecessor and distance arrays
	for (int i = 0; i < vCount; i++) {
		PQInsert(vSet, i, INT_MAX);
		spsPointer->pred[i] = NULL;
		spsPointer->dist[i] = INT_MAX;
	}
	PQUpdate(vSet, src, 0);
	spsPointer->dist[src] = 0;
	findShortestPaths(g, vSet, spsPointer);
	PQFree(vSet);
	ShortestPaths sps = *spsPointer;
	free(spsPointer);
	return sps;
}

// This function returns a pointer to a ShortestPaths struct
// with the fields initialised
static ShortestPaths *newShortestPath(Vertex src, int vCount) {
	ShortestPaths *sps = malloc(sizeof(struct ShortestPaths));
	sps->numNodes = vCount;
	sps->src = src;
	sps->dist = malloc(sizeof(int) * vCount);
	sps->pred = malloc(sizeof(struct PredNode) * vCount);
	return sps;
}

// This function finds the shortest paths from the src vertex, using a modified
// form of Dijkstra's algorithm
static void findShortestPaths(Graph g, PQ vSet, ShortestPaths *sps) {
	while (!PQIsEmpty(vSet)) {
		Vertex w = PQDequeue(vSet);
		// List of vertices adjacent and reachable from w
		AdjList temp = GraphOutIncident(g, w);
		while (temp != NULL) {
			// Iterating through the vertices to determine if there is a shorter path
			if (sps->dist[w] + temp->weight <= sps->dist[temp->v] 
			&& sps->dist[w] + temp->weight > 0) {
				// If we have found a suitable path, we update the ShortestPaths
				updateShortestPath(sps, temp->weight + sps->dist[w], temp->v, w);
				PQUpdate(vSet, temp->v, sps->dist[w] + temp->weight);
			}
			temp = temp->next;
		}
	}
}

// This function updates the ShortestPaths sps when a suitable path is discovered between
// vertex w and vertex prev
static void updateShortestPath(ShortestPaths *sps, int minWeight, Vertex dest, Vertex prev) {
	if (sps->dist[dest] > minWeight) {
		// If the new path is shorter than the current path
		sps->dist[dest] = minWeight;
		// we free the previously stored path
		freePredNode(sps->pred[dest]);
		sps->pred[dest] = newPredNode(prev);
	} else {
		// If the new path if of equal distance to the existing shortest path
		PredNode *temp = sps->pred[dest];
		// Iterating the find the end of the predecessor list
		if (temp != NULL) {
			while (temp->next != NULL) {
				temp = temp->next;
			}
			// Inserting the new path
			temp->next = newPredNode(prev);
		}
	}
}

// This function returns a pointer to a PredNode struct with all fields
// initialised
static PredNode *newPredNode(Vertex prev) {
	PredNode *newPredNode = malloc(sizeof(struct PredNode));
	newPredNode->v = prev;
	newPredNode->next = NULL;
	return newPredNode;
}

// This function frees all malloced fields of the given ShortestPaths struct
void freeShortestPaths(ShortestPaths sps) {
	free(sps.dist);
	for (int i = 0; i < sps.numNodes; i++) {
		freePredNode(sps.pred[i]);
	}
	free(sps.pred);
}

// This function frees the given pointer to the PredNode and all associated memory
static void freePredNode(PredNode *pred) {
	PredNode *temp = pred;
	while (temp != NULL) {
		// Creating a temporary PredNode to retain the pointer we aim to free
		// while still able to traverse the list
		PredNode *discard = temp;
		temp = temp->next;
		free(discard);
	}
}

void showShortestPaths(ShortestPaths sps) {
}
