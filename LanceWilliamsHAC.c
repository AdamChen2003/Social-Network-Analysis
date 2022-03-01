// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// Adam Chen (z5363476)
// This program contains the functions which execute the Lance William algorithm
// for clusters. The merging processes are stored in a DNode struct.
// There are two different techniques used for this program, single or complete
// linkage. Note that for complete linkage, only complete graphs are accepted
// as input.
// COMP2521 Assignment 2

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY DBL_MAX
#define ROOT -1

static void populateArrays(int vNum, double dist[vNum][vNum], 
double distVertex[vNum][vNum], Dendrogram *dendArray, Graph g);
static void closestClusters(int vNum, double dist[vNum][vNum], 
int *indexI, int *indexJ, int clusterCount);
static Dendrogram newDendrogram(Vertex v);
static Dendrogram merge(Dendrogram minI, Dendrogram minJ);
static void removeDendrogram(Dendrogram *dendArray, int index, int clusterCount);
static double distance(Dendrogram a, Dendrogram b, int vNum, 
double distVertex[vNum][vNum], int method);
static int DendrogramToArray(Dendrogram root, Vertex *allNodes, int index);
static void updateArrays(int vNum, double dist[vNum][vNum], 
double distVertex[vNum][vNum], Dendrogram *dendArray, int method, int clusterCount);

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {
	int vNum = GraphNumVertices(g);
	double dist[vNum][vNum];
	double distVertex[vNum][vNum];
	// distVertex stores the distance between each vertex and will be used to
	// find the distance between clusters
	Dendrogram dendArray[vNum];
	populateArrays(vNum, dist, distVertex, dendArray, g);
	Dendrogram cluster = NULL;
	int clusterCount = vNum;
	// Loop continues until there is only one cluster remaining
	while (clusterCount > 1) {
		// Starting indexI at 0 and IndexJ at 1 ensures that disconnected
		// clusters can be merged if the remaining clusters are not connected
		int minIndexI = 0;
		int minIndexJ = 1;
		closestClusters(vNum, dist, &minIndexI, &minIndexJ, clusterCount);
		cluster = merge(dendArray[minIndexI], dendArray[minIndexJ]);
		// Removing the merged Dendrograms from dendArray
		removeDendrogram(dendArray, minIndexI, clusterCount);
		removeDendrogram(dendArray, minIndexJ - 1, clusterCount - 1);
		clusterCount--;
		// Inserting the merged Dendrogram into dendArray
		dendArray[clusterCount - 1] = cluster;
		updateArrays(vNum, dist, distVertex, dendArray, method, clusterCount);
	}
	return cluster;
}

// This function populates the given arrays
static void populateArrays(int vNum, double dist[vNum][vNum], 
double distVertex[vNum][vNum], Dendrogram *dendArray, Graph g) {
	// Populating the dist 2D array
	for (Vertex i = 0; i < vNum; i++) {
		AdjList edgesIn = GraphInIncident(g, i);
		AdjList edgesOut = GraphOutIncident(g, i);
		// Initialising every index of distance array with DBL_MAX
		for (int j = 0; j < vNum; j++) {
			dist[i][j] = DBL_MAX;
			distVertex[i][j] = DBL_MAX;
		}
		// Filling the distance array with the distances between i
		// and adjacent vertices
		for (AdjList temp = edgesIn; temp != NULL; temp = temp->next) {
			dist[i][temp->v] = 1.0 / temp->weight;
			distVertex[i][temp->v] = 1.0 / temp->weight;
		}
		for (AdjList temp = edgesOut; temp != NULL; temp = temp->next) {
			// Checking if the edge weight is less than the existing edge weight
			if (1.0 / temp->weight < dist[i][temp->v]) {
				dist[i][temp->v] = 1.0 / temp->weight;
				distVertex[i][temp->v] = 1.0 / temp->weight;
			}
		}
		// Populating dendArray
		dendArray[i] = newDendrogram(i);
	}
}

// This function determines which clusters are the closest together
// The index of the closest clusters are then stored in variables pointed to by
// indexI and indexJ
static void closestClusters(int vNum, double dist[vNum][vNum], 
int *indexI, int *indexJ, int clusterCount) {
	double minDist = DBL_MAX;
	for (int i = 0; i < clusterCount; i++) {
		for (int j = i; j < clusterCount; j++) {
			// If a new minimum distance has been discovered
			if (dist[i][j] < minDist) {
				minDist = dist[i][j];
				*indexI = i;
				*indexJ = j;
			}
		}
	}
}

// This function returns a malloced Dendrogram with all fields initalised
static Dendrogram newDendrogram(Vertex v) {
	Dendrogram new = malloc(sizeof(struct DNode));
	new->left = NULL;
	new->right = NULL;
	new->vertex = v;
	return new;
}

// This function returns a malloced Dendrogram which contains the two
// argument Dendrograms as the children
static Dendrogram merge(Dendrogram minI, Dendrogram minJ) {
	Dendrogram new = newDendrogram(ROOT);
	new->left = minI;
	new->right = minJ;
	return new;
}

// This function removes a Dendrogram from the given array
static void removeDendrogram(Dendrogram *dendArray, int index, int clusterCount) {
	for (int i = index; i < clusterCount - 1; i++) {
		dendArray[i] = dendArray[i + 1];
	}
}

// This function returns the distance between two clusters through either the
// single or complete linkage specified by the method
static double distance(Dendrogram a, Dendrogram b, int vNum, 
double distVertex[vNum][vNum], int method) {
	double dist = DBL_MAX;
	if (method == COMPLETE_LINKAGE) {
		dist = 0;
	}
	// Arrays which will store the verticies in the the given dendrograms, 'a' and 'b'
	// This will make the direct comparison between two Dendrograms easier
	Vertex allNodesA[vNum];
	Vertex allNodesB[vNum];
	int aNum = DendrogramToArray(a, allNodesA, 0);
	int bNum = DendrogramToArray(b, allNodesB, 0);
	// Finding the minimum/maximum distance between any two vertices from each cluster
	for (int i = 0; i < aNum; i++) {
		if (method == SINGLE_LINKAGE) {
			// For single linkage, find minimum distance
			for (int j = 0; j < bNum; j ++) {
				if (allNodesA[i] != allNodesB[j] && 
				distVertex[allNodesA[i]][allNodesB[j]] < dist) {
					dist = distVertex[allNodesA[i]][allNodesB[j]];
				}
			}
		} else if (method == COMPLETE_LINKAGE) {
			// For complete linkage, find maximum distance
			for (int j = 0; j < bNum; j ++) {
				if (allNodesA[i] != allNodesB[j] && 
				distVertex[allNodesA[i]][allNodesB[j]] > dist) {
					dist = distVertex[allNodesA[i]][allNodesB[j]];
				}
			}
		}
	}
	return dist;
}

// This function populates the given array with verticies from the 
// given Dendrogram
static int DendrogramToArray(Dendrogram root, Vertex *allNodes, int index) {
	if (root == NULL) {
		return index;
	}
	if (root->vertex != ROOT) {
		allNodes[index] = root->vertex;
		index++;
	}
	if (root->left != NULL) {
		index = DendrogramToArray(root->left, allNodes, index);
	}
	if (root->right != NULL) {
		index = DendrogramToArray(root->right, allNodes, index);
	}
	return index;
}

// This function updates the dist array with the distance between new clusters
static void updateArrays(int vNum, double dist[vNum][vNum], 
double distVertex[vNum][vNum], Dendrogram *dendArray, int method, int clusterCount) {
	for (int i = 0; i < clusterCount; i++) {
		for (int j = i; j < clusterCount; j++) {
			if (i == j) {
				dist[i][j] = DBL_MAX;
			} else {
				dist[i][j] = distance(dendArray[i], dendArray[j], vNum, distVertex, method);
			}
		}
	}
}

/**
 * Frees all memory associated with the given Dendrogram structure.
 */
void freeDendrogram(Dendrogram d) {
	if (d != NULL) {
		freeDendrogram(d->left);
		freeDendrogram(d->right);
		free(d);
	}
}

