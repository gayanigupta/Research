#pragma once
// Library to find the perfect elimination order using the Lex-M algorithm
#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

using namespace std;

// Function to calculate the lowest parent of a vertex in the graph
int findLowestParent(const vector<vector<int>>& graph, const vector<bool>& eliminated, int vertex) {
    int lowestParent = vertex;

    for (int neighbor : graph[vertex]) {
        if (!eliminated[neighbor] && neighbor < lowestParent) {
            lowestParent = neighbor;
        }
    }

    return lowestParent;
}

bool findPerfectEliminationOrder(vector<vector<int>>& graph, vector<int>& eliminationOrder)
{
    int n = graph.size();
    vector<bool> eliminated(n, false);

    for (int i = 0; i < n; ++i) {
        int lowestParent = -1;

        // Find the lowest parent that has not been eliminated
        for (int v = 0; v < n; ++v) {
            if (!eliminated[v]) {
                int parent = findLowestParent(graph, eliminated, v);
                if (lowestParent == -1 || parent < lowestParent) {
                    lowestParent = parent;
                }
            }
        }

        // Check if removing the lowest parent makes the remaining graph a clique
        bool isClique = true;
        for (int v = 0; v < n; ++v) {
            if (!eliminated[v] && v != lowestParent) {
                if (find(graph[v].begin(), graph[v].end(), lowestParent) == graph[v].end()) {
                    isClique = false;
                    break;
                }
            }
        }

        if (!isClique) {
            cerr << "The input graph is not chordal. Please ensure it is a chordal graph before using this algorithm." << endl;
            //exit(1);
            return false;
        }

        eliminationOrder.push_back(lowestParent);
        eliminated[lowestParent] = true;
    }

    return true;
}