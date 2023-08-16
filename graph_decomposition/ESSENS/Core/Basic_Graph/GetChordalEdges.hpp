#pragma once
// Library to find the chordal edges with a given elimination order
#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

using namespace std;

void findChordalEdgesWithEliminationOrder(const vector<vector<int>>& graph, const vector<int>& eliminationOrder) {
    int n = graph.size();
    vector<unordered_set<int>> neighbors(n);

    // Store neighbors of each vertex
    for (int u = 0; u < n; ++u) {
        for (int v : graph[u]) {
            neighbors[u].insert(v);
        }
    }

    for (int i = 0; i < n; ++i) {
        int u = eliminationOrder[i];
        for (int j = i + 1; j < n; ++j) {
            int v = eliminationOrder[j];

            // Check if v is a neighbor of u
            if (neighbors[u].count(v) == 1) {
                // Check if v is adjacent to all the previous vertices in the elimination order
                bool isChordalEdge = true;
                for (int k = i + 1; k < j; ++k) {
                    int w = eliminationOrder[k];
                    if (neighbors[v].count(w) == 0) {
                        isChordalEdge = false;
                        break;
                    }
                }

                if (isChordalEdge) {
                    cout << "Chordal Edge: (" << u << ", " << v << ")" << endl;
                }
            }
        }
    }
}