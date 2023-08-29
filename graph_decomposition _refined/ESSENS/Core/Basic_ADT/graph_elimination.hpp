#include <iostream>
#include <vector>
#include <unordered_set>

//#include "graph_core.hpp"

using namespace std;

// Function to find the perfect elimination order using the Lex-M algorithm
vector<int> findPerfectEliminationOrder(const vector<vector<int>>& graph) {
    int n = graph.size();
    vector<bool> eliminated(n, false);
    vector<int> eliminationOrder;

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
            exit(1);
        }

        eliminationOrder.push_back(lowestParent);
        eliminated[lowestParent] = true;
    }

    return eliminationOrder;
}

// Function to find the chordal edges with a given elimination order
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