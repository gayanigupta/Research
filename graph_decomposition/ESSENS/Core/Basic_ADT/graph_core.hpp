#include <vector>

using namespace std;

void addEdge(vector<vector<int>>& graph, int u, int v) {
    graph[u].push_back(v);
    graph[v].push_back(u);
}

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