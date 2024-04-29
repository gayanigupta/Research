#include <vector>

using namespace std;


bool operator==(const Edge &e1, const Edge &e2)
{
    return ((e1.node1 == e2.node1 && e1.node2 == e2.node2) || (e1.node1 == e2.node2 && e1.node2 == e2.node1)) && e1.edge_wt == e2.edge_wt;
}

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




