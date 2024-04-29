#include <iostream>
#include <limits.h>

using namespace std;

const int MAX_VERTICES = 100; // Maximum number of vertices
const int INF = INT_MAX; // Infinity for initializing distances

// Function to find neighbors of a vertex in CSR graph
void findNeighbors(int vertex, const int *graph, int numVertices, int *neighbors, int &numNeighbors)
{
    numNeighbors = 0;

    int start = graph[vertex];
    int end = (vertex == numVertices - 1) ? (numVertices + graph[numVertices - 1]) : graph[vertex + 1];

    for (int i = start; i < end; ++i)
    {
        neighbors[numNeighbors++] = graph[i];
    }
}

// Function to create bags based on elimination order and connect them using Prim's algorithm
void createBags(const int *eliminationOrder, const int *graph, int *bags, int numVertices)
{
    int bagIndex = 0;

    // Initialize MST set
    bool inMST[MAX_VERTICES] = {false};
    int parent[MAX_VERTICES];
    int key[MAX_VERTICES];
    for (int i = 0; i < numVertices; ++i)
    {
        key[i] = INF;
    }

    // Include first vertex in MST
    key[eliminationOrder[0]] = 0;
    parent[eliminationOrder[0]] = -1;

    for (int count = 0; count < numVertices - 1; ++count)
    {
        // Find vertex with minimum key value
        int minKey = INF;
        int minIndex;
        for (int v = 0; v < numVertices; ++v)
        {
            if (!inMST[v] && key[v] < minKey)
            {
                minKey = key[v];
                minIndex = v;
            }
        }

        int u = minIndex;
        inMST[u] = true;

        // Update key values and parent index of adjacent vertices of the picked vertex
        int neighbors[MAX_VERTICES]; // Assuming a maximum of numVertices neighbors (adjust as needed)
        int numNeighbors = 0;
        findNeighbors(u, graph, numVertices, neighbors, numNeighbors);
        for (int j = 0; j < numNeighbors; ++j)
        {
            int v = neighbors[j];
            if (graph[u] != graph[v] && !inMST[v] && graph[u] != numVertices && graph[v] != numVertices && graph[u + 1] != graph[v + 1])
            {
                if (graph[v] < key[v])
                {
                    parent[v] = u;
                    key[v] = graph[v];
                }
            }
        }
    }

    // Connect bags based on MST edges
    for (int i = 1; i < numVertices; ++i)
    {
        // For simplicity, printing the connected bags here
        std::cout << "Bags " << parent[i] << " and " << i << " are connected with edge weight: " << key[i] << std::endl;
    }
}

int main() {
    // Example usage
    int eliminationOrder[MAX_VERTICES] = {0, 1, 2, 3}; // Example elimination order
    int graph[MAX_VERTICES] = {0, 2, 3, 5, 7}; // Example graph in CSR format
    int bags[MAX_VERTICES]; // Array to store bags
    int numVertices = 4; // Number of vertices in the graph

    createBags(eliminationOrder, graph, bags, numVertices);

    return 0;
}
