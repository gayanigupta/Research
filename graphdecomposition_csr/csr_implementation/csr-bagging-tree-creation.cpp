#include <iostream>
#include <limits.h>

using namespace std;

const int MAX_VERTICES = 100; // Maximum number of vertices
const int INF = INT_MAX;      // Infinity for initializing distances

// Function to find the minimum spanning tree (MST) connecting bags using Prim's algorithm
void primMST(const int *graph, int numVertices, int *parent)
{
    // Initialize MST set
    bool inMST[MAX_VERTICES] = {false};
    int key[MAX_VERTICES];
    for (int i = 0; i < numVertices; ++i)
    {
        key[i] = INF;
    }

    // Include first vertex in MST
    key[0] = 0;
    parent[0] = -1;

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
        for (int v = 0; v < numVertices; ++v)
        {
            if (!inMST[v] && graph[u * numVertices + v] && graph[u * numVertices + v] < key[v])
            {
                parent[v] = u;
                key[v] = graph[u * numVertices + v];
            }
        }
    }
}

// Function to create bags based on elimination order and connect them using Prim's algorithm
void createBags(const int *eliminationOrder, const int *graph, int *bags, int numVertices)
{
    // Initialize parent array for MST
    int parent[MAX_VERTICES];

    // Compute MST using Prim's algorithm
    primMST(graph, numVertices, parent);

    // Debug output to check parent array
    //cout << "Parent array: ";
    for (int i = 0; i < numVertices; ++i)
    {
        cout << parent[i] << " ";
    }
    cout << endl;

    // Map parent indices to bag indices based on elimination order
    for (int i = 1; i < numVertices; ++i)
    {
        // Get the indices of the bags corresponding to the parent and current vertex
        int bagIndex1 = -1, bagIndex2 = -1;
        for (int j = 0; j < numVertices; ++j)
        {
            if (eliminationOrder[j] == parent[i])
            {
                bagIndex1 = j;
            }
            if (eliminationOrder[j] == i)
            {
                bagIndex2 = j;
            }
        }

        // Ensure both bags are found and then print the connection
        if (bagIndex1 != -1 && bagIndex2 != -1)
        {
            // Ensure the bags are connected only once (avoid duplicates)
            if (bagIndex1 < bagIndex2)
            {
                cout << "Bags " << bagIndex1 << " and " << bagIndex2 << " are connected" << endl;
            }
        }
    }
}
