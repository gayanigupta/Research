#include <iostream>

using namespace std;

// Function to create a tree in CSR notation
void treeCreation(int *eliminationOrder, int *graph, int numVertices, int *&treeIndices, int *&treeNeighbors, int &numEdges)
{
    // Create an array to store neighbors for each node
    int *nodeNeighbors = new int[numVertices];

    // Initialize CSR indices
    treeIndices = new int[numVertices + 1];
    treeIndices[0] = 0;

    // Initialize the number of edges
    numEdges = 0;

    // Iterate through the elimination order
    for (int i = 0; i < numVertices; i++)
    {
        int node = eliminationOrder[i];
        int neighborsCount = 0;

        // Find neighbors of the current node
        for (int j = 0; j < numVertices; j++)
        {
            if (graph[node * numVertices + j])
            {
                nodeNeighbors[neighborsCount++] = j;
                numEdges++; // Increment the number of edges
            }
        }

        // Update CSR indices
        treeIndices[i + 1] = treeIndices[i] + neighborsCount;

        // Debug: Print neighbors for the current node
        cout << "Node " << node << " Neighbors: ";
        for (int j = 0; j < neighborsCount; j++)
        {
            cout << nodeNeighbors[j] << " ";
        }
        cout << endl;
    }

    // Create the CSR neighbors array and copy neighbors into it
    treeNeighbors = new int[numEdges];
    int neighborIndex = 0;

    for (int i = 0; i < numVertices; i++)
    {
        for (int j = treeIndices[i]; j < treeIndices[i + 1]; j++)
        {
            treeNeighbors[neighborIndex++] = nodeNeighbors[j - treeIndices[i]];
        }
    }

    // Debug: Print CSR Indices
    cout << "CSR Indices: ";
    for (int i = 0; i <= numVertices; i++)
    {
        cout << treeIndices[i] << " ";
    }
    cout << endl;

    // Debug: Print CSR Neighbors
    cout << "CSR Neighbors: ";
    for (int i = 0; i < numEdges; i++)
    {
        cout << treeNeighbors[i] << " ";
    }
    cout << endl;

    // Clean up memory
    delete[] nodeNeighbors;
}

int main()
{
    int numVertices = 5; // Adjust the number of vertices as needed

    // Define the elimination order and graph (adjacency matrix)
    int eliminationOrder[] = {2, 3, 1, 4, 0};
    int graph[] = {
        0, 1, 0, 0, 0,
        1, 0, 1, 1, 0,
        0, 1, 0, 0, 0,
        0, 1, 0, 0, 1,
        0, 0, 0, 1, 0
    };

    // Variables to store CSR representation
    int *treeIndices;
    int *treeNeighbors;
    int numEdges;

    // Call the treeCreation function to create the tree in CSR notation
    treeCreation(eliminationOrder, graph, numVertices, treeIndices, treeNeighbors, numEdges);

    return 0;
}
