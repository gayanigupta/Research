#include <iostream>

using namespace std;

// Function to create a tree in CSR notation
void treeCreation(int *eliminationOrder, int *graph, int numVertices, int *&treeIndices, int *&treeNeighbors, int &numEdges, int *&csrBagging)
{
    // Create an array to store neighbors for each node
    int **nodeNeighbors = new int *[numVertices];
    for (int i = 0; i < numVertices; i++)
    {
        nodeNeighbors[i] = new int[numVertices];
    }

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
                nodeNeighbors[node][neighborsCount++] = j;
                numEdges++; // Increment the number of edges
            }
        }

        // Update CSR indices
        treeIndices[i + 1] = treeIndices[i] + neighborsCount;

        // Debug: Print neighbors for the current node
        cout << "Node " << node << " Neighbors: ";
        for (int j = 0; j < neighborsCount; j++)
        {
            cout << nodeNeighbors[node][j] << " ";
        }
        cout << endl;

        // Update CSR Bagging
        csrBagging[i] = treeIndices[i];
    }

    // Add the final CSR Bagging entry
    csrBagging[numVertices] = numEdges;

    // Create the CSR neighbors array and copy neighbors into it
    treeNeighbors = new int[numEdges];
    int neighborIndex = 0;
    for (int i = 0; i < numVertices; i++)
    {
        for (int j = treeIndices[i]; j < treeIndices[i + 1]; j++)
        {
            treeNeighbors[neighborIndex++] = nodeNeighbors[i][j - treeIndices[i]];
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
    for (int i = 0; i < numVertices; i++)
    {
        delete[] nodeNeighbors[i];
    }
    delete[] nodeNeighbors;
}

// Function to create bags of vertices and link bags based on common vertices
void createAndLinkBags(int *graph, int numVertices)
{
    // Define the size of each bag (number of vertices per bag)
    int bagSize = 2; // You can adjust the bag size as needed

    int maxCommonVertices = 0;
    int bestBag1 = -1;
    int bestBag2 = -1;

    for (int i = 0; i < numVertices; i += bagSize)
    {
        for (int j = i + bagSize; j < numVertices; j += bagSize)
        {
            int commonVertices = 0;

            for (int v1 = i; v1 < i + bagSize; v1++)
            {
                for (int v2 = j; v2 < j + bagSize; v2++)
                {
                    if (graph[v1 * numVertices + v2])
                        commonVertices++;
                }
            }

            if (commonVertices > maxCommonVertices)
            {
                maxCommonVertices = commonVertices;
                bestBag1 = i;
                bestBag2 = j;
            }
        }
    }

    // Implement your logic to link bags based on the bestBag1 and bestBag2
    // Modify the graph or create a new graph to represent bag connections
    if (bestBag1 >= 0 && bestBag2 >= 0)
    {
        // Link bags by modifying the graph
        for (int v1 = bestBag1; v1 < bestBag1 + bagSize; v1++)
        {
            for (int v2 = bestBag2; v2 < bestBag2 + bagSize; v2++)
            {
                graph[v1 * numVertices + v2] = 1;
                graph[v2 * numVertices + v1] = 1;
            }
        }
    }

    // Debug: Print the bags with the maximum common vertices
    cout << "Bags with Max Common Vertices: " << bestBag1 << " " << bestBag2 << endl;
}

int main()
{
    int numVertices = 5; // Adjust the number of vertices as needed

    // Define the elimination order
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
    int *csrBagging = new int[numVertices + 1];

    // Call the treeCreation function to create the tree in CSR notation
    treeCreation(eliminationOrder, graph, numVertices, treeIndices, treeNeighbors, numEdges, csrBagging);

    // Create and link bags based on common vertices
    createAndLinkBags(graph, numVertices);

    // Clean up memory
    delete[] treeIndices;
    delete[] treeNeighbors;
    delete[] csrBagging;

    return 0;
}
