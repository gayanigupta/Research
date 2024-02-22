#include <iostream>

using namespace std;


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

// Function to create bags based on elimination order
void createBags(const int *eliminationOrder, const int *graph, int *bags, int numVertices)
{
    int bagIndex = 0;

    for (int i = 0; i < numVertices; ++i)
    {
        int vertex = eliminationOrder[i];

        // Array to store neighbors of the current vertex
        int neighbors[100]; // Assuming a maximum of 100 neighbors (adjust as needed)
        int numNeighbors = 0;

        findNeighbors(vertex, graph, numVertices, neighbors, numNeighbors);

        // Add the vertex and its neighbors to the bag
        bags[bagIndex++] = vertex;
        for (int j = 0; j < numNeighbors; ++j)
        {
            bags[bagIndex++] = neighbors[j];
        }

        // Link bags based on common vertices
        for (int k = 0; k < i; ++k)
        {
            int prevBagIndex = k;
            int commonVertices[100]; // Assuming a maximum of 100 common vertices (adjust as needed)
            int numCommon = 0;

            // Find common vertices between current bag and previous bag
            for (int l = i + 1; l < bagIndex; ++l)
            {
                int currentVertex = bags[l];
                if (bags[prevBagIndex] == currentVertex)
                {
                    commonVertices[numCommon++] = currentVertex;
                }
            }

            // Update bag links
            if (numCommon > 0)
            {
                // Connect bags with common vertices
                // For simplicity, printing the common vertices here
                std::cout << "Bags " << k << " and " << i << " are linked with common vertices: ";
                for (int m = 0; m < numCommon; ++m)
                {
                    std::cout << commonVertices[m] << " ";
                }
                std::cout << std::endl;
            }
        }
    }
}
