#include <iostream>
#include <fstream>
#include "csr-bagging-tree-creation.cpp"
using namespace std;

// Function to find chordal edges with elimination order.
void findChordalEdgesWithEliminationOrder(int *graph, int numVertices, int *eliminationOrder, int numEdges)
{
    std::cout << "I am in function" << std::endl;

    // Initialize lowestParents with -1 for all vertices
    int *lowestParents = new int[numVertices];
    for (int i = 0; i < numVertices; i++)
    {
        lowestParents[i] = -1;
    }
    // Initialize chordalEdges with -1 for all vertices
    bool *chordalEdges = new bool[numEdges];
    for (int i = 0; i < numEdges; i++)
    {
        chordalEdges[i] = false;
    }
    // go through the neighbors

    for (int i = 0; i < numVertices; i++)
    {
        int stopIndex = graph[i + 1];

        if (i == numVertices - 1)
        {
            stopIndex = numVertices + numEdges;
        }

        int minNeghbor = numVertices;

        for (int j = graph[i]; j < stopIndex; j++)
        {
            int neighbor = graph[j];
            if (neighbor < minNeghbor)
            {
                minNeghbor = neighbor;
            }
            std::cout << i << " " << neighbor << std::endl;
        }

        if (minNeghbor < i)
        {
            lowestParents[i] = minNeghbor;

            // Check if the neighbors of the lowest parent are a subset of the current node's neighbors
            int size1 = graph[i + 1] - graph[i];
            int size2 = graph[lowestParents[i] + 1] - graph[lowestParents[i]];
            bool isSubset = true;

            std::cout << "i " << i << std::endl;
            std::cout << "array 1" << std::endl;
            for (int a = graph[i]; a < (graph[i + 1]); a++)
            {
                std::cout << graph[a] << " ";
            }
            std::cout << "\n"
                      << std::endl;
            std::cout << "lowestParent " << lowestParents[i] << std::endl;
            std::cout << "array 2 :: " << std::endl;
            int k = lowestParents[i];
            for (int b = graph[k]; b < graph[k + 1]; b++)
            {
                std::cout << graph[b] << " ";
            }
            std::cout << "\n"
                      << std::endl;

            // Check if the neighbors of the lowest parent are a subset of the current node's neighbors then mark all 
            // consider the edges the second array(current node's neighbors) as chordal
            int indexArr1 = graph[i];
            int indexArr2 = graph[lowestParents[i]];

            for (int a = graph[i]; a < (graph[i + 1]); a++)
            {
                int neighbor1 = graph[a];
                bool found = false;
                for (int b = graph[lowestParents[i]]; b < graph[lowestParents[i] + 1]; b++)
                {
                    int neighbor2 = graph[b];
                    if (neighbor1 == neighbor2)
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    isSubset = false;
                    break;
                }
            }

            if (isSubset)
            {
                // Mark edges as chordal
                for (int a = graph[i]; a < (graph[i + 1]); a++)
                {
                    int neighbor1 = graph[a];
                    for (int b = graph[lowestParents[i]]; b < graph[lowestParents[i] + 1]; b++)
                    {
                        int neighbor2 = graph[b];
                        if (neighbor1 == neighbor2)
                        {
                            int edgeIndex = -1;
                            for (int j = graph[i]; j < graph[i + 1]; j++)
                            {
                                if (graph[j] == neighbor2)
                                {
                                    edgeIndex = j - graph[i];
                                    break;
                                }
                            }
                            if (edgeIndex != -1)
                            {
                                chordalEdges[graph[i] + edgeIndex] = true;
                            }
                        }
                    }
                }
            }
        }
    }

    // Print chordal edges
    for (int i = 0; i < numEdges; i++)
    {
        if (chordalEdges[i])
        {
            std::cout << "Chordal Edge: " << i << std::endl;
        }
    }


    

    // Clean up dynamic memory
    delete[] lowestParents;
    delete[] chordalEdges;

    int *bags = new int[numVertices * (numVertices + 1)]; // Assuming a maximum of numVertices bags

    // Call the function to create bags and link them
   //createBags(eliminationOrder, graph, bags, numVertices);

    // Clean up dynamic memory
    delete[] bags;
}
