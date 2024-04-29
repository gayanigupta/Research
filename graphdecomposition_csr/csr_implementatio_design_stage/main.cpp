#include <iostream>
#include <fstream>
#include "elimination_order_csr.hpp"

// Function to read the input file and process the data.
void readInputFile(const std::string &filename, int &numEdges, int *graph, int &numVertices)
{
    std::ifstream inputFile(filename);
    if (!inputFile.is_open())
    {
        std::cerr << "Failed to open the input file: " << filename << std::endl;
        return;
    }
    // std::cout<<"I am here \n";
    // std::cout<<filename<<"\n";

    int EIndex = numVertices;
    int currentNode = -1;

    // Read each line of the input file.
    std::cout << "printing the csr format vertex position in the graph:" << std::endl;
    while (!inputFile.eof())
    {
        int vertex1, vertex2;
        double weight;
        inputFile >> vertex1 >> vertex2 >> weight;
        // std::cout<< "vertex1:  "<<vertex1<< std::endl;
        // std::cout<< "vertex2: "<<vertex2<< std::endl;
        // std::cout<< vertex2<< std::endl;

        // If vertex1 is not equal to currentNode, update currentNode to vertex1.
        if (vertex1 != currentNode)
        {
            currentNode = vertex1;
            graph[currentNode] = EIndex;
            std::cout << currentNode << "::" << EIndex << "\n";
        }
        graph[EIndex] = vertex2;

        EIndex++;
    }
    std::cout << "******************************************************" << std::endl;
    for (int i = 0; i <= numVertices + numEdges; i++)
    {
        std::cout << graph[i] << " " << i << std::endl;
    }

    inputFile.close();
}

/*int main(int argc, char *argv[])*/
int main()
{
    // TODO: Graph size is not right fix it.
    //  std::string filename ="/Users/gayanigupta/Documents/GitHub/Research/graph_decomposition _refined/Example/Tests/core_3.txt";
    int *graph = nullptr;
    int *eliminationOrder = nullptr;
    //int numVertices, numEdges;

     int numVertices = 3;
     int numEdges = 3;
    // std::cout << argc << std::endl;
   /* if (argc < 4)
    {
        std::cerr << "Usage:" << argv[1] << "file, numVertices,numEdges" << std::endl;
        return 1;
    }*/
   // const char *filename = argv[1];
   // numVertices = std::atoi(argv[2]);
   // numEdges = std::atoi(argv[3]);
    const char *filename = "data/core_3.txt";
    // Create an array to hold the graph data.

    // TODO: Initialize with -1 automatically
    graph = new int[numVertices + numEdges + 4];
    //std::fill_n(eliminationOrder, numVertices, -1);

   // eliminationOrder = new int[numEdges]{5,4};

   /*for (int i = 0; i <= numVertices + numEdges; ++i)
    {
        graph[i] = -1;
    }*/

    // Call the function to read and process the input file.
    readInputFile(filename, numEdges, graph, numVertices);

    // DEBUG :: Printing the graph

   /* for (int i = 0; i <= numVertices + numEdges; i++)
    {
        if (i <= numVertices + numEdges)
        {
            std::cout << graph[i] << " " << i << std::endl;
        }
    }*/
    // Print the number of vertices and edges.
    std::cout << "Number of Vertices: " << numVertices << std::endl;
    std::cout << "Number of Edges: " << numEdges << std::endl;

    std::cout << "size of graph: " << (numVertices + numEdges + 4) << " elements" << std::endl;

    //findChordalEdgesWithEliminationOrder(graph, numVertices, eliminationOrder, numEdges);

    int *bags = new int[numVertices * (numVertices + 1)]; // Assuming a maximum of numVertices bags

    // Clean up allocated memory.
    // delete[] graph;

    return 0;
}
