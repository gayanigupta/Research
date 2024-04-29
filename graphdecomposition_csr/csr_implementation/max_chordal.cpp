#include <iostream>
#include <fstream>

const int MAX_VERTICES = 100;//adjust the maximum number of vertices

int LP[MAX_VERTICES];
int rows[MAX_VERTICES + 1];
int cols[MAX_VERTICES * MAX_VERTICES];
int chordalNeighbors[MAX_VERTICES][MAX_VERTICES];


void processVertices(int vertex)
{
    // Process each vertex to find its chordal neighbors
    for (int i = rows[vertex]; i < rows[vertex + 1]; ++i)
    {
        int neighbor = cols[i];
        if (LP[neighbor] == vertex)
        {
            if (chordalNeighbors[vertex][neighbor] == 0)
            {
                chordalNeighbors[vertex][neighbor] = 1;
                LP[neighbor] = 0; // Set LP of w to its next lowest parent
            }
        }
    }
}

void maximalChordalSubgraph(int numVertices, int *graph)
{
    int Q1[MAX_VERTICES], Q2[MAX_VERTICES];
    int Q1_size = 0, Q2_size = 0;

    // Add vertices to Q1 that are LP to at least one other vertex
    for (int v = 0; v < numVertices; ++v)
    {
        LP[v] = -1; // assign default value

        if (v > 0)
        {
            int min_neigh = v;
            std::cout << "v: " << v << std::endl;
            std::cout << "graph[v]: " << graph[v] << std::endl;
            std::cout << "graph[v + 1]: " << graph[v + 1] << std::endl;
            for (int i = graph[v]; i < graph[v + 1]; ++i)
            {
                int neighbor = graph[i];
                std::cout << "i: " << i << ", neighbor: " << neighbor << std::endl;
                if (neighbor < min_neigh) // update minimum neighbor
                {
                    min_neigh = neighbor;
                }
                if (neighbor >= min_neigh) // becos vertices in increasing order, thus break at first largest
                {
                    break;
                }
            }
            LP[v] = min_neigh;
        } // end of if
    }
    std::cout << "LP: ";
    for (int v = 0; v < numVertices; ++v)
    {
        std::cout << LP[v] << " ";
    }
    std::cout << std::endl;

    while (Q1_size > 0)
    {
        for (int i = 0; i < Q1_size; ++i)
        {
            int vertex = Q1[i];
            processVertices(vertex);
            Q2[Q2_size++] = vertex; // Move to Q2 after processing
        }

        Q1_size = Q2_size;
        Q2_size = 0;
    }
}

void readInputFile(const std::string &filename, int &numEdges, int *graph, int &numVertices)
{
    std::ifstream inputFile(filename);
    if (!inputFile.is_open())
    {
        std::cerr << "Failed to open the input file: " << filename << std::endl;
        return;
    }

    int EIndex = numVertices + 1;
    
    int currentNode = -1;
    while (!inputFile.eof())
    {
        int vertex1, vertex2;
        double weight;
        inputFile >> vertex1 >> vertex2 >> weight;
        if (vertex1 != currentNode)
        {
            currentNode = vertex1;
            graph[currentNode] = EIndex;
            std::cout << currentNode << "::" << EIndex << "\n";
        }
        graph[EIndex] = vertex2;

        EIndex++;
    }
    graph[numVertices] = EIndex;
    /*for (int i = 0; i <= numVertices + numEdges; i++)
    {
        std::cout << graph[i] << " " << i << std::endl;
    }*/

    inputFile.close();
}

int main(int argc, char *argv[])
{
    // TODO: Graph size is not right fix it.
    //  std::string filename ="/Users/gayanigupta/Documents/GitHub/Research/graph_decomposition _refined/Example/Tests/core_3.txt";
    int *graph = nullptr;
    int *eliminationOrder = nullptr;
    int numVertices, numEdges;
    // numVertices = 6;
    // numEdges = 8;
    // std::cout << argc << std::endl;
    if (argc < 3)
    {
        std::cerr << "Usage:" << argv[1] << "file, numVertices,numEdges" << std::endl;
    }
    const char *filename = argv[1];
    numVertices = std::atoi(argv[2]);
    numEdges = std::atoi(argv[3]);

    // Create an array to hold the graph data.

    // TODO: Initialize with -1 automatically
    graph = new int[numVertices + numEdges + 4];
    eliminationOrder = new int[numEdges];

    for (int i = 0; i <= numVertices + numEdges; ++i)
    {
        graph[i] = -1;
    }

    // Call the function to read and process the input file.
    readInputFile(filename, numEdges, graph, numVertices);

    // DEBUG :: Printing the graph

    for (int i = 0; i <= numVertices + numEdges; i++)
    {
        std::cout << graph[i] << " " << i << std::endl;
    }
    std::cout << "Number of Vertices: " << numVertices << std::endl;
    std::cout << "Number of Edges: " << numEdges << std::endl;

    maximalChordalSubgraph(numVertices, graph);

    return 0;
}