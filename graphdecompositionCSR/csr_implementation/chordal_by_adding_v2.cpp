#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

const int MAX_NODES = 1000;  // Adjust the maximum number of nodes as needed

void readInputFile(const std::string &filename, int &numEdges, int *graph, int &numVertices);

const int MAX_NODES = 1000;  // Adjust the maximum number of nodes as needed

void readInputFile(const std::string &filename, int &numEdges, int *graph, int &numVertices);

const int MAX_NODES = 1000;  // Adjust the maximum number of nodes as needed

void readInputFile(const std::string &filename, int &numEdges, int *graph, int &numVertices);

// Function to calculate clustering coefficient of a node
float clusteringCoefficient(int adjacencyList[MAX_NODES][MAX_NODES], int node, int nodeDegrees[MAX_NODES]) {
    int edges = 0;
    int neighborsCount = nodeDegrees[node];

    for (int i = 0; i < neighborsCount; ++i) {
        int neighbor = adjacencyList[node][i];

        for (int j = 0; j < neighborsCount; ++j) {
            int anotherNeighbor = adjacencyList[node][j];

            if (neighbor != anotherNeighbor && find(adjacencyList[neighbor], adjacencyList[neighbor] + nodeDegrees[neighbor], anotherNeighbor) != adjacencyList[neighbor] + nodeDegrees[neighbor]) {
                edges++;
            }
        }
    }

    int possibleEdges = neighborsCount * (neighborsCount - 1) / 2;
    return (possibleEdges == 0) ? 0.0 : static_cast<float>(edges) / static_cast<float>(possibleEdges);
}

// Function to make the graph chordal by adding edges
void makeChordal(int adjacencyList[MAX_NODES][MAX_NODES], int nodeDegrees[MAX_NODES], int totalNodes) {
    for (int node = 0; node < totalNodes; ++node) {
        // Calculate clustering coefficient for each node
        float cc = clusteringCoefficient(adjacencyList, node, nodeDegrees);

        // Find a neighbor with the highest clustering coefficient
        int maxCCNeighbor = -1;
        float maxCC = -1.0;

        for (int i = 0; i < nodeDegrees[node]; ++i) {
            int neighbor = adjacencyList[node][i];
            float neighborCC = clusteringCoefficient(adjacencyList, neighbor, nodeDegrees);

            if (neighborCC > maxCC) {
                maxCC = neighborCC;
                maxCCNeighbor = neighbor;
            }
        }

        // Add an edge between the node and its highest clustering coefficient neighbor
        if (maxCCNeighbor != -1) {
            adjacencyList[node][nodeDegrees[node]++] = maxCCNeighbor;
            adjacencyList[maxCCNeighbor][nodeDegrees[maxCCNeighbor]++] = node;
        }
    }
}

int main(int argc, char *argv[]) {
    int adjacencyList[MAX_NODES][MAX_NODES] = {0};  // Adjacency list represented using 2D array
    int nodeDegrees[MAX_NODES] = {0};               // Array to store the degree of each node
    int totalNodes = 0;                             // Number of nodes in the graph

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " filename numVertices numEdges" << std::endl;
        return 1;
    }

    const char *filename = argv[1];
    int numVertices = std::atoi(argv[2]);
    int numEdges = std::atoi(argv[3]);

    // Create an array to hold the graph data.
    int graph[MAX_NODES + MAX_NODES + 4];

    // Initialize the graph array with -1
    std::fill_n(graph, MAX_NODES + MAX_NODES + 4, -1);

    // Call the function to read and process the input file.
    readInputFile(filename, numEdges, graph, numVertices);

    // DEBUG :: Printing the graph
    for (int i = 0; i <= numVertices + numEdges; i++) {
        std::cout << graph[i] << " " << i << std::endl;
    }

    // Print the number of vertices and edges.
    std::cout << "Number of Vertices: " << numVertices << std::endl;
    std::cout << "Number of Edges: " << numEdges << std::endl;

    std::cout << "size of graph: " << sizeof(graph) / sizeof(int) << std::endl;

    // Make the graph chordal by adding edges
    makeChordal(adjacencyList, nodeDegrees, totalNodes);

    // Print the final chordal graph
    for (int node = 0; node < totalNodes; ++node) {
        std::cout << "Node " << node << " connects to: ";
        for (int i = 0; i < nodeDegrees[node]; ++i) {
            std::cout << adjacencyList[node][i] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}

void readInputFile(const std::string &filename, int &numEdges, int *graph, int &numVertices) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the input file: " << filename << std::endl;
        return;
    }

    int EIndex = numVertices;
    int currentNode = -1;

    // Read each line of the input file.
    while (!inputFile.eof()) {
        int vertex1, vertex2;
        double weight;
        inputFile >> vertex1 >> vertex2 >> weight;

        // If vertex1 is not equal to currentNode, update currentNode to vertex1.
        if (vertex1 != currentNode) {
            currentNode = vertex1;
            graph[currentNode] = EIndex;
            std::cout << currentNode << "::" << EIndex << "\n";
        }
        graph[EIndex] = vertex2;

        EIndex++;
    }

    for (int i = 0; i <= numVertices + numEdges; i++) {
        std::cout << graph[i] << " " << i << std::endl;
    }

    inputFile.close();
}
