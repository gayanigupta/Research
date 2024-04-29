#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

const int MAX_NODES = 1000;  // Adjust the maximum number of nodes as needed

void readInputFile(const string &filename, int &numEdges, int *graph, int &numVertices);
float clusteringCoefficient(int *adjacencyList, int *nodeDegrees, int node, int *eliminationOrder, int &count);
void makeChordal(int *adjacencyList, int *nodeDegrees, int totalNodes, int *eliminationOrder);

int main(int argc, char *argv[]) {
    int adjacencyList[MAX_NODES * MAX_NODES] = {0};  // Adjacency list represented using 1D array
    int nodeDegrees[MAX_NODES] = {0};               // Array to store the degree of each node
    int eliminationOrder[MAX_NODES];                 // Array to store the elimination order of nodes
    int totalNodes = 0;                             // Number of nodes in the graph

    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " filename numVertices numEdges" << endl;
        return 1;
    }

    const char *filename = argv[1];
    int numVertices = atoi(argv[2]);
    int numEdges = atoi(argv[3]);

    // Initialize totalNodes
    totalNodes = numVertices;

    // Create an array to hold the graph data.
    int graph[MAX_NODES + MAX_NODES + 4];

    // Initialize the graph array with -1
    fill_n(graph, MAX_NODES + MAX_NODES + 4, -1);

    // Call the function to read and process the input file.
    readInputFile(filename, numEdges, graph, numVertices);

    // Make the graph chordal by adding edges
    makeChordal(adjacencyList, nodeDegrees, totalNodes, eliminationOrder);

    // Print the final chordal graph
    for (int node = 0; node < totalNodes; ++node) {
        cout << "Node " << node << " connects to: ";
        for (int i = 0; i < nodeDegrees[node]; ++i) {
            cout << adjacencyList[node * MAX_NODES + i] << " ";
        }
        cout << endl;
    }

    // Print the elimination order
    cout << "Elimination Order: ";
    for (int i = 0; i < totalNodes; ++i) {
        cout << eliminationOrder[i] << " ";
    }
    cout << endl;

    return 0;
}

void readInputFile(const string &filename, int &numEdges, int *graph, int &numVertices) {

    ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        cerr << "Failed to open the input file: " << filename << endl;
        return;
    }

    int EIndex = numVertices;
    cout<<"EIndex: "<<EIndex<<endl;

    // Read each line of the input file.
    while (inputFile >> graph[EIndex] >> graph[EIndex + 1] >> graph[EIndex + 2]) {
        EIndex += 3;
    }

    inputFile.close();
}

float clusteringCoefficient(int *adjacencyList, int *nodeDegrees, int node, int *eliminationOrder, int &count) {
    int edges = 0; // Counter to store the number of edges between neighbors
    int neighborsCount = nodeDegrees[node]; // Get the number of neighbors for the given node

    // Iterate through each pair of neighbors
    for (int i = 0; i < neighborsCount; ++i) {
        int neighbor = adjacencyList[node * MAX_NODES + i]; // Get the current neighbor

        // Start the inner loop from i + 1 to iterate over the upper triangle of the adjacency list
        for (int j = i + 1; j < neighborsCount; ++j) {
            int anotherNeighbor = adjacencyList[node * MAX_NODES + j]; // Get another neighbor

            // Check if there is an edge between the current pair of neighbors
            if (find(adjacencyList + neighbor * MAX_NODES, adjacencyList + neighbor * MAX_NODES + nodeDegrees[neighbor], anotherNeighbor) != adjacencyList + neighbor * MAX_NODES + nodeDegrees[neighbor]) {
                edges++; // Increment the edge count
            }
        }
    }

    // Calculate the maximum possible number of edges between neighbors
    int possibleEdges = neighborsCount * (neighborsCount - 1) / 2;

    // Calculate the clustering coefficient using the formula: (number of edges between neighbors * 2) / (maximum possible number of edges between neighbors)
    // Multiply by 2 because each edge between neighbors is counted twice (once for each neighbor)
    float cc = (possibleEdges == 0) ? 0.0f : static_cast<float>(edges * 2) / static_cast<float>(possibleEdges);

    // Update the clustering coefficient in the elimination order array
    eliminationOrder[count++] = node;

    return cc; // Return the clustering coefficient
}


void makeChordal(int *adjacencyList, int *nodeDegrees, int totalNodes, int *eliminationOrder) {
    int count = 0; // Counter for elimination order

    for (int node = 0; node < totalNodes; ++node) {
        // Calculate clustering coefficient for each node
        float cc = clusteringCoefficient(adjacencyList, nodeDegrees, node, eliminationOrder, count);

        // Find a neighbor with the highest clustering coefficient
        int maxCCNeighbor = -1;
        float maxCC = -1.0f;

        for (int i = 0; i < nodeDegrees[node]; ++i) {
            int neighbor = adjacencyList[node * MAX_NODES + i];
            float neighborCC = clusteringCoefficient(adjacencyList, nodeDegrees, neighbor, eliminationOrder, count);

            if (neighborCC > maxCC) {
                maxCC = neighborCC;
                maxCCNeighbor = neighbor;
            }
        }

        // Add an edge between the node and its highest clustering coefficient neighbor
        if (maxCCNeighbor != -1) {
            adjacencyList[node * MAX_NODES + nodeDegrees[node]++] = maxCCNeighbor;
            adjacencyList[maxCCNeighbor * MAX_NODES + nodeDegrees[maxCCNeighbor]++] = node;
        }
    }
}
