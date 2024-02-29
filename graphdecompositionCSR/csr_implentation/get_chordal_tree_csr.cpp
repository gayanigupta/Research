#include <iostream>
#include <fstream>
#include <queue>

class ChordalGraph {
private:
    int vertices;
    int* rowPtr;
    int* colIndex;
    int* eliminationOrder;
    int* bags;
    int* bagParents;

public:
    ChordalGraph(int v) : vertices(v) {
        rowPtr = new int[v + 1]();
        bagParents = new int[v]();
    }

    ~ChordalGraph() {
        delete[] rowPtr;
        delete[] colIndex;
        delete[] eliminationOrder;
        delete[] bags;
        delete[] bagParents;
    }

    void addEdge(int v1, int v2) {
        rowPtr[v1 + 1]++;
    }

    void buildCSR() {
        for (int i = 1; i <= vertices; ++i) {
            rowPtr[i] += rowPtr[i - 1];
        }
        colIndex = new int[rowPtr[vertices]];
    }

void extractMaximalChordalSubgraph() {
    // Create an array to mark vertices as eliminated or not
    bool* eliminated = new bool[vertices]();

    // Assuming a greedy elimination order
    for (int i = 0; i < vertices; ++i) {
        int v = eliminationOrder[i];

        if (!eliminated[v]) {
            // Mark v and its non-neighbor vertices as eliminated
            eliminated[v] = true;

            for (int w = 0; w < vertices; ++w) {
                if (!eliminated[w] && rowPtr[w] > 0 && colIndex[rowPtr[w] - 1] == v) {
                    eliminated[w] = true;
                }
            }
        }
    }

    // Create a new CSR representation for the chordal subgraph
    int* newColIndex = new int[vertices * vertices]();
    int* newRowPtr = new int[vertices + 1]();

    int newIndex = 0;
    for (int i = 0; i < vertices; ++i) {
        if (!eliminated[i]) {
            for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
                int w = colIndex[j];
                if (!eliminated[w]) {
                    newColIndex[newIndex++] = w;
                }
            }
            newRowPtr[i + 1] = newIndex;
        }
    }

    // Cleanup old CSR representation
    delete[] colIndex;
    delete[] rowPtr;

    // Assign the new CSR representation
    colIndex = newColIndex;
    rowPtr = newRowPtr;

    delete[] eliminated;
}
    void findBags() {
        bags = new int[vertices * vertices]();
        int* bagsSize = new int[vertices]();

        for (int v = 0; v < vertices; ++v) {
            int bagSize = 1; // Initialize with the vertex itself
            bags[v * vertices] = v;

            for (int i = rowPtr[v]; i < rowPtr[v + 1]; ++i) {
                int w = colIndex[i];

                if (bagParents[w] == -1) {
                    // If w does not have a parent bag, create a new bag
                    for (int j = 0; j < bagSize; ++j) {
                        bags[v * vertices + bagSize + j] = bags[w * vertices + j];
                        bagParents[w] = v;
                    }
                    bagSize += bagSize;
                } else {
                    // If w has a parent bag, link the bags using the common vertex
                    int commonVertex = bags[v * vertices]; // Assuming the first vertex in each bag is the common vertex
                    for (int j = 0; j < bagSize; ++j) {
                        bags[bagParents[w] * vertices + bagSize + j] = bags[v * vertices + j];
                    }
                    bagParents[v] = bagParents[w];
                    bagSize += bagSize;
                }
            }

            bagsSize[v] = bagSize;
        }

        // Update bags to trim excess space
        int* trimmedBags = new int[vertices * vertices]();
        for (int v = 0; v < vertices; ++v) {
            for (int j = 0; j < bagsSize[v]; ++j) {
                trimmedBags[v * vertices + j] = bags[v * vertices + j];
            }
        }

        delete[] bags;
        bags = trimmedBags;

        delete[] bagsSize;
    }

    void printBags() const {
        std::cout << "Bags:\n";
        for (int v = 0; v < vertices; ++v) {
            std::cout << v << ": ";
            for (int j = 0; j < vertices && bags[v * vertices + j] != -1; ++j) {
                std::cout << bags[v * vertices + j] << " ";
            }
            std::cout << std::endl;
        }
    }
};

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <file> <numVertices> <numEdges>" << std::endl;
        return 1; // Exit with an error code
    }

    const char *filename = argv[1];
    int numVertices = std::atoi(argv[2]);
    int numEdges = std::atoi(argv[3]);

    // Create ChordalGraph instance
    ChordalGraph myChordalGraph(numVertices);

    // Read data from the file and populate the graph
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return 1; // Exit with an error code
    }

    for (int i = 0; i < numEdges; ++i) {
        int v1, v2;
        inputFile >> v1 >> v2;
        // Assuming an undirected graph
        myChordalGraph.addEdge(v1, v2);
        myChordalGraph.addEdge(v2, v1);
    }

    // Close the input file
    inputFile.close();

    // Build CSR representation
    myChordalGraph.buildCSR();

    // Extract maximal chordal subgraph
    myChordalGraph.extractMaximalChordalSubgraph();

    // Find bags and link them
    myChordalGraph.findBags();

    // Print the bags
    myChordalGraph.printBags();

    return 0;
}
