#include <iostream>
#include <fstream>
#include <vector>
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
    ChordalGraph(int v) : vertices(v), rowPtr(new int[v + 1]()), colIndex(nullptr), eliminationOrder(nullptr), bags(nullptr), bagParents(new int[v]()) {}

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
        colIndex = new int[rowPtr[vertices]]();
    }

    void extractMaximalChordalSubgraph() {
        eliminationOrder = new int[vertices];
        for (int i = 0; i < vertices; ++i) {
            int maxDegree = -1;
            int maxVertex = -1;

            for (int v = 0; v < vertices; ++v) {
                if (rowPtr[v + 1] - rowPtr[v] > maxDegree) {
                    maxDegree = rowPtr[v + 1] - rowPtr[v];
                    maxVertex = v;
                }
            }

            eliminationOrder[i] = maxVertex;
            rowPtr[maxVertex + 1] = rowPtr[maxVertex];

            std::cout << "Elimination Order[" << i << "]: " << maxVertex << std::endl;
        }

        bool* eliminated = new bool[vertices]();
        int* newColIndex = new int[rowPtr[vertices]];
        int* newRowPtr = new int[vertices + 1]();

        int newIndex = 0;
        for (int i = 0; i < vertices; ++i) {
            int v = eliminationOrder[i];

            if (!eliminated[v]) {
                std::cout << "Processing vertex: " << v << std::endl;

                for (int j = rowPtr[v]; j < rowPtr[v + 1]; ++j) {
                    int w = colIndex[j];
                    if (!eliminated[w]) {
                        newColIndex[newIndex++] = w;
                        std::cout << "Adding neighbor: " << w << std::endl;
                    }
                }
                newRowPtr[v + 1] = newIndex;

                eliminated[v] = true;
                for (int w = 0; w < vertices; ++w) {
                    if (!eliminated[w] && rowPtr[w] > 0 && colIndex[rowPtr[w] - 1] == v) {
                        eliminated[w] = true;
                        std::cout << "Marking as eliminated: " << w << std::endl;
                    }
                }
            }
        }

        delete[] colIndex;
        colIndex = newColIndex;
        rowPtr = newRowPtr;

        delete[] eliminated;
        delete[] eliminationOrder;
    }

    void findBags() {
        bags = new int[vertices * vertices]();
        int* bagsSize = new int[vertices]();
        int* visited = new int[vertices]();

        for (int v = 0; v < vertices; ++v) {
            if (visited[v] == 0) {
                int bagSize = 0;
                int current = v;

                while (current != -1 && visited[current] == 0) {
                    visited[current] = 1;
                    bags[v * vertices + bagSize++] = current;

                    int nextNeighbor = -1;
                    for (int i = rowPtr[current]; i < rowPtr[current + 1]; ++i) {
                        int w = colIndex[i];
                        if (visited[w] == 0 && bagParents[w] == current) {
                            nextNeighbor = w;
                            break;
                        }
                    }

                    current = nextNeighbor;
                }

                bagParents[current] = v;
                bagsSize[v] = bagSize;

                for (int i = 0; i < bagSize; ++i) {
                    visited[bags[v * vertices + i]] = 0;
                }
            }
        }

        int* trimmedBags = new int[vertices * vertices]();
        for (int v = 0; v < vertices; ++v) {
            for (int j = 0; j < bagsSize[v]; ++j) {
                trimmedBags[v * vertices + j] = bags[v * vertices + j];
            }
        }

        delete[] bags;
        bags = trimmedBags;

        delete[] bagsSize;
        delete[] visited;
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
        return 1;
    }

    const char *filename = argv[1];
    int numVertices = std::atoi(argv[2]);
    int numEdges = std::atoi(argv[3]);

    ChordalGraph myChordalGraph(numVertices);

    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return 1;
    }

    for (int i = 0; i < numEdges; ++i) {
        int v1, v2;
        inputFile >> v1 >> v2;
        myChordalGraph.addEdge(v1, v2);
        myChordalGraph.addEdge(v2, v1);
    }

    inputFile.close();

    myChordalGraph.buildCSR();

    myChordalGraph.extractMaximalChordalSubgraph();

    myChordalGraph.findBags();

    myChordalGraph.printBags();

    return 0;
}
