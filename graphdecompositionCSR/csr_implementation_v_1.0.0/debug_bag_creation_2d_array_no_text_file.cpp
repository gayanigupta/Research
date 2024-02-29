#include <iostream>

class ChordalGraph {
private:
    int vertices;
    int* lowestParent;
    int* rowPtr;
    int* colIndices;
    int* chordalSubgraph;

public:
    ChordalGraph(int v, int edges) : vertices(v) {
        lowestParent = new int[v];
        rowPtr = new int[v + 1]();
        colIndices = new int[edges];
        chordalSubgraph = new int[edges]();
    }

    ~ChordalGraph() {
        delete[] lowestParent;
        delete[] rowPtr;
        delete[] colIndices;
        delete[] chordalSubgraph;
    }

    void addEdge(int v1, int v2, int& edgeCount) {
        if (v1 != v2) {
            colIndices[edgeCount] = v2;
            chordalSubgraph[edgeCount] = 1;
            edgeCount++;
        }
        rowPtr[v1 + 1]++;
    }

    void buildCSR() {
        for (int i = 1; i <= vertices; ++i) {
            rowPtr[i] += rowPtr[i - 1];
        }
    }

    void extractMaximalChordalSubgraph() {
        int* queue1 = new int[vertices];
        int* queue2 = new int[vertices];
        int queue1Size = 0, queue2Size = 0;

        for (int v = 0; v < vertices; ++v) {
            lowestParent[v] = findLowestParent(v);
            if (lowestParent[v] != -1) {
                queue1[queue1Size++] = lowestParent[v];
                chordalSubgraph[v] = 1;
            }
        }

        while (queue1Size > 0) {
            while (queue1Size > 0) {
                int v = queue1[--queue1Size];

                for (int i = rowPtr[v]; i < rowPtr[v + 1]; ++i) {
                    int w = colIndices[i];

                    if (lowestParent[w] == v && isSubset(chordalSubgraph, rowPtr, w)) {
                        chordalSubgraph[i] = 1;
                        queue2[queue2Size++] = w;
                        lowestParent[w] = findNextLowestParent(w);
                    }
                }
            }

            // Swap queues
            int* temp = queue1;
            queue1 = queue2;
            queue2 = temp;
            queue1Size = queue2Size;
            queue2Size = 0;
        }

        delete[] queue1;
        delete[] queue2;
    }

    int findLowestParent(int v) const {
        int lowestParent = -1;
        for (int i = rowPtr[v]; i < rowPtr[v + 1]; ++i) {
            int neighbor = colIndices[i];
            if (lowestParent == -1 || neighbor < lowestParent) {
                lowestParent = neighbor;
            }
        }
        return lowestParent;
    }

    int findNextLowestParent(int v) const {
        int nextLowest = -1;
        for (int i = rowPtr[v]; i < rowPtr[v + 1]; ++i) {
            int neighbor = colIndices[i];
            if (nextLowest == -1 || (neighbor > lowestParent[v] && neighbor < nextLowest)) {
                nextLowest = neighbor;
            }
        }
        return nextLowest;
    }

    bool isSubset(const int* chordalSubgraph, const int* rowPtr, int w) const {
        for (int i = rowPtr[w]; i < rowPtr[w + 1]; ++i) {
            if (!chordalSubgraph[i]) {
                return false;
            }
        }
        return true;
    }

    void printChordalSubgraph() const {
        std::cout << "Chordal Subgraph:\n";
        for (int v = 0; v < vertices; ++v) {
            std::cout << v << ": ";
            for (int i = rowPtr[v]; i < rowPtr[v + 1]; ++i) {
                if (chordalSubgraph[i] == 1) {
                    std::cout << colIndices[i] << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    void findBags() const {
        int** bags = new int*[vertices];
        int* bagsSize = new int[vertices]();

        for (int i = 0; i < vertices; ++i) {
            bags[i] = new int[vertices];
            for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
                int w = colIndices[j];
                if (chordalSubgraph[j] == 1) {
                    bags[i][bagsSize[i]++] = w;
                }
            }
        }

        std::cout << "Bags:\n";
        for (int v = 0; v < vertices; ++v) {
            std::cout << v << ": ";
            for (int i = 0; i < bagsSize[v]; ++i) {
                std::cout << bags[v][i] << " ";
            }
            std::cout << std::endl;
        }

        linkBags(bags, bagsSize);

        // Cleanup
        for (int i = 0; i < vertices; ++i) {
            delete[] bags[i];
        }
        delete[] bags;
        delete[] bagsSize;
    }

    void linkBags(int** bags, const int* bagsSize) const {
        std::cout << "Linked Bags:\n";

        for (int v = 0; v < vertices; ++v) {
            std::cout << v << ": ";
            for (int i = 0; i < bagsSize[v]; ++i) {
                std::cout << bags[v][i] << " ";
            }

            for (int u = 0; u < vertices; ++u) {
                if (u != v) {
                    // Check if there is a common vertex to link bags
                    bool commonVertexExists = false;
                    for (int j = 0; j < bagsSize[u]; ++j) {
                        for (int k = rowPtr[v]; k < rowPtr[v + 1]; ++k) {
                            if (colIndices[k] == bags[u][j]) {
                                commonVertexExists = true;
                                break;
                            }
                        }
                        if (commonVertexExists) {
                            break;
                        }
                    }

                    if (commonVertexExists) {
                        std::cout << "-> " << u << " ";
                    }
                }
            }

            std::cout << std::endl;
        }
    }
};

int main() {
    // Example usage
    int numVertices = 5;
    int numEdges = 5;
    ChordalGraph myChordalGraph(numVertices, numEdges);

    int edgeCount = 0;
    myChordalGraph.addEdge(0, 1, edgeCount);
    myChordalGraph.addEdge(0, 2, edgeCount);
    myChordalGraph.addEdge(1, 3, edgeCount);
    myChordalGraph.addEdge(3, 4, edgeCount);

    myChordalGraph.buildCSR();
    myChordalGraph.extractMaximalChordalSubgraph();
    myChordalGraph.printChordalSubgraph();
    myChordalGraph.findBags();

    return 0;
}
