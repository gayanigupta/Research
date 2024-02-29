#include <iostream>
#include <vector>
#include <queue>
#include <unordered_set>

class ChordalGraph {
private:
    int vertices;
    int* lowestParent;
    int* chordalNeighbors;
    int* chordalSubgraph;

public:
    ChordalGraph(int v) : vertices(v) {
        lowestParent = new int[v];
        chordalNeighbors = new int[v * v]();  // Assuming a maximum of v * v neighbors
        chordalSubgraph = new int[v * v]();  // Assuming a maximum of v * v edges in chordal subgraph
    }

    ~ChordalGraph() {
        delete[] lowestParent;
        delete[] chordalNeighbors;
        delete[] chordalSubgraph;
    }

    void addEdge(int v1, int v2) {
        chordalNeighbors[v1 * vertices + v2] = 1;
        chordalNeighbors[v2 * vertices + v1] = 1;  // Assuming an undirected graph
    }

    void extractMaximalChordalSubgraph() {
        std::queue<int> queue1, queue2;

        for (int v = 0; v < vertices; ++v) {
            lowestParent[v] = findLowestParent(v);
            if (lowestParent[v] != -1) {
                queue1.push(lowestParent[v]);
                chordalSubgraph[v * vertices + lowestParent[v]] = 1;
                chordalSubgraph[lowestParent[v] * vertices + v] = 1;
            }
        }

        while (!queue1.empty()) {
            while (!queue1.empty()) {
                int v = queue1.front();
                queue1.pop();

                for (int w = 0; w < vertices; ++w) {
                    if (chordalNeighbors[v * vertices + w] && lowestParent[w] == v &&
                        isSubset(chordalSubgraph + v * vertices, chordalSubgraph + w * vertices)) {
                        chordalSubgraph[w * vertices + v] = 1;
                        chordalSubgraph[v * vertices + w] = 1;
                        queue2.push(w);
                        lowestParent[w] = findNextLowestParent(w);
                    }
                }
            }

            std::swap(queue1, queue2);
        }
    }

    int findLowestParent(int v) const {
        int lowestParent = -1;
        for (int neighbor = 0; neighbor < vertices; ++neighbor) {
            if (chordalNeighbors[v * vertices + neighbor] == 1) {
                if (lowestParent == -1 || neighbor < lowestParent) {
                    lowestParent = neighbor;
                }
            }
        }
        return lowestParent;
    }

    int findNextLowestParent(int v) const {
        int nextLowest = -1;
        for (int neighbor = 0; neighbor < vertices; ++neighbor) {
            if (chordalNeighbors[v * vertices + neighbor] == 1) {
                if (nextLowest == -1 || (neighbor > lowestParent[v] && neighbor < nextLowest)) {
                    nextLowest = neighbor;
                }
            }
        }
        return nextLowest;
    }

    bool isSubset(const int* set1, const int* set2) const {
        for (int i = 0; i < vertices; ++i) {
            if (set1[i] && !set2[i]) {
                return false;
            }
        }
        return true;
    }

    void printChordalSubgraph() const {
        std::cout << "Chordal Subgraph:\n";
        for (int v = 0; v < vertices; ++v) {
            std::cout << v << ": ";
            for (int w = 0; w < vertices; ++w) {
                if (chordalSubgraph[v * vertices + w] == 1) {
                    std::cout << w << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    void findBags() const {
        std::vector<std::unordered_set<int>> bags(vertices);

        for (int v = 0; v < vertices; ++v) {
            for (int w = 0; w < vertices; ++w) {
                if (chordalSubgraph[v * vertices + w] == 1) {
                    bags[v].insert(w);
                }
            }
        }

        std::cout << "Bags:\n";
        for (int v = 0; v < vertices; ++v) {
            std::cout << v << ": ";
            for (int w : bags[v]) {
                std::cout << w << " ";
            }
            std::cout << std::endl;
        }

        linkBags(bags);
    }

    void linkBags(const std::vector<std::unordered_set<int>>& bags) const {
        std::cout << "Linked Bags:\n";

        for (int v = 0; v < vertices; ++v) {
            std::cout << v << ": ";
            for (int w : bags[v]) {
                std::cout << w << " ";
            }

            for (int u = 0; u < vertices; ++u) {
                if (u != v && !bags[v].count(u)) {
                    // Check if there is a common vertex to link bags
                    bool commonVertexExists = false;
                    for (int x : bags[u]) {
                        if (chordalSubgraph[v * vertices + x] == 1) {
                            commonVertexExists = true;
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
    ChordalGraph myChordalGraph(numVertices);

    myChordalGraph.addEdge(0, 1);
    myChordalGraph.addEdge(0, 2);
    myChordalGraph.addEdge(1, 3);
    myChordalGraph.addEdge(3, 4);

    myChordalGraph.extractMaximalChordalSubgraph();
    myChordalGraph.printChordalSubgraph();
    myChordalGraph.findBags();

    return 0;
}
