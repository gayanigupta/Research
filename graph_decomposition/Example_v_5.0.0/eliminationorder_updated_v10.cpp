#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

// INPUT HEADERS
#include "translate_from_input.hpp"
#include "input_to_network.hpp"
#include "structure_defs.hpp"


// OUTPUT HEADERS
#include "printout_network.hpp"
#include "printout_others.hpp"

using namespace std;

// Function to add an edge to the graph
void addEdge(vector<vector<int>>& graph, int u, int v) {
    graph[u].push_back(v);
    graph[v].push_back(u);
}

// Function to calculate the lowest parent of a vertex in the graph
int findLowestParent(const vector<vector<int>>& graph, const vector<bool>& eliminated, int vertex) {
    int lowestParent = vertex;

    for (int neighbor : graph[vertex]) {
        if (!eliminated[neighbor] && neighbor < lowestParent) {
            lowestParent = neighbor;
        }
    }

    return lowestParent;
}

// Function to find the perfect elimination order using the Lex-M algorithm
vector<int> findPerfectEliminationOrder(const vector<vector<int>>& graph) {
    int n = graph.size();
    vector<bool> eliminated(n, false);
    vector<int> eliminationOrder;

    for (int i = 0; i < n; ++i) {
        int lowestParent = -1;

        // Find the lowest parent that has not been eliminated
        for (int v = 0; v < n; ++v) {
            if (!eliminated[v]) {
                int parent = findLowestParent(graph, eliminated, v);
                if (lowestParent == -1 || parent < lowestParent) {
                    lowestParent = parent;
                }
            }
        }

        // Check if removing the lowest parent makes the remaining graph a clique
        bool isClique = true;
        for (int v = 0; v < n; ++v) {
            if (!eliminated[v] && v != lowestParent) {
                if (find(graph[v].begin(), graph[v].end(), lowestParent) == graph[v].end()) {
                    isClique = false;
                    break;
                }
            }
        }

        if (!isClique) {
            cerr << "The input graph is not chordal. Please ensure it is a chordal graph before using this algorithm." << endl;
            exit(1);
        }

        eliminationOrder.push_back(lowestParent);
        eliminated[lowestParent] = true;
    }

    return eliminationOrder;
}

// Function to find the chordal edges with a given elimination order
void findChordalEdgesWithEliminationOrder(const vector<vector<int>>& graph, const vector<int>& eliminationOrder) {
    int n = graph.size();
    vector<unordered_set<int>> neighbors(n);

    // Store neighbors of each vertex
    for (int u = 0; u < n; ++u) {
        for (int v : graph[u]) {
            neighbors[u].insert(v);
        }
    }

    for (int i = 0; i < n; ++i) {
        int u = eliminationOrder[i];
        for (int j = i + 1; j < n; ++j) {
            int v = eliminationOrder[j];

            // Check if v is a neighbor of u
            if (neighbors[u].count(v) == 1) {
                // Check if v is adjacent to all the previous vertices in the elimination order
                bool isChordalEdge = true;
                for (int k = i + 1; k < j; ++k) {
                    int w = eliminationOrder[k];
                    if (neighbors[v].count(w) == 0) {
                        isChordalEdge = false;
                        break;
                    }
                }

                if (isChordalEdge) {
                    cout << "Chordal Edge: (" << u << ", " << v << ")" << endl;
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    clock_t q, q1, q2, t;
    vector<Edge> edges;

    q = clock();

    if (argc < 3)
    {
        cout << "INPUT ERROR:: At least 2 inputs required. First: filename \n Second: Filetypes: 1:node_node_wt 2:node_wt_node 3:node_node 4:node_node (Only option 1 is active now) \n Third. Name of new file \n Fourth. Name of Map file\n";
        return 0;
    }

    ifstream the_file(argv[1]);
    if (!the_file.is_open())
    {
        cout << "INPUT ERROR:: Could not open file\n";
        return 0;
    }

    A_Network x1;
    int nodes = -1;
    map_int_st revmap;
    int type = atoi("1");
    translate_input(argv[1], type, argv[3], argv[4]);

    q = clock() - q;
    cout << "Total Time for Preprocessing: " << ((float)q) / CLOCKS_PER_SEC << "\n";

    q = clock();
    readin_network(&x1, argv[3], nodes);
    nodes = x1.size();
    create_map(argv[4], &revmap);
    q = clock() - q;
    cout << "Total Time for Reading Network: " << ((float)q) / CLOCKS_PER_SEC << "\n";

    // Populate the edges vector
    ifstream dataFile(argv[3], ios::in); // Open input file
    string line;
    stringstream linestream;
    int vertex_N = 0;
    while (getline(dataFile, line))
    {
        linestream.clear();
        linestream << line;
        int node1, node2;
        double weight;
        linestream >> node1 >> node2 >> weight;
        vertex_N = max(vertex_N, max(node1, node2)) ;
        Edge edge;
        edge.node1 = node1;
        edge.node2 = node2;
        edge.edge_wt = weight;
        if (edge.node1 < edge.node2){
        edges.push_back(edge);
        }
    }
    //dataFile.close();

    //print(edges);

    vertex_N = vertex_N + 1;
    int n = 0; // Number of vertices in the graph
    vector<vector<int>> graph;
    set<int> verticesSet;

    while (getline(dataFile, line)) {
        int u, v;
        istringstream iss(line);
        if (!(iss >> u >> v)) {
            cerr << "Error reading the edge from line: " << line << endl;
            continue;
        }

        verticesSet.insert(u);
        verticesSet.insert(v);
    }

    n = verticesSet.size();
    graph.resize(n);

    // Reindex vertices to avoid gaps in the vertex indices
    unordered_map<int, int> reindexMap;
    int newIndex = 0;
    for (int vertex : verticesSet) {
        reindexMap[vertex] = newIndex;
        newIndex++;
    }

    // Rewind the input file and re-read edges, now adding them to the graph
    dataFile.clear();
    dataFile.seekg(0);

    while (getline(dataFile, line)) {
        int u, v;
        istringstream iss(line);
        if (!(iss >> u >> v)) {
            continue;
        }

        int reindexedU = reindexMap[u];
        int reindexedV = reindexMap[v];

        addEdge(graph, reindexedU, reindexedV);
    }

    dataFile.close();

    // Calculate the perfect elimination order using the Lex-M algorithm
    vector<int> eliminationOrder = findPerfectEliminationOrder(graph);

    // Find the chordal edges using the elimination order
    findChordalEdgesWithEliminationOrder(graph, eliminationOrder);

    for (int i = 0; i < eliminationOrder.size(); i++)
        cout << eliminationOrder[i] << endl;

    return 0;
}
