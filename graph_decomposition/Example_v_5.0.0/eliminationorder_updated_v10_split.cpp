#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <set>
#include <unordered_map>
#include "graph_core.hpp"
#include "graph_elimination.hpp"

// INPUT HEADERS
#include "translate_from_input.hpp"
#include "input_to_network.hpp"
#include "structure_defs.hpp"


// OUTPUT HEADERS
#include "printout_network.hpp"
#include "printout_others.hpp"
using namespace std;

// Function to add an edge to the graph


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
    // Populate the edges vector
     ifstream inputFile(argv[3]);
    if (!inputFile) {
        cerr << "Error opening the file." << endl;
        return 1;
    }

    int n = 0; // Number of vertices in the graph
    vector<vector<int>> graph;
    set<int> verticesSet;

    string line;
    while (getline(inputFile, line)) {
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
    inputFile.clear();
    inputFile.seekg(0);

    while (getline(inputFile, line)) {
        int u, v;
        istringstream iss(line);
        if (!(iss >> u >> v)) {
            continue;
        }

        int reindexedU = reindexMap[u];
        int reindexedV = reindexMap[v];

        addEdge(graph, reindexedU, reindexedV);
    }

    inputFile.close();
    cout<<graph.size();
    // Calculate the perfect elimination order using the Lex-M algorithm
    vector<int> eliminationOrder = findPerfectEliminationOrder(graph);

    // Find the chordal edges using the elimination order
    findChordalEdgesWithEliminationOrder(graph, eliminationOrder);


    for (int i = 0; i < eliminationOrder.size(); i++)
        cout << eliminationOrder[i] << endl;

    return 0;
}
