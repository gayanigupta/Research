#pragma once
/* This library take a files as  input, reads all the vertices and outputs the graph*/
#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

using namespace std;

// Function to add an edge to the graph
void addEdge(vector<vector<int>>& graph, int u, int v) {
    graph[u].push_back(v);
    graph[v].push_back(u);
}

bool BuildGraph(string strFileName, vector<vector<int>>& graph)
{
    ifstream inputFile(strFileName);
    if (!inputFile) {
        cerr << "Error opening the file." << endl;
        return false;
    }

    int n = 0; // Number of vertices in the graph
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
    return true;
}