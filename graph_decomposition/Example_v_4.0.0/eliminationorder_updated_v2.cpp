#include <iostream>
#include <list>
#include <set>

// INPUT HEADERS
#include "translate_from_input.hpp"
#include "input_to_network.hpp"
#include "structure_defs.hpp"

// OUTPUT HEADERS
#include "printout_network.hpp"
#include "printout_others.hpp"

using namespace std;

bool isSubset(const std::set<int>& set1, const std::set<int>& set2);
std::pair<std::vector<Edge>, std::vector<int> > findChordalEdgesWithEliminationOrder(const std::vector<Edge>& edges, int numVertices);
/**
 * - Description: Checks if one set is a subset of another set.
   - Parameters: `set1` and `set2` are the sets to be compared.
   - Return Type: `bool`
   - Returns: `true` if `set1` is a subset of `set2`, `false` otherwise.

*/
bool isSubset(const std::set<int>& set1, const std::set<int>& set2)
{
    std::cout << "Checking if set1 is a subset of set2..." << std::endl;
    std::cout << "set1: ";
    for (int element : set1)
    {
        std::cout << element << " ";
    }
    std::cout << std::endl;
    
    std::cout << "set2: ";
    for (int element : set2)
    {
        std::cout << element << " ";
    }
    std::cout << std::endl;
    
    for (int element : set1)
    {
        if (set2.find(element) == set2.end())
        {
            return false;
        }
    }
    
    return true;
}


/**
 *    - Description: Finds the chordal edges and elimination order of a graph.
   - Parameters:
     - `edges`: A vector of edges in the graph.
     - `numVertices`: The number of vertices in the graph.
   - Return Type: `std::pair<std::vector<Edge>, std::vector<int>>`
   - Returns: A pair containing the chordal edges and elimination order.
   - Steps:
     1. Initialize variables for chordal edges (`chordalEdges`), elimination order (`eliminationOrder`), lowest parent (`LP`), sets (`C`), processed vertices (`processed`), and Q1 array (`Q1Array`).
     2. Find the lowest parent (`LP`) and populate the Q1 array.
     3. Process the vertices until the Q1 array is empty:
        - For each vertex in the Q1 array, check its neighbors and update the sets (`C`), chordal edges (`chordalEdges`), and processed vertices (`processed`) accordingly.
     4. Find the nodes without parents and eliminate chordal nodes first by adding them to the elimination order.
     5. Add the remaining non-chordal nodes to the elimination order.
     6. Return the pair containing the chordal edges and elimination order.

*/


std::pair<std::vector<Edge>, std::vector<int> > findChordalEdgesWithEliminationOrder(const std::vector<Edge>& edges, int numVertices)
{
    std::vector<Edge> chordalEdges;
    std::vector<int> eliminationOrder;

    std::vector<int> LP(numVertices, -1);
    std::vector<std::set<int> > C(numVertices, std::set<int>());

    std::vector<bool> processed(numVertices, false);
    std::vector<bool> Q1Array(numVertices, false);

    // Step 1: Initialize LP, C, processed, and Q1Array based on input edges
    for (const Edge& edge : edges)
    {
        int v = edge.node2;
        int w = edge.node1;
        double wt = edge.edge_wt;

        if (wt != 0 && LP[v] == -1)
        {
            LP[v] = w;
            processed[v] = false;  // Make all vertices initially unprocessed

            if (!Q1Array[w])
            {
                Q1Array[w] = true;
                C[v].clear();
            }
        }
    }

    // Step 2: Find nodes without parents (nodes not present in LP)
    std::vector<int> nodesWithoutParents;
    for (int v = 0; v < numVertices; ++v)
    {
        if (LP[v] == -1)
        {
            nodesWithoutParents.push_back(v);
        }
    }

    // Step 3: Process nodes without parents and eliminate chordal nodes
    for (int v : nodesWithoutParents)
    {
        eliminationOrder.push_back(v);
        processed[v] = true;

        // Check if the node has children and update their chordal sets
        for (const Edge& edge : edges)
        {
            int w = edge.node2;
            int u = edge.node1;
            double wt = edge.edge_wt;

            if (LP[w] == v && std::includes(C[v].begin(), C[v].end(), C[w].begin(), C[w].end()) && !processed[w])
            {
                C[v].insert(w);  // Insert to child's chordal set
                Edge chordalEdge = {v, w, wt};  // Edge should be v -> w
                chordalEdges.push_back(chordalEdge);
                Q1Array[w] = true;
                processed[w] = true;  // Set w as processed

                // Print additional information
                std::cout << "Found chordal edge: " << chordalEdge.node1 << " - " << chordalEdge.node2 << std::endl;
                std::cout << "Updated C[" << v << "]: ";
                for (int node : C[v])
                {
                    std::cout << node << " ";
                }
                std::cout << std::endl;
            }
        }
    }

    bool allProcessed = false;

    // Step 4: Iterate until all vertices are processed
    while (!allProcessed)
    {
        allProcessed = true;

        for (int v = 0; v < numVertices; ++v)
        {
            if (Q1Array[v] && !processed[v])
            {
                allProcessed = false;

                // Step 5: Check edges for chordal properties
                for (const Edge& edge : edges)
                {
                    int w = edge.node2;
                    int u = edge.node1;
                    double wt = edge.edge_wt;

                    if (LP[w] == v && std::includes(C[v].begin(), C[v].end(), C[w].begin(), C[w].end()) && !processed[w])
                    {
                        C[v].insert(w);  // Insert to child's chordal set
                        Edge chordalEdge = {v, w, wt};  // Edge should be v -> w
                        chordalEdges.push_back(chordalEdge);
                        Q1Array[w] = true;
                        processed[v] = true; 
                        processed[w] = true;  // Set w as processed

                        // Print additional information
                        std::cout << "Found chordal edge: " << chordalEdge.node1 << " - " << chordalEdge.node2 << std::endl;
                        std::cout << "Updated C[" << v << "]: ";
                        for (int node : C[v])
                        {
                            std::cout << node << " ";
                        }
                        std::cout << std::endl;
                    }
                }

                Q1Array[v] = false;
            }
        }
    }

    // Print variable contents for debugging purposes
    std::cout << "LP: ";
    for (int value : LP)
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    std::cout << "C:" << std::endl;
    for (int i = 0; i < numVertices; ++i)
    {
        std::cout << "C[" << i << "]: ";
        for (int node : C[i])
        {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Processed: ";
    for (bool value : processed)
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    std::cout << "Q1Array: ";
    for (bool value : Q1Array)
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    // Return the chordal edges and elimination order
    return std::make_pair(chordalEdges, eliminationOrder);
}


/**
 *   The main driver program to test the methods of the graph class.
   - Parameters: `argc` and `argv` are command-line arguments.
   - Return Type: `int`
   - Returns: Exit status of the program.
   - Steps:
     1. Perform input validation and check if valid input is given.
     2. Open the input file and check if it was successfully opened.
     3. Initialize variables for clock times (`q`, `q1`, `q2`, `t`), vector of edges (`edges`), network (`x1`), number of nodes (`nodes`), reverse map (`revmap`), and type of file (`type`).
     4. Perform preprocessing steps, including file translation and creating a reverse map.
     5. Read the network from the translated file and create the reverse map.
     6. Find the chordal edges with elimination order using the `findChordalEdgesWithEliminationOrder` method.
     7. Output the chordal edges and elimination order to the console.
     8. Return the exit status of the program.

*/
int main(int argc, char* argv[])
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
    ifstream dataFile(argv[1], ios::in); // Open input file
    string line;
    stringstream linestream;
    while (getline(dataFile, line))
    {
        linestream.clear();
        linestream << line;
        int node1, node2;
        double weight;
        linestream >> node1 >> node2 >> weight;
        Edge edge;
        edge.node1 = node1;
        edge.node2 = node2;
        edge.edge_wt = weight;
        edges.push_back(edge);
    }
    dataFile.close();

    std::pair<std::vector<Edge>, std::vector<int> > result = findChordalEdgesWithEliminationOrder(edges, nodes);
    std::vector<Edge> chordalEdges = result.first;
    std::vector<int> eliminationOrder = result.second;

    cout << "Chordal Edges: " << endl;
    if (chordalEdges.empty()) {
        cout << "No chordal edges found." << endl;
    } else {
        for (const Edge& edge : chordalEdges)
        {
            cout << edge.node1 << " - " << edge.node2 << endl;
        }
    }

    cout << "Elimination Order: " << endl;
    if (eliminationOrder.empty()) {
        cout << "No elimination order found." << endl;
    } else {
        for (int node : eliminationOrder)
        {
            cout << node << endl;
        }
    }

    cout << "Edges vector: " << endl;
    for (const Edge& edge : edges)
    {
        cout << edge.node1 << " - " << edge.node2 << " : " << edge.edge_wt << endl;
    }

    return 0;
}