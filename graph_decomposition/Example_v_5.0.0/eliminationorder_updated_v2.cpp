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

/**
 *  START OF FUNCTION
 *
 */
void print(vector<Edge> edges)
{
    for (Edge e : edges)
    {
        cout << e.node1 << " " << e.node2 << " " << e.edge_wt << endl;
    }
}

/**
 *  Below function returns true if all entries of 'Q1Array' are false else returns false
 */
bool isBoolVectorEmpty(const std::vector<bool> &Q1Array)
{
    for (bool b : Q1Array)
    {
        if (b)
        {
            return false;
        }
    }
    return true;
}

/**
 *  END OF FUNCTION
 *
 */

/**
 * - Description: Checks if one set is a subset of another set.
   - Parameters: `set1` and `set2` are the sets to be compared.
   - Return Type: `bool`
   - Returns: `true` if `set1` is a subset of `set2`, `false` otherwise.

*/
bool isSubset(const std::set<int> &set1, const std::set<int> &set2)
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

void getNeighbors(int vertex, const std::vector<Edge>& edges, std::vector<int>& neighbors) {
    neighbors.clear(); // Clear the vector before populating it with new neighbors

    for (const Edge& edge : edges) {
        if (edge.node1 == vertex) {
            neighbors.push_back(edge.node2);
        } else if (edge.node2 == vertex) {
            neighbors.push_back(edge.node1);
        }
    }
}

std::pair<std::vector<Edge>, std::vector<int> > findChordalEdgesWithEliminationOrder(const std::vector<Edge>& edges, int numVertices) {
    std::vector<Edge> chordalEdges;
    std::vector<int> eliminationOrder;

    std::vector<int> LP(numVertices, 0); // Lowest parent initialization set to 0
    std::vector<std::set<int> > C(numVertices, std::set<int>());

    std::vector<bool> Q1Array(numVertices, false);
    std::vector<bool> newQ1Array(numVertices, false);

    // Algorithm 1 - Loop over all vertices to identify the neighbors of each vertice
    vector<int> neighbors;
    for (int v = 0; v < numVertices; ++v) {
        // Step 1: Get the neighbors of vertex v
        getNeighbors(v,edges,neighbors);
        int lowest = std::numeric_limits<int>::max();
        int lowestParent =0;
        cout << "Neigbors of selected vertex "<< v << " : "<<endl;
        for (int neighbor : neighbors) {
        
                // Step 2: Find the neighbor with the smallest id and is smaller than v
                 if(neighbor< lowest && neighbor < v ) {
                    lowest=neighbor;
                    //cout<< "lowest: " << lowest <<endl;
                    lowestParent = lowest;
                   // cout<< "lowest parent: " << lowestParent <<endl;
                    //cout << "Lowest parent of vertex " << v << ": " << lowest << endl;
                 } else {
                // Step 3: If a lowest parent is found, print it; otherwise, set lowest parent to 0
                    //cout << "No lowest parent exists for vertex " << v << ", setting lowest parent to 0" << endl;
                    lowestParent = 0;

                 }
                 
                cout << neighbor << "  ";

              
          }
 std::cout << endl;

        // Step 4: Set LP of v to its lowest parent
       
        LP[v] = lowest==std::numeric_limits<int>::max()?0:lowest; // Set LP of v to its next lowest parent

        Q1Array[v] = true;
    
}

        // Print the lowestParent and LP at each iteration
      int i =1;
      std::cout << "\n After calculating the Lowest parent are : "<<endl;
        for (i=1; i < numVertices; i++) {
          
              cout << "Vertices : " << i <<" Lowest Parent : " << LP[i] <<endl;
            
        }
        cout<<endl;


    // Algorithm 1 - Iterations
    int iteration = 2;
    while (!isBoolVectorEmpty(Q1Array)) {
        std::fill(newQ1Array.begin(), newQ1Array.end(), false);

        for (int v = 0; v < numVertices; ++v) {
            if (Q1Array[v]) {
                int lowestParent = LP[v]; // Initialize the lowestParent to LP[v]

                // Algorithm 1 - Check edges for chordal properties
                for (const Edge& edge : edges) {
                    int w = edge.node2;
                    int u = edge.node1;
                    double wt = edge.edge_wt;

                    // Check if v is the lowest parent of w and v is lower than w
                    if (LP[w] == v && v < w && isSubset(C[w], C[v])) {
                        // Update lowest parent if id is lower than the current lowestParent
                        if (w < lowestParent) {
                            lowestParent = w;
                        }

                        C[v].insert(w);
                        chordalEdges.push_back(edge);
                        newQ1Array[w] = true;
                    }
                }

                if (lowestParent != LP[v]) {
                    LP[v] = lowestParent; // Algorithm 1 - Set LP of v to its next lowest parent
                }

                // Print the lowestParent and LP at each iteration
                std::cout << "Iteration " << iteration << ": Lowest Parent: " << lowestParent << ", LP: ";
                for (int i = 0; i < numVertices; ++i) {
                    std::cout << LP[i] << " ";
                }
                std::cout << std::endl;
            }
        }

        Q1Array = newQ1Array;
        iteration++;
    }

    // Algorithm 1 - Building the maximal chordal subgraph
    for (int v = 0; v < numVertices; ++v) {
        if (LP[v] == 0) {
            eliminationOrder.push_back(v);
        }
    }
    


    return std::make_pair(chordalEdges, eliminationOrder);
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
        if (edge.node1 < edge.node2){
        edges.push_back(edge);
        }
    }
    dataFile.close();

    // print(edges);

    std::pair<std::vector<Edge>, std::vector<int> > result = findChordalEdgesWithEliminationOrder(edges, nodes);
    std::vector<Edge> chordalEdges = result.first;
    std::vector<int> eliminationOrder = result.second;

    cout << "Chordal Edges: " << endl;
    if (chordalEdges.empty())
    {
        cout << "No chordal edges found." << endl;
    }
    else
    {
        for (const Edge &edge : chordalEdges)
        {
            cout << edge.node1 << " - " << edge.node2 << endl;
        }
    }

    cout << "Elimination Order: " << endl;
    if (eliminationOrder.empty())
    {
        cout << "No elimination order found." << endl;
    }
    else
    {
        for (int node : eliminationOrder)
        {
            cout << node << endl;
        }
    }

    cout << "Edges vector: " << endl;
    for (const Edge &edge : edges)
    {
        cout << edge.node1 << " - " << edge.node2 << " : " << edge.edge_wt << endl;
    }

    return 0;
}
