
// INPUT HEADERS
#include "translate_from_input.hpp"
#include "input_to_network.hpp"
#include "structure_defs.hpp"
#include "graph_core.hpp"
#include "basic_array.hpp"

#include "struct_def.hpp"



// OUTPUT HEADERS
#include "printout_network.hpp"
#include "printout_others.hpp"


using namespace std;

int main(int argc, char *argv[])
{
    clock_t q, q1, q2, t;
    vector<Edge> temp_edges;

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

    print_network(x1);
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
        vertex_N = max(vertex_N, max(node1, node2));
        Edge edge;
        edge.node1 = node1;
        edge.node2 = node2;
        edge.edge_wt = weight;
        if (edge.node1 < edge.node2)
        {
            temp_edges.push_back(edge);
        }
    }
    dataFile.close();

    // print(edges);
    Array<Edge> edges(temp_edges.size());
    for(int i = 0; i < temp_edges.size(); i++){
        edges.push_back(temp_edges[i]);
    }

    vertex_N = vertex_N + 1;

std::pair<Array<Edge>, Array<int>> result = findChordalEdgesWithEliminationOrder(edges, nodes);

    Array<Edge> chordalEdges = result.first;
    Array<int> eliminationOrder = result.second;

    // cout << "Generating chordal graph" << endl;
    Graph chordalGraph(vertex_N, chordalEdges);
    // cout << "Fine till here" << endl;

    Tree *result_tree = generateTree(chordalGraph, eliminationOrder);
    // cout << "Chordal Edges: " << endl;

    cout << "Chordal Edges: " << endl;
    if (chordalEdges.empty())
    {
        cout << "No chordal edges found." << endl;
    }
    else
    {
        for(int i = 0; i < chordalEdges.size(); i++){
            auto edge = chordalEdges[i];
            cout << edge.node1 << " - " << edge.node2 << endl;
        }
        /*for (const Edge &edge : chordalEdges)
        {
            cout << edge.node1 << " - " << edge.node2 << endl;
        }*/
    }

    cout << "Elimination Order: " << endl;
    if (eliminationOrder.empty())
    {
        cout << "No elimination order found." << endl;
    }
    else
    {
        for(int i = 0; i < eliminationOrder.size(); i++){
            int node = eliminationOrder[i];
            cout << node << endl;
        }
        
    }

    cout << "Edges vector: " << endl;

    for(int i = 0; i < edges.size(); i++){
        auto edge = edges[i];
        cout << edge.node1 << " - " << edge.node2 << " : " << edge.edge_wt << endl;
    }

    cout << "Tree output" << endl;
    for(int i = 0; i < result_tree->nodes.size(); i++){
        auto n = result_tree->nodes[i];
        if(n->bag_nodes.size() > 0){
            cout << "node " << n->data << " has " << n->bag_nodes.size() << " vertexes:" << endl;
            for (int j = 0; j < n->bag_nodes.size(); j++)
                cout << n->bag_nodes[j] << " ";
            cout << endl;
        }
    }

    cout << "Printing edges between Tree Nodes below: " << endl;
    for(int i = 0; i < result_tree->nodes.size(); i++){
        auto n = result_tree->nodes[i];
        if (n->neighbours.size() > 0)
        {
            cout << "Tree Node " << n->data << " had edges to following tree nodes: " << endl;
            for (int j = 0; j < n->neighbours.size(); j++)
            {
                cout << n->neighbours[j]->data << " ";
            }
            cout << endl;
        }
    }


    return 0;
}