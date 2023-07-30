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
    //std::cout << "Checking if set1 is a subset of set2..." << std::endl;
    //std::cout << "set1: ";
    /*
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

    */

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


/**
 *  ADDING CODE for TREE DECOMPOSITION BELOW
 * 
*/

/**
 * - Description: Creates a tree using the three properties mentioned.
   - Parameters:
     - `chordalGraph`: The chordal graph.
     - `eliminationOrder`: The elimination order of the graph.
   - Return Type: `Tree`
   - Returns: The generated tree.
   - Steps:
     1. Create an empty tree object.
     2. Create a vector to store the mapping between vertices in `chordalGraph` and `tree`.
     3. Process the vertices in the given elimination order.
        - For each vertex `v` in the elimination order:
          - Create a new node in the tree.
          - Add the new node to the mapping vector.
          - For each neighbor `w` of `v` that has not been processed yet:
            - Create a new node in the tree.
            - Add an edge between `v` and `w` in the tree.
            - Add the new node to the mapping vector.
     4. Return the generated tree.
*/


struct Tree {
    int V;

    struct TreeNode {
        int data;
        vector<int> bag_nodes;
        vector<TreeNode*> children;
        bool flg = 0;

        TreeNode(int k=0) {
            data = k;
            //bag_nodes.push_back(k);
        }

        bool isVertexInTreeNode(int x){
            // returns true if x in current TreeNode else false
            for(int i = 0; i < bag_nodes.size(); i++){
                if(bag_nodes[i] == x)
                return true;
            }
            return false;
        }

        bool isNodeChild(TreeNode *p){
            // return true if p is already a child of current TreeNode else false
            for(int i = 0; i < children.size(); i++){
                if(children[i]->data == p->data){
                    return true;
                }
            }
            return false;
        }

        void addVertInTreeNode(int x){
            // adds the vertex x to the bag of tree nodes
            bag_nodes.push_back(x);
        }

        void addChild(TreeNode *p){
            if(isNodeChild(p)){
                return;
            }
            //cout << "ADDED child " << endl;
            children.push_back(p);
        }

        int getNumCommonVertices(TreeNode *p){
            // finds and returns the number of common vertices between this node and 
            // treeNode p
            int numIntersection = 0;
            for(int i = 0; i < bag_nodes.size(); i++){
                if(find(p->bag_nodes.begin(), p->bag_nodes.end(), bag_nodes[i]) != p->bag_nodes.end()){
                    numIntersection += 1;
                }
            }
            return numIntersection;
        }
        
    };

    TreeNode* root;

    std::vector<TreeNode*> nodes;
    
    Tree(int vv = 0) {
        V = vv;
        nodes.reserve(vv);
        root = nodes[0];
        for (int i = 0; i < V; ++i) {
            //nodes[i] = new TreeNode;
            nodes.push_back(new TreeNode);
            nodes[i]->data = i;
            nodes[i]->flg = false;
        }

    }


    void addVertexInTreeNode(int treeNodeIndex, int x){
        nodes[treeNodeIndex]->addVertInTreeNode(x);
    }

    int addNode(TreeNode p) {
        for (int i = 0; i < V; i++) if (nodes[i]->data == p.data) {
            nodes[i]->flg = true;
            return i;
        }
        return -1;
    }

    /*void addEdge(int i, int k) {
        nodes[i]->children.push_back(nodes[k]);
    }*/

    int getNodeIndex(int w) {
        for (int i = 0; i < V; i++) if (nodes[i]->data == w) return i;
    }

    void print(){
        for(int i = 0; i < V; i++){
            if(nodes[i]->flg){
                cout << "Tree Node  " << nodes[i]->data << " - ";
                for(int j = 0; j < nodes[i]->bag_nodes.size(); j++) {
                    cout << nodes[i]->bag_nodes[j] << " , " ;
                } 
                cout << endl;
            }
        }
    }

    void addEdge(int parentIdx, int childIdx){
        // if nodes[childIdex] is not already present as a child of nodes[parentIdx] then add it in children
        //.. list of nodes[parentIdx]
        nodes[parentIdx]->addChild(nodes[childIdx]);
    }


    void addEdgesBetweenTreeNodes(){
        for(int i = 0; i < V; i++){
            if(nodes[i]->flg){
                int maxIntersection = 0;
                int nodeWithMaxIntersection = -1;
                for(int j = i+1; j < V; j++){
                    if(nodes[j]->flg){
                        // check the number of vertices common between treeNode i and j
                        if(maxIntersection < nodes[i]->getNumCommonVertices(nodes[j])){
                            maxIntersection = nodes[i]->getNumCommonVertices(nodes[j]);
                            nodeWithMaxIntersection = j;
                        }
                    }
                }

                if(maxIntersection != 0){
                    // add edge with TreeNode i and TreeNode[nodeWithMaxIntersection]
                    addEdge(i, nodeWithMaxIntersection);
                }
            }
        }
    }



    /*void checkAndAddEdge(int x, int y){
        for(int i = 0; i < V; i++){
            for(int j = i+1 ; j < V; j++){
                if(nodes[i]->flg && nodes[j]->flg){
                    // check if edge(x,y) is present among two treeNodes i and j
                    if(nodes[i]->isVertexInTreeNode(x) && nodes[i]->isVertexInTreeNode(y)){
                        if(nodes[j]->isVertexInTreeNode(x) && nodes[j]->isVertexInTreeNode(y)){
                            // include nodes[j] as a children of nodes[i]
                            //cout << "ADDING EDGE" << endl;
                            addEdge(i, j);
                        }
                    }

                }
            }
        }
    }*/
};


struct Graph {
    int V;
    vector<vector<int> > adj;

    Graph(int v,vector<Edge> Edges) {
        //cout <<" entered grapph" << endl;
        V = v;
        //adj.reserve(V+3);
        
        for(int i = 0; i < V+3; i++){
            vector<int> temp_vec;
            adj.push_back(temp_vec);
        }

        for (auto it : Edges) {
            //cout << " edge from " << it.node1 << " TO " << it.node2 << endl; 
            //cout << "Printing vect contents " << adj[it.node2].size() << endl;
            adj[it.node1].push_back(it.node2);
            adj[it.node2].push_back(it.node1);
            
        }
        //cout << "Done with graph " << endl;
    }
};


//wikepedia: https://en.wikipedia.org/wiki/Tree_decomposition
//Create a structure for Tree
Tree *generateTree(const Graph& chordalGraph, const vector<int>& eliminationOrder) {
    int V = chordalGraph.V;
    cout << "number of vertices in chordal graph " << V << endl;
    cout << "number of vertices in elimination order " << eliminationOrder.size() << endl;
    Tree *tree = new Tree(V);

    // Create a map to store the mapping between vertices in chordalGraph and tree
    vector<int> mapping(V, -1);

    // Process vertices in the given elimination order
    for (int i = 0; i < eliminationOrder.size(); ++i) {
        int v = eliminationOrder[i];
        Tree::TreeNode newNode(v);
        int treeNodeIndex = tree->addNode(newNode);
        if(treeNodeIndex != -1){
             mapping[v] = treeNodeIndex;
             tree->addVertexInTreeNode(mapping[v], v);
        }
        cout << " mapping of " << v << " is - " << mapping[v] << endl;
        // Process neighbors of v that haven't been processed yet
        for (int w : chordalGraph.adj[v]) {
            if (mapping[w] == -1) {
                // add w to the bag of vertices in TreeNode that has vertex v
                tree->addVertexInTreeNode(mapping[v], w);
                //Tree::TreeNode newNeighbor(w);
                //tree.addNode(newNeighbor);
                //tree.addEdge(i, tree.getNodeIndex(w));
                //mapping[w] = tree.getNodeIndex(w);
            }
        }
    }

    // print bag of verttices (or treeNodes)
    tree->print();

    // add edges now  ---- UPDATE THIS AS PER NEW ALGORITHM
    /** AS per the algorithm, a tree node x will have edges with only those tree nodes with which 
       it shares the maximum number of vertices. 
    **/
    tree->addEdgesBetweenTreeNodes();
    /*for(int x = 0; x < V; x++){
        for(int y: chordalGraph.adj[x]){
            // check if edge (x,y) is commong among two tree Nodes
            tree->checkAndAddEdge(x,y);
        }
    }*/

    cout << "Size of tree nodes " << tree->nodes.size() << endl;

    return tree;
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
        vertex_N = max(vertex_N, max(node1, node2));
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

    //cout << "Generating chordal graph" << endl;
    Graph chordalGraph(vertex_N, chordalEdges);
    //cout << "Fine till here" << endl;

    Tree *result_tree = generateTree(chordalGraph, eliminationOrder);
    //cout << "Chordal Edges: " << endl;


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

    //cout << result_tree->nodes.size();
    cout << "Tree output" << endl;
    for (auto node : result_tree->nodes) {
        cout << "node " << node->data << " has " << node->bag_nodes.size()<<" vertexes:" << endl;;
        for (int i = 0; i < node->bag_nodes.size(); i++) cout << node->bag_nodes[i] << " " ;;
        cout << endl;
    }

    cout << "Printing edges between Tree Nodes below: " << endl;
    for (auto node : result_tree->nodes) {
        if (node->children.size() > 0){
            cout << "Tree Node " << node->data << " had edges to following tree nodes: " << endl;
            for(int i = 0; i < node->children.size(); i++){
                cout << node->children[i]->data << " "; 
            }
            cout << endl;
        }
        
    }


    return 0;
}