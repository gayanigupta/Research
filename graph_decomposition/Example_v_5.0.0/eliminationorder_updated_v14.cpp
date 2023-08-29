#include <iostream>
#include <list>
#include <set>
#include <map>
#include <array>

// INPUT HEADERS
#include "translate_from_input.hpp"
#include "input_to_network.hpp"
#include "structure_defs.hpp"

// OUTPUT HEADERS
#include "printout_network.hpp"
#include "printout_others.hpp"

using namespace std;


bool operator==(const Edge &e1, const Edge &e2)
{
    return ((e1.node1 == e2.node1 && e1.node2 == e2.node2) || (e1.node1 == e2.node2 && e1.node2 == e2.node1)) && e1.edge_wt == e2.edge_wt;
}

/**
 * Making a generic Array Class which will be substitute vector in the program
 * 
*/
template<typename T>
class Array {
private:
    T* arr;
    int max_size = 0;
    int cur_size = 0;

public:
        Array(){
            cur_size = 0;
            max_size = 0;
        }

        Array(int n){
            arr = new T[n];
            max_size = n;
            cur_size = 0;
        }

        Array(T ar[], int n){
            arr = new T[n];
            max_size = n;
            for(int i = 0; i < n; i++){
                arr[i] = ar[i];
            }
            cur_size = n;
        }

        Array(int n, T val){
            arr = new T[n];
            max_size = n;
            for(int i = 0; i < n; i++){
                arr[i] = val;
            }
            cur_size = n;
        }

// overriding default copy constructor for deep copy
Array(const Array<T> &obj) {
    cur_size = obj.cur_size;
    max_size = obj.max_size;
    arr = new T[cur_size];
    for (int i = 0; i < cur_size; i++) {
        arr[i] = obj.arr[i];
    }
}

        // overloading []
        T& operator[](int index){
            if(index >= max_size){
                cout << "Index out of bound" << endl;
                exit(0);
            }
            return arr[index];
        }

        // this functions allots a memory of size s to an array which has no space
        void makeArray(int s){
            if(max_size == 0){
                arr = new T[s];
                max_size = s;
            }
        }

        void push_back(T val){
            if(cur_size >= max_size){
                cout << "No more space to add an extra element" << endl;
                exit(0);
            }
            arr[cur_size] = val;
            cur_size += 1;
        }

        // remove 'val' from the array if present
        void erase(T val){
            for(int i = 0; i < cur_size; i++){
                if(arr[i] == val){
                    // replace the current value with last element 
                    arr[i] = arr[cur_size-1];
                    //arr[cur_size - 1] = -1;
                    cur_size = cur_size - 1;
                    break;
                }
            }
        }

        // remove 'val' from array if present and shift the subsequent elements by 1 position to the left
        void eraseWithOrder(T val){
            for(int i = 0; i < cur_size; i++){
                if(arr[i] == val){

                    for(int j = i; j < cur_size-1; j++){
                        arr[j] = arr[j+1];
                    }
                    //arr[cur_size - 1] = -1;
                    cur_size = cur_size - 1;
                    break;
                }
            }
        }

        // reverse the entries of the array 
        void reverse(){
            for(int i = 0; i <= cur_size/2 - 1; i++){
                T temp; 
                temp = arr[i];
                arr[i] = arr[cur_size - 1 - i];
                arr[cur_size - 1- i] = temp;
            }
        }

        // clear the elements of the array
        void clear(){
            cur_size = 0;
        }

        int size(){
            return cur_size;
        }

        // returns the index of the 'val' if present, else returns -1
        int getIndex(T val){
            for(int i = 0; i < cur_size; i++){
                if(arr[i] == val){
                    return i;
                }
            }
            return -1;
        }

        // returns the value at stored at the 'idx' element
        T at(int idx){
            if(idx >= cur_size){
                cout << "Index out of bounds" << endl;
                exit(0);
            }
            return arr[idx];
        }

        // returns true if the array is empty
        bool empty(){
            if(cur_size == 0){
                return true;
            }
            return false;
        }
};


bool isSubset(const std::set<int> &set1, const std::set<int> &set2)
{

    for (int element : set1)
    {
        if (set2.find(element) == set2.end())
        {
            return false;
        }
    }

    return true;
}

void getNeighbours(int vertex, Array<Edge> &edges, map<int, Array<int>> &adjList, int numVertices)
{
    if (adjList.find(vertex) == adjList.end())
    {
        Array<int> temp(numVertices);
        adjList[vertex] = temp;
    }

    for(int i = 0; i < edges.size(); i++){
        Edge edge = edges[i];
        if (edge.node1 == vertex)
        {
            adjList[vertex].push_back(edge.node2);
        }
        else if (edge.node2 == vertex)
        {
            adjList[vertex].push_back(edge.node1);
        }

    }

}

int getLeftMostElem(Array<int> **setPartition, int &setPartSize)
{
    int elem = setPartition[0]->at(0);
    // remove the 'elem' from the set we got
    setPartition[0]->eraseWithOrder(elem);
    //setPartition[0].erase(setPartition[0].begin());

    // if the first set  has gone empty, remove it too
    if(setPartition[0]->size() == 0){
        for(int i = 0; i < setPartSize - 1; i++){
            setPartition[i]->clear();
            for(int j = 0; j < setPartition[i+1]->size(); j++){
                setPartition[i]->push_back(setPartition[i+1]->at(j));
            }
            //setPartition[i] = setPartition[i+1];
        }
        setPartSize -= 1;
    }

    return elem;
}


void updateSet(Array<int> neighbours, Array<int> **setPartition, int &setPartSize, int V)
{
    Array<int> **tempSets = new Array<int>*[V];
    for(int i = 0; i < V; i++){
        tempSets[i] = new Array<int>(V);
    }
    int tempSetSize = 0;
    
    for (int i = 0; i < setPartSize; i++)
    {
        Array<int> left(V);
        Array<int> right(V);

        for (int j = 0; j < setPartition[i]->size(); j++)
        {
            if( neighbours.getIndex(setPartition[i]->at(j)) != -1){
                // put that elem in the left
                left.push_back(setPartition[i]->at(j));
            }else{
                right.push_back(setPartition[i]->at(j));
            }
     
        }

        if (left.size() != 0){
            for(int k = 0; k < left.size(); k++){
                //cout << "  left #"  << k << " = " << left[k] << endl;
                tempSets[tempSetSize]->push_back(left[k]);
            }
            tempSetSize += 1;
        }

        //cout << " tempSetSize = " << tempSetSize << endl;

        if (right.size() != 0){
            for(int k = 0; k < right.size(); k++){
                //cout << "  right #"  << k << " = " << right[k] << endl;
                tempSets[tempSetSize]->push_back(right[k]);
            }
            tempSetSize += 1;
        }


    }


    // copy the tempSets back to setPartition
    //setPartition.clear();
    setPartSize = tempSetSize;
    for(int i = 0; i < tempSetSize; i++){
        setPartition[i]->clear();
        for(int j = 0; j < tempSets[i]->size(); j++){
            //cout << "     tempSets[" << i << "][" << j << "] = " << tempSets[i]->at(j) ;
            setPartition[i]->push_back(tempSets[i]->at(j));

        }
    }

}

std::pair< Array<Edge>, Array<int>> findChordalEdgesWithEliminationOrder(Array<Edge> &edges, int numVertices)
{
    Array<Edge> chordalEdges(edges.size());
    Array<int> eliminationOrder(numVertices);
    map<int, int> visited;

    map<int, Array<int>> adjList;

    Array<int> ** setPartition = new Array<int>*[numVertices];
    for(int i = 0; i < numVertices; i++){
        setPartition[i] = new Array<int>(numVertices);
    }
    //vector<vector<int>> setPartition;

    Array<int> initSet(numVertices);
    //vector<int> initSet;

    int setPartSize = 0;
    for (int v = 0; v < numVertices; v++)
    {
        // get the neighbours of each vertex and store them in an adjacency list
        getNeighbours(v, edges, adjList, numVertices);
        // initially set visited of all vertices to False(0)
        visited[v] = 0;
        initSet.push_back(v);
    }


    // get chordal edges
    for (int i = 0; i < edges.size(); i++){
        Edge edge = edges[i];
        chordalEdges.push_back(edge);
    }
    /*for (const Edge &edge : edges)
    {
        chordalEdges.push_back(edge);
    }*/

    for(int i = 0; i < initSet.size(); i++){
        setPartition[0]->push_back(initSet[i]);
    }
    setPartSize += 1;
    //setPartition.push_back(initSet);


    while(setPartSize != 0){
        // get the leftmost element from the leftmost set
        int u = getLeftMostElem(setPartition, setPartSize);
        // cout << "Current element in Lex BFS = " << u << endl;
        // mark u as visited
        visited[u] = 1;
        // put u in the partially created elimination order
        eliminationOrder.push_back(u);

        // check u's neighbours and update the sets
        Array<int> neighbours(numVertices);
        for(int i = 0 ; i < adjList[u].size(); i++){
            int neighbour = adjList[u][i];
            if (visited[neighbour] == 0)
            {
                neighbours.push_back(neighbour);
            }
        }

        updateSet(neighbours, setPartition, setPartSize, numVertices);
    }

    // reverse the elements in the 'eliminationOrder' to get PEO
    eliminationOrder.reverse();
    //reverse(eliminationOrder.begin(), eliminationOrder.end());

    // print the elements of the elimination order
    cout << "Elimination order below \n";
    for (int i = 0; i < eliminationOrder.size(); i++)
    {
        cout << eliminationOrder[i] << endl;
    }

    return std::pair<Array<Edge>, Array<int>>(chordalEdges, eliminationOrder);

}

struct node
{
    int data;
    Array<int> bag_nodes;
    Array<node *> neighbours;
    bool flg = 0;

    node(int k , int V)
    {
        data = k;
        bag_nodes.makeArray(V);
        neighbours.makeArray(V);
        // bag_nodes.push_back(k);
    }

    node(node *nodeToCopy, int V)
    {
        data = nodeToCopy->data;
        bag_nodes.makeArray(V);
        neighbours.makeArray(V);
        // copy the internal contents
        for (int i = 0; i < nodeToCopy->bag_nodes.size(); i++)
        {
            bag_nodes.push_back(nodeToCopy->bag_nodes[i]);
        }
        flg = nodeToCopy->flg;
    }

    bool isVertexInNode(int x)
    {
        // returns true if x in current TreeNode else false
        for (int i = 0; i < bag_nodes.size(); i++)
        {
            if (bag_nodes[i] == x)
                return true;
        }
        return false;
    }

    bool isNodeNeighbour(node *p)
    {
        // return true if p is already a child of current TreeNode else false
        for (int i = 0; i < neighbours.size(); i++)
        {
            if (neighbours[i]->data == p->data)
            {
                return true;
            }
        }
        return false;
    }

    void addVertInNode(int x)
    {
        // adds the vertex x to the bag of tree nodes
        bag_nodes.push_back(x);
    }

    void addNeighbour(node *p)
    {
        if (isNodeNeighbour(p))
        {
            return;
        }
        // cout << "ADDED child " << endl;
        neighbours.push_back(p);
    }

    int getNumCommonVertices(node *p)
    {
        // finds and returns the number of common vertices between this node and
        // treeNode p
        int numIntersection = 0;
        for (int i = 0; i < bag_nodes.size(); i++)
        {
            if( p->bag_nodes.getIndex(bag_nodes[i]) != -1){
                numIntersection += 1;
            }
            /*
            if (find(p->bag_nodes.begin(), p->bag_nodes.end(), bag_nodes[i]) != p->bag_nodes.end())
            {
                numIntersection += 1;
            }*/
        }
        return numIntersection;
    }
};

struct Tree
{
    int V;

    node *root;

    Array<node *> nodes;

    Tree(int vv = 0)
    {
        V = vv;
        nodes.makeArray(vv);
        //nodes.reserve(vv);
        root = nodes[0];
        for (int i = 0; i < V; ++i)
        {
            // nodes[i] = new TreeNode;
            nodes.push_back(new node(i, V));
            nodes[i]->data = i;
            nodes[i]->flg = false;
        }
    }

    Tree(Tree *treeToCopy)
    {
        V = treeToCopy->V;
        nodes.makeArray(V);
        //nodes.reserve(V);
        root = nodes[0];
        for (int i = 0; i < V; i++)
        {
            nodes.push_back(new node(treeToCopy->nodes[i], V));
        }
    }

    void addVertexInTreeNode(int treeNodeIndex, int x)
    {
        nodes[treeNodeIndex]->addVertInNode(x);
    }

    int addNode(node p)
    {
        for (int i = 0; i < V; i++)
            if (nodes[i]->data == p.data)
            {
                nodes[i]->flg = true;
                return i;
            }
        return -1;
    }

    int getNodeIndex(int w)
    {
        for (int i = 0; i < V; i++)
            if (nodes[i]->data == w)
                return i;
    }

    void print()
    {
        for (int i = 0; i < V; i++)
        {
            if (nodes[i]->flg)
            {
                cout << "Tree Node  " << nodes[i]->data << " - ";
                for (int j = 0; j < nodes[i]->bag_nodes.size(); j++)
                {
                    cout << nodes[i]->bag_nodes[j] << " , ";
                }
                cout << endl;
            }
        }
    }

    void addEdge(int parentIdx, int childIdx)
    {
        // add edge
        nodes[parentIdx]->addNeighbour(nodes[childIdx]);

        /** Add a back neighbour too for it too be undirected graph  */
        nodes[childIdx]->addNeighbour(nodes[parentIdx]);
    }

    void addEdge(node *x, node *y)
    {
        int p = -1, q = -1;
        for (int i = 0; i < V; i++)
        {
            if (nodes[i]->flg && nodes[i]->data == x->data)
            {
                p = i;
            }

            if (nodes[i]->flg && nodes[i]->data == y->data)
            {
                q = i;
            }
        }

        if (p != -1 && q != -1)
        {
            addEdge(p, q);
        }
    }

    void addEdgesBetweenNodes()
    {
        // We add edge between two tree nodes if they share a common vertex v in the original chordal graph.
        // The resultant structure may be a cyclic graph instead of tree, however, that's not an issue since
        // we will build a minimum spanning tree from clique graph.
        for (int i = 0; i < V; i++)
        {
            if (nodes[i]->flg)
            {
                for (int j = i + 1; j < V; j++)
                {
                    if (nodes[j]->flg)
                    {
                        // if node i and j have the common vertex then add an edge between them
                        if (nodes[i]->getNumCommonVertices(nodes[j]) > 0)
                        {
                            addEdge(i, j);
                        }
                    }
                }
            }
        }
    }
        // Method to clear a node and its associated vertices
void clearNode(int nodeIndex)
{
    if (nodeIndex >= 0 && nodeIndex < nodes.size())
    {
        nodes[nodeIndex]->bag_nodes.clear(); // Clear the bag of vertices for this node
        // Set flg to false for all nodes within this tree node
        for(int i = 0; i < nodes[nodeIndex]->neighbours.size(); i++){
            node *n = nodes[nodeIndex]->neighbours[i];
            n->flg = false;
        }
    }
}
};

struct Graph
{
    int V;
    Array<int> **adj;
    //vector<vector<int>> adj;

    Graph(int v, Array<Edge> Edges)
    {
        // cout <<" entered grapph" << endl;
        V = v;
        // adj.reserve(V+3);
        adj = new Array<int>*[V+3];

        for(int i = 0; i < V+3; i++){
            adj[i] = new Array<int>(V+3);
        }

        /*
        for (int i = 0; i < V + 3; i++)
        {
            vector<int> temp_vec;
            adj.push_back(temp_vec);
        }*/
        for(int i = 0; i < Edges.size(); i++){
            auto it = Edges[i];
            adj[it.node1]->push_back(it.node2);
            adj[it.node2]->push_back(it.node1);
        }

    }

    bool isEdge(int x, int y)
    {
        if(adj[x]->getIndex(y) != -1){
            return true;
        }
        return false;
    }
};


// Below function retuns the maximal clique from set uncheckedNeighbours
Array<int> getMaximalClique(Graph &chordalGraph, Array<int> uncheckedNeighbours)
{

    if (uncheckedNeighbours.size() <= 1)
    {
        return uncheckedNeighbours;
    }
    else
    {
        // check if all the vertices in the set form a clique
        bool isClique = true;
        int n = uncheckedNeighbours.size();
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (!chordalGraph.isEdge(uncheckedNeighbours[i], uncheckedNeighbours[j]))
                {
                    isClique = false;
                    break;
                }
            }
            if (!isClique)
            {
                break;
            }
        }

        if (isClique)
        {
            return uncheckedNeighbours;
        }
        else
        {

            // recursively find the maximal clique from uncheckedNeighbours
            Array<int> maxClique(uncheckedNeighbours.size());
            int maxCliqueSize = 1;
            Array<int> neighbourSubset = uncheckedNeighbours;
            // try removing one element at a time from neighbourSubset and recursively check for maximal clique
            for (int i = 0; i < n; i++)
            {
                int it = neighbourSubset.getIndex(uncheckedNeighbours[i]);
                neighbourSubset.erase( neighbourSubset[it] );
                // recursively find the maximal clique from the smaller vector
                Array<int> result = getMaximalClique(chordalGraph, neighbourSubset);

                if (maxCliqueSize < result.size())
                {
                    maxCliqueSize = result.size();
                    maxClique = result;
                }
                // push back the element that we just deleted
                neighbourSubset.push_back(uncheckedNeighbours[i]);
            }

            return maxClique;
        }
    }
}

// returns the index where the node pointed by n exists in 'maxEdge' else returns -1
int findIfNodeExists(node *n, Array<pair<node *, int>> maxEdge)
{
    for (int i = 0; i < maxEdge.size(); i++)
    {
        if (n->data == maxEdge[i].first->data)
        {
            return i;
        }
    }
    return -1;
}
// The below function returns the maximum spanning tree of the clique graph
// Prim's alogorithm is used to find MST
Tree *getMaxSpanTree(Tree *cliqueGraph){
    Tree *mst = new Tree(cliqueGraph);
    map<node*, node*> parent;
    Array< pair<node*, int> > maxEdge(cliqueGraph->V);
    //vector< pair<node*, int> > maxEdge;

    bool isFirst = true;
    for(int i = 0; i < cliqueGraph->nodes.size(); i++){
        if(cliqueGraph->nodes[i]->flg){
            
            parent[cliqueGraph->nodes[i]] = NULL;
            //not.push_back(cliqueGraph->nodes[i]);

            pair<node*, int> p;
            p.first = cliqueGraph->nodes[i];
            if(isFirst){
                p.second = 0;
                isFirst = false;
            }else{
                p.second = -1;
            }
            maxEdge.push_back(p);
        }
    }

    // run the prim's algorithm
    //vector< pair<node*, int> >::iterator it, maxIt;
    int maxIt;
    int maxEdgeWeight ;
    while(maxEdge.size() != 0){

        maxEdgeWeight = -1;
        maxIt = 0;

        //get the node with max edge weight
        for(int i = 0; i < maxEdge.size(); i++){
            if(maxEdgeWeight < maxEdge[i].second){
                maxEdgeWeight = maxEdge[i].second;
                maxIt = i;
            }
        }


        if(parent.find(maxEdge[maxIt].first) != parent.end() && parent[maxEdge[maxIt].first] != NULL){
            // add an edge between curNode and it's parent
            mst->addEdge(maxEdge[maxIt].first, parent[maxEdge[maxIt].first]);
        }

        // update the edge weight to the neighbours if needed
        for(int i = 0; i < maxEdge[maxIt].first->neighbours.size(); i++){
            //edgeweight between tree nodes = number of common vertices in their bags
            int edgeWeight = maxEdge[maxIt].first->getNumCommonVertices(maxEdge[maxIt].first->neighbours[i]);
            // update the edge weight in the 'maxEdge' vector if needed
            int index = findIfNodeExists(maxEdge[maxIt].first->neighbours[i] , maxEdge);
            if(index != -1){
                // if current edgeWeight is more then update the weight and parent to the curNode
                if(edgeWeight > maxEdge[index].second){
                    maxEdge[index].second = edgeWeight;
                    //update the parent
                    if(parent.find(maxEdge[maxIt].first->neighbours[i]) != parent.end()){
                        parent[maxEdge[maxIt].first->neighbours[i]] = maxEdge[maxIt].first;
                    }
                }
            } 
        }
        // remove the element pointed to by 'maxIt' iterator from the maxEdge vector
        maxEdge.erase(maxEdge[maxIt]);
        //maxEdge.erase(maxIt);

    }
    return mst;
}


// Below function returns True if some neighbour of vertex v is unvisited
bool isSomeNeighbourUnvisited(Graph &chordalGraph, Array<bool> visited, int v){
    for(int i = 0; i < chordalGraph.adj[v]->size(); i++){
        int w = chordalGraph.adj[v]->at(i);
        if(!visited[w]){
            return true;
        }
    }
    return false;

}


// The below function returns the maximum spanning tree of the clique graph
// Prim's alogorithm is used to find MST
Tree *generateTree(Graph &chordalGraph,  Array<int> &eliminationOrder)
{
    int V = chordalGraph.V;  // Get the number of vertices in the chordalGraph
    Tree *tree = new Tree(V);  // Create a new Tree object with V vertices

    // Create a map to store the mapping between vertices in chordalGraph and tree nodes
    Array<int> mapping(V, -1);

    // Create a vector to keep track of visited vertices
    Array<bool> visited(V, false);
    //vector<bool> visited(V, false);

    // Process vertices in the given elimination order
    for(int i = 0; i < eliminationOrder.size(); ++i){
        int v = eliminationOrder[i];
        // make a new Tree Node for this vertex if and only some of its neighbour is still unvisited
        if(isSomeNeighbourUnvisited(chordalGraph, visited, v)){

            node newNode(v, V);
            int treeNodeIndex = tree->addNode(newNode);
            if(treeNodeIndex != -1){
                mapping[v] = treeNodeIndex;
                tree->addVertexInTreeNode(mapping[v], v);
                visited[v] = true;
            }

            cout << " mapping of " << v << " is - " << mapping[v] << endl;

            // we add those neighbours of v which form a clique 
                //vector<int> uncheckedNeighbours;
            Array<int> uncheckedNeighbours(V);
            for(int j = 0; j < chordalGraph.adj[v]->size(); j++){
                int w = chordalGraph.adj[v]->at(j);
                if(mapping[w] == -1){
                    uncheckedNeighbours.push_back(w);
                }
            }

            // Now check how many of the vertices in 'uncheckedNeighbours' form a clique
            Array<int> maxClique = getMaximalClique(chordalGraph, uncheckedNeighbours);
            
            // insert all the elements of the 'maxClique into the current treeNode'
            for(int j = 0; j < maxClique.size(); j++){
                int w = maxClique[j];
                visited[w] = true;
                tree->addVertexInTreeNode(mapping[v],w);
            }

        }

    }


    // Generate edges between the clique Graph 
    tree->addEdgesBetweenNodes();

    // Build the maximum spanning tree from the generated tree
    Tree *mst = getMaxSpanTree(tree);  // Call a function to construct the maximum spanning tree

    return mst;  // Return the maximum spanning tree
}




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