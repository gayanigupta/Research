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
int eliminationOrder[6]= {0,3,2,4,5,1};

    // Create a map to store the mapping between vertices in chordalGraph and tree nodes
    Array<int> mapping(V, -1);

    // Create a vector to keep track of visited vertices
    Array<bool> visited(V, false);
    //vector<bool> visited(V, false);

    // Process vertices in the given elimination order
    for(int i = 0; i < eliminationOrder.size(); ++i){
        int v = eliminationOrder[i];
        cout <<"my V:" << v <<"\n";
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

