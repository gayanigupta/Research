// C++ program to find K-Cores of a graph
#include <iostream>
#include<list>
#include <algorithm>
#include <set>
//INPUT HEADERS
#include "input_to_network.hpp"
#include"structure_defs.hpp"
#include "translate_from_input.hpp"

int vertex_N;

void heapify(vector<int>& arr, int n, int i)
{
    int largest = i; // Initialize largest as root Since we are using 0 based indexing
    int l = 2 * i + 1; // left = 2*i + 1
    int r = 2 * i + 2; // right = 2*i + 2

    // If left child is larger than root
    if (l < n && arr[l] > arr[largest])
        largest = l;

    // If right child is larger than largest so far
    if (r < n && arr[r] > arr[largest])
        largest = r;

    // If largest is not root
    if (largest != i) {
        swap(arr[i], arr[largest]);

        // Recursively heapify the affected sub-tree
        heapify(arr, n, largest);
    }
}

// main function to do heap sort
void heapSort(vector<int>& arr, int n)
{
    // Build heap (rearrange array)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    // One by one extract an element from heap
    for (int i = n - 1; i >= 0; i--) {
        // Move current root to end
        swap(arr[0], arr[i]);

        // call max heapify on the reduced heap
        heapify(arr, i, 0);
    }
}

void print_network(A_Network& x, const char* title)
{
    cout << title << endl;
    for (int i = 0; i < x.size(); ++i)
    {
        ADJ_Bundle& adj = x[i];

        cout << "Vertex " << i << "'s linked vertices: ";
        for (int j = 0; j < adj.ListW.size(); ++j)
        {
            cout << adj.ListW.at(j).first;
            if (j < adj.ListW.size() - 1)
                cout << ", ";
        }
        cout << endl;
    }
    cout << endl;
}

int get_max_no_of_network(A_Network& x)
{
    int max_no = -1;
    for (int i = 0; i < x.size(); ++i)
    {
        vector <int_double>& listW = x[i].ListW;
        for (int j = 0; j < listW.size(); ++j)
        {
            if (listW[j].first > max_no)
                max_no = listW[j].first;
        }
    }
    return (max_no < x.size() - 1) ? x.size() - 1 : max_no;
}

bool checkEdge(ADJ_Bundle& bundle, int vertex_no)
{
    vector<int_double>& adj = bundle.ListW;
    for (int i = 0; i < adj.size(); ++i)
    {
        if (adj[i].first == vertex_no)
            return true;
    }

    return false;
}

bool addUndirectedEdge(ADJ_Bundle& bundle, int vertex_no)
{
    bool found = false;
    vector<int_double>& adj = bundle.ListW;
    vector<int_double>::iterator it = adj.begin();
    for (; it != adj.end(); ++it)
    {
        if ((*it).first >= vertex_no)
        {
            found = true;
            break;
        }
    }

    int_double newEdge;
    newEdge.first = vertex_no;
    newEdge.second = 1;
    if (adj.size() > 0 && adj[adj.size() - 1].first < vertex_no)
        adj.push_back(newEdge);
    else if (adj.size() == 0 || (found && (*it).first != vertex_no))
        adj.insert(it, newEdge);

    return false;
}

bool removeVertex(A_Network& x, int vertex_no)
{
    for (int i = 0; i < x.size(); ++i)
    {
        vector<int_double>& listW = x[i].ListW;
        for (vector<int_double>::iterator it = listW.begin(); it != listW.end(); ++it)
        {
            if ((*it).first == vertex_no)
            {
                listW.erase(it);
                break;
            }
        }
    }

    if (vertex_no < x.size())
        x[vertex_no].ListW.clear();

    return true;
}

bool removeDirectedEdge(A_Network& x, int start_no, int end_no)
{
    if (start_no > x.size() - 1)
        return false;

    vector<int_double>& listW = x[start_no].ListW;
    vector<int_double>::iterator it = listW.begin();
    for (; it != listW.end(); ++it)
    {
        if ((*it).first == end_no)
        {
            listW.erase(it);
            return true;
        }
    }
    return false;
}

bool removeUndirectedEdge(A_Network& x, int start_no, int end_no)
{
    return removeDirectedEdge(x, start_no, end_no) && removeDirectedEdge(x, end_no, start_no);
}

void convertDirected2UndirectedNetwork(A_Network& a, A_Network& b)
{
    int count = get_max_no_of_network(a) + 1;
    int i = 0;
    for (i = 0; i < a.size(); ++i)
        b.push_back(a[i]);

    ADJ_Bundle adj;
    for (; i < count; ++i)
    {
        adj.Row = i;
        b.push_back(adj);
    }

    for (i = 0; i < a.size(); ++i)
    {
        vector<int_double>& adj = a[i].ListW;
        for (int j = 0; j < adj.size(); ++j)
            addUndirectedEdge(b[adj[j].first], i);
    }

}

void print_network_kcores(A_Network& x, int count, int k, const char* title)
{
    cout << title << endl;

    vector<vector<int> > groups;
    for (int i = 0; i < count; ++i)
    {
        vector<int> group;
        groups.push_back(group);
    }

    for (int i = 0; i < count; ++i)
    {
        if (i >= x.size())
            groups[0].push_back(i);
        else
            groups[x[i].ListW.size()].push_back(i);
    }

    cout << k << " cores " << ": ";
    for (int i = 1; i < count; ++i)
    {
        if (groups[i].size() == 0)
            continue;

        for (int j = 0; j < groups[i].size(); ++j)
            cout << groups[i][j] << ", ";
    }
    cout << endl;
}

/*
* Calculate kcore graph and vertices removed.
* src:  input graph
* x:    kcore graph as result
* removedVertices : vertices removed as they have samller degree than k.
* outputDegrees:    check if output of degrees' changing is displayed or not.
*/
void calculateGraphCore(A_Network& src, A_Network& x, int k, vector<int>& removedVertices, bool outputDegrees)
{
    // Convert directed graph into undirected graph in order to calculate easily
    convertDirected2UndirectedNetwork(src, x);

    int count = x.size();

    // Store degrees of all vertices
    if (outputDegrees)
        cout << "---------previous degrees---------" << endl;

    vector<int> degrees(count, 0);
    vector<bool> processed(count, false);

    // Get the degrees of vertices
    for (int i = 0; i < count; i++)
    {
        degrees[i] = x[i].ListW.size();
        if (degrees[i] == 0)
            processed[i] = true;

        if (outputDegrees)
            cout << "Vertex " << i << "'s degree: " << degrees[i] << endl;
    }

    // Remove the vertices smaller than k.
    while (true)
    {
        int removed_count = 0;
        for (int i = 0; i < count; i++)
        {
            if (processed[i])
                continue;

            if (degrees[i] < k)
            {   // Remove a vertex
                vector<int_double>& listW = x[i].ListW;
                for (int j = 0; j < listW.size(); ++j)
                {
                    int end_vertex = listW[j].first;
                    if (removeDirectedEdge(x, end_vertex, i))
                    {
                        degrees[listW[j].first]--;
                    }
                }

                listW.clear();
                degrees[i] = 0;

                processed[i] = true;
                removed_count++;

                // insert a removed vertex in store.
                removedVertices.push_back(i);
            }
        }

        // finish if there is no vertex to process
        if (removed_count == 0)
            break;
    }

    if (outputDegrees)
    {
        cout << "---------post degrees---------" << endl;
        for (int i = 0; i < count; i++)
            cout << "Vertex " << i << "'s degree: " << degrees[i] << endl;
    }
}

void swapNonAdjacentNodes(A_Network& x, int node1, int node2)
{
    // Check if node1 and node2 are valid vertices
    if (node1 >= x.size() || node2 >= x.size()) {
        cout << "Invalid vertices!" << endl;
        return;
    }

    // Check if node1 and node2 are non-adjacent
    if (checkEdge(x[node1], node2) || checkEdge(x[node2], node1)) {
        cout << "The nodes are already adjacent!" << endl;
        return;
    }

    // Remove all edges of node1
    removeVertex(x, node1);

    // Remove all edges of node2
    removeVertex(x, node2);

    // Add an undirected edge between node1 and node2
    addUndirectedEdge(x[node1], node2);
}

vector<int> getDegreeDistribution(const A_Network& network) {
    
    vector<int> degreeDistribution(network.size(), 0);
    for (const auto& node : network) {
        int degree = node.ListW.size();
        //cout<<degree<<endl;
        degreeDistribution[degree]++;
    }
    return degreeDistribution;
}

void printDegreeDistribution(const vector<int>& degreeDistribution) {
    cout << "Degree Distribution:\n";
    for (int i = 0; i < degreeDistribution.size(); ++i) {
        if(degreeDistribution[i] !=0){
        cout << i << "-th degree: " << degreeDistribution[i] << endl;
    }
    }
}

vector<int_int> Vedges;
vector<int_int> VedgesBefore;
vector<int_int> degree;
vector<int_int> degreeBefore;

void network_swap(A_Network& x) {

    for (int i = 0; i < x.size(); ++i)
    {
        ADJ_Bundle& adj = x[i];

        for (int j = 0; j < adj.ListW.size(); ++j)
        {
            Vedges.push_back(int_int(i, adj.ListW.at(j).first));
        }
    }
    for (int i = 0; i < x.size(); ++i)
    {
        ADJ_Bundle& adj = x[i];

        for (int j = 0; j < adj.ListW.size(); ++j)
        {
            VedgesBefore.push_back(int_int(i, adj.ListW.at(j).first));
        }
    }

    int preBefore = -1;
    int degBefore = 0;

    int degree[1000];
    vector<int> V[1000];
    set<int_int> S;
    memset(degree, 0, sizeof(degree));

    for (int i = 0; i < VedgesBefore.size(); i++)
    {
        int first = VedgesBefore[i].first;
        int second = VedgesBefore[i].second;
        degree[first]++; V[first].push_back(second);
        degree[second]++; V[second].push_back(first);
        S.insert(int_int(first, second));
        S.insert(int_int(second, first));
    }
    vector<int_int> nonadj;
    int root = 0;
    for (int i = 0; i < vertex_N; i++)
        if (degree[i]) for (int j = i + 1; j < vertex_N; j++) if (degree[j]) {
            if (!root) root = i;
            if (S.count(int_int(i, j))) continue;
            nonadj.push_back(int_int(i, j));
        }


    queue<int> Q;
    vector<int_int> rlt;
    Q.push(root);
    bool vis[1000];
    memset(vis, 0, sizeof(vis));
    vis[root] = 1;
    int tmp;
    while (!Q.empty()) {
        tmp = Q.front(); Q.pop();
        for (int i = 0; i < V[tmp].size(); i++) if (!vis[V[tmp][i]]) {
            Q.push(V[tmp][i]);
            rlt.push_back(int_int(tmp, V[tmp][i]));
            vis[V[tmp][i]] = 1;
        }
    }
    for (int i = 0; i < vertex_N; i++) degree[i] = 0;

    cout << endl;
    cout << "--Final Result--" << endl;
    for (int i = 0; i < rlt.size();) {
        int node = rlt[i].first;
        cout << "Node" << node << "'s vertexs: ";
        while (rlt[i].first == node) {
            degree[rlt[i].first]++;
            degree[rlt[i].second]++;
            cout << rlt[i].second << " ";
            i++;
            if (i == rlt.size()) break;
        }
        cout << endl;
    }

    cout << endl;
    cout << "--Degree Distribution After Swap--" << endl;

    for (int i = 0; i < vertex_N; i++)if (degree[i]) {
        cout << "Node " << i << " has " << degree[i] << " vertexs " << endl;
    }

}

// Driver program to test methods of graph class
int main(int argc, char *argv[]) {

clock_t q, q1, q2,t;
    vector<Edge> edges;

    // read edges from file.
    /*fstream fin;
    fin.open("Tests/core_2.txt", ios::in);
    if (fin.is_open() == false)
        {return 0;}

        cout << "read file went through "<<endl;

    while (!fin.eof())
    {
        int start_no, end_no, weight;
        fin >> start_no >> end_no >> weight;
        edges.push_back(create(start_no, end_no, 1));
    }

    fin.close();
    cout << "read file "<<endl; */

    /*

    // Create a network with edges.
    A_Network x1;
    create_Network(&edges, 0, &x1, -1);
    cout << "**** \n";
    */

// Preprocess Nodes to Numbers
        //Stores file in argv[3]: store map in argv[4]
        //Vertices start from 0
    q=clock();
    //Check if valid input is given
    if ( argc < 2) { cout << "INPUT ERROR:: At least 2 inputs required. First: filename \n Second: Filetypes: 1:node_node_wt 2:node_wt_node 3:node_node 4:node_node (Only option 1 is active now) \n Third. Name of new file \n Fourth. Name of Map file\n"; return 0;}
    //Check to see if file opening succeeded
    ifstream the_file ( argv[1] ); if (!the_file.is_open() ) { cout<<"INPUT ERROR:: Could not open file\n";}
    
     A_Network x1,x2;
    int nodes=-1;
    map_int_st revmap;
        int type=atoi("1");
        translate_input(argv[1],type,argv[3],argv[4]);
        
        //Remove Duplicate Edges and Self Loops; Create Undirected Graphs
        // process_to_simple_undirected();
        q=clock()-q;
        cout << "Total Time for Preprocessing"<< ((float)q)/CLOCKS_PER_SEC <<"\n";
        
        /***** Preprocessing to Graph (GUI) ***********/
        
        
        /******* Read Graph (GUI) and Create Reverse Map*****************/
        //Obtain the list of edges.
        q=clock();
        readin_network(&x1,argv[3],nodes);
          readin_network(&x2,argv[3],nodes);
        
        //Create Reversemap
        
        nodes=x1.size();
        create_map(argv[4],&revmap);
        
        q=clock()-q;
        cout << "Total Time for Reading Network"<< ((float)q)/CLOCKS_PER_SEC <<"\n";
    
        int count = get_max_no_of_network(x1) + 1;
            int start_core = 2, end_core = 4; // count - 1;

        cout << count << endl;
    // ... Your code ...

 vector<int> degreeDistributionBefore1 = getDegreeDistribution(x1);
    cout << "Before K Core Degree Distribution Before Swapping:\n";
    printDegreeDistribution(degreeDistributionBefore1);


    // Calculate the graph per kcore.
    vector<int> kcoresPerVertex(count, 0);
    for (int k = start_core; k <= end_core; k++) {
        A_Network kcoreGraph;
        vector<int> removedVertices;

        cout << "---------" << (k - 1) << "-core graph------------" << endl;
        // Calculate kcore graph and vertices removed.
        calculateGraphCore(x1, kcoreGraph, k, removedVertices, false);
        
        cout << "\tvertices with " << (k - 1) << "-core: ";

        // Sort the removed vertices by using heap sort
        heapSort(removedVertices, removedVertices.size());
        for (int i = 0; i < removedVertices.size(); ++i) {
            // Remove the vertex with lower degree than k.
            removeVertex(x1, removedVertices[i]);

            // Store the kcore value of vertex
            kcoresPerVertex[removedVertices[i]] = k - 1;

            cout << removedVertices[i] << ", ";
        }
        cout << endl;
    }

    cout << endl;
    // Output vertices' kcore.
    cout << "---------" << "vertices' kcore------" << endl;
    for (int i = 0; i < count; ++i)
        cout << "\t" << i << "-th vertex's kcore: " << kcoresPerVertex[i] << endl;

    // Calculate and print degree distribution before and after swapping
    vector<int> degreeDistributionBefore = getDegreeDistribution(x2);
    cout << "Degree Distribution Before Swapping:\n";
    printDegreeDistribution(degreeDistributionBefore);

    // Perform swapping or other operations

    network_swap(x2);

    vector<int> degreeDistributionAfter = getDegreeDistribution(x2);
    cout << "Degree Distribution After Swapping:\n";
   printDegreeDistribution(degreeDistributionAfter);

    cout << endl << endl;
    return 0;
}
