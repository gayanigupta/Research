// C++ program to find K-Cores of a graph
#include <iostream>
#include <list>

// INPUT HEADERS
#include "translate_from_input.hpp"
#include "input_to_network.hpp"
#include "structure_defs.hpp"
#include "find_neighbors.hpp"

// OUTPUT HEADERS
#include "printout_network.hpp"
#include "printout_others.hpp"
using namespace std;

void heapify(vector<pair<int, float>>& arr, int n, int i)
{
	int largest = i; // Initialize largest as root Since we are using 0 based indexing
	int l = 2 * i + 1; // left = 2*i + 1
	int r = 2 * i + 2; // right = 2*i + 2

	// If left child is larger than root
	if (l < n && arr[l].second > arr[largest].second)
		largest = l;

	// If right child is larger than largest so far
	if (r < n && arr[r].second > arr[largest].second)
		largest = r;

	// If largest is not root
	if (largest != i) {
		swap(arr[i], arr[largest]);

		// Recursively heapify the affected sub-tree
		heapify(arr, n, largest);
	}
}

// main function to do heap sort
void heapSort(vector<pair<int, float>>& arr, int n)
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

void print_network(A_Network& x, const char* network_name)
{
	cout << "--" << network_name << "--" << endl;
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

int get_max_no_of_network(A_Network &x)
{
	int max_no = -1;
	for (int i = 0; i < x.size(); ++i)
	{
		vector<int_double> &listW = x[i].ListW;
		for (int j = 0; j < listW.size(); ++j)
		{
			if (listW[j].first > max_no)
				max_no = listW[j].first;
		}
	}
	return max_no;
}

bool checkEdge(ADJ_Bundle &bundle, int vertex_no)
{
	vector<int_double> &adj = bundle.ListW;
	for (int i = 0; i < adj.size(); ++i)
	{
		if (adj[i].first == vertex_no)
			return true;
	}

	return false;
}

bool addUndirectedEdge(ADJ_Bundle &bundle, int vertex_no)
{
	bool found = false;
	vector<int_double> &adj = bundle.ListW;
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

bool removeDirectedEdge(A_Network &x, int start_no, int end_no)
{
	if (start_no > x.size() - 1)
		return false;

	vector<int_double> &listW = x[start_no].ListW;
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

bool removeUndirectedEdge(A_Network &x, int start_no, int end_no)
{
	return removeDirectedEdge(x, start_no, end_no) && removeDirectedEdge(x, end_no, start_no);
}

void convertDirected2UndirectedNetwork(A_Network &a, A_Network &b)
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
		vector<int_double> &adj = a[i].ListW;
		for (int j = 0; j < adj.size(); ++j)
			addUndirectedEdge(b[adj[j].first], i);
	}
}

// Calculate the local degree coefficient of a vertex in undirected graph.
float calculateLocalDegreeCoefficient(A_Network& x, int vertex)
{
	int localLinkCount = 0;
	 
	vector<int_double>& neighbors = x[vertex].ListW;
	int neighborCount = neighbors.size();

	// For all neighbors of the vertex
	for (int i = 0; i < neighborCount; ++i)
	{
		int neighbor = neighbors[i].first;
		vector<int_double>& neighborsOfNeighbor = x[neighbor].ListW;

		// For all neighbors of the neighbor
		for (int j = 0; j < neighborsOfNeighbor.size(); ++j)
		{
			int neighborOfNeighbor = neighborsOfNeighbor[j].first;
			if (neighbor == neighborOfNeighbor)
				continue;

			// Count if the vertex is linked with neighor of its neighbor.
			for (int k = 0; k < neighborCount; ++k)
			{
				if (neighborOfNeighbor == neighbors[k].first)
					localLinkCount++;
			}

		}
	}

	// Return the local degree coefficient.
	return neighborCount < 2 ? 0 : (float)localLinkCount / (float)(neighborCount * (neighborCount - 1));
}

enum DEGREE_MODE
{
	PURE_DEGREE,    
	LOCAL_DEGREE_COEFFICIENT
};

// Run the elimination ordering of a graph with an specific degree mode.
void runEliminationOrdering(A_Network &src, vector<int> &orderList, vector<vector<int> > &tree, DEGREE_MODE degreeMode)
{
	A_Network x;
	convertDirected2UndirectedNetwork(src, x);
	print_network(x, "Undirected Network");

	int count = x.size();

	// Store degrees of all vertices
	vector<pair<int, float>> verticesWithDegree(count, pair<int, int>(0, 0));
	for (int i = 0; i < count; i++)
	{
		verticesWithDegree[i].first = i;
		// Get the degree value according to degree mode.
		verticesWithDegree[i].second = (degreeMode == PURE_DEGREE) ? x[i].ListW.size() : calculateLocalDegreeCoefficient(x, i);
	}

	// Sort according to degree of vertex.
	heapSort(verticesWithDegree, verticesWithDegree.size());

	// Mark all vertices not visited.
	vector<bool> processed(count, false);

	orderList.clear();
	tree.clear();

	int start = 0;

	// For all vertices of the graph with accending order
	while (start != count)
	{
		int sel_vertex = verticesWithDegree[start].first;
		// Check if the current vertex is processed.
		if (processed[sel_vertex])
		{
			start++;
			continue;
		}

		// Create a tree of the current vertex.
		vector<int> treeVertices;
		// Register oneself into the tree.
		treeVertices.push_back(sel_vertex);

		// Register all the neighbors of the current vertex into the tree
		// and mark them visited status.
		vector<int_double>& listW = x[sel_vertex].ListW;
		for (int i = 0; i < listW.size(); ++i)
		{
			int neighbor = listW[i].first;
			if (processed[neighbor])
				continue;

			treeVertices.push_back(neighbor);
			processed[neighbor] = true;
		}

		processed[sel_vertex] = true;
		orderList.push_back(sel_vertex);

		tree.push_back(treeVertices);
		start++;

		// print_network(x, "Undirected Network");
	}
}

// Driver program to test methods of graph class
int main(int argc, char *argv[])
{
	/*vector<Edge> edges;
	edges.push_back(create(0, 1, 1));
	edges.push_back(create(0, 2, 1));
	edges.push_back(create(1, 2, 1));
	edges.push_back(create(1, 5, 1));
	edges.push_back(create(2, 3, 1));
	edges.push_back(create(2, 4, 1));
	edges.push_back(create(2, 5, 1));
	edges.push_back(create(2, 6, 1));
	edges.push_back(create(3, 4, 1));
	edges.push_back(create(3, 6, 1));
	edges.push_back(create(3, 7, 1));
	edges.push_back(create(4, 6, 1));
	edges.push_back(create(4, 7, 1));
	edges.push_back(create(5, 6, 1));
	edges.push_back(create(5, 8, 1));
	edges.push_back(create(6, 7, 1));
	edges.push_back(create(6, 8, 1));*/

	clock_t q, q1, q2, t;
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

	

	// Create a network with edges.
	/*A_Network x1;
	create_Network(&edges, 0, &x1, -1);
	cout << "**** \n";*/
	

	// Preprocess Nodes to Numbers
	// Stores file in argv[3]: store map in argv[4]
	// Vertices start from 0
	q = clock();
	// Check if valid input is given
	if (argc < 3)
	{
		cout << "INPUT ERROR:: At least 2 inputs required. First: filename \n Second: Filetypes: 1:node_node_wt 2:node_wt_node 3:node_node 4:node_node (Only option 1 is active now) \n Third. Name of new file \n Fourth. Name of Map file\n";
		return 0;
	}
	// Check to see if file opening succeeded
	ifstream the_file(argv[1]);
	if (!the_file.is_open())
	{
		cout << "INPUT ERROR:: Could not open file\n";
	}

	A_Network x1;
	int nodes = -1;
	map_int_st revmap;
	int type = atoi("1");
	translate_input(argv[1], type, argv[3], argv[4]);

	// Remove Duplicate Edges and Self Loops; Create Undirected Graphs
	//  process_to_simple_undirected();
	q = clock() - q;
	cout << "Total Time for Preprocessing" << ((float)q) / CLOCKS_PER_SEC << "\n";

	/***** Preprocessing to Graph (GUI) ***********/

	/******* Read Graph (GUI) and Create Reverse Map*****************/
	// Obtain the list of edges.
	q = clock();
	readin_network(&x1, argv[3], nodes);

	// Create Reversemap

	nodes = x1.size();
	create_map(argv[4], &revmap);

	q = clock() - q;
	cout << "Total Time for Reading Network" << ((float)q) / CLOCKS_PER_SEC << "\n";


	vector<int> orderList;
	vector<vector<int> > tree;

	// Get the order list and its responding tree 
	// with an specific degree mode(PURE DEGREE or LOCAL_DEGREE_COEFFICIENT).
	//runEliminationOrdering(x1, orderList, tree, PURE_DEGREE);
	runEliminationOrdering(x1, orderList, tree, LOCAL_DEGREE_COEFFICIENT);

	// Output the result of elimination ordering.
	cout << "Tree decomposition result:" << endl;
	for (int i = 0; i < orderList.size(); ++i)
	{
		cout << orderList[i] << ":";
		for (int j = 0; j < tree[i].size(); ++j)
		{
			cout << tree[i][j];
			if (j < tree[i].size() - 1)
				cout << "-";
		}
		cout << endl;
	}

	return 0;
}
