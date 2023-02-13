// C++ program to find K-Cores of a graph
#include <iostream>
#include<list>

//INPUT HEADERS
#include "translate_from_input.hpp"
#include "input_to_network.hpp"
#include"structure_defs.hpp"
#include "find_neighbors.hpp"

//OUTPUT HEADERS
#include "printout_network.hpp"
#include "printout_others.hpp"
using namespace std;

class Graph;

void print_network(A_Network& x, const char* network_name)
{
	cout << "--" << network_name << "--" << endl;
	for (int i = 0; i < x.size(); ++i)
	{
		ADJ_Bundle& adj = x[i];

		cout << "Core " << i << "'s vertexs: ";
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

void convert_network_2_graph(A_Network& x, Graph& g);

class Graph
{
public:
	int V; // No. of vertices

	// Pointer to an array containing adjacency lists
	vector< vector<int> > adj;

	Graph(int V);
	Graph(A_Network& x); // Constructor

	// function to add an edge to graph
	void addEdge(int u, int v);

	// A recursive function to print DFS starting from v
	bool DFSUtil(int, vector<bool>&, vector<int>&, int k);

	// generate k-Cores degree information of given graph
	vector<int>  genKCores(int k);

	// prints k-Cores of given graph
	void printKCores(int k);
};

bool Graph::DFSUtil(int v, vector<bool>& visited, vector<int>& vDegree, int k)
{
	// Mark the current node as visited and print it
	visited[v] = true;

	// Recur for all the vertices adjacent to this vertex
	vector<int>& links = adj[v];

	for (int i = 0; i < links.size(); ++i)
	{
		// degree of v is less than k, then degree of adjacent
		// must be reduced
		if (vDegree[v] < k)
			vDegree[links[i]]--;

		// If adjacent is not processed, process it
		if (!visited[links[i]])
		{
			// If degree of adjacent after processing becomes
			// less than k, then reduce degree of v also.
			DFSUtil(links[i], visited, vDegree, k);
		}
	}

	// Return true if degree of v is less than k
	return (vDegree[v] < k);
}

Graph::Graph(int V)
{
	this->V = V;
	vector<int> links;
	for (int i = 0; i < V; ++i)
	{
		adj.push_back(links);
	}
}

Graph::Graph(A_Network& x)
{
	convert_network_2_graph(x, *this);
}

void Graph::addEdge(int u, int v)
{
	adj[u].push_back(v);
	adj[v].push_back(u);
}

vector<int> Graph::genKCores(int k)
{
	// INITIALIZATION
	// Mark all the vertices as not visited and not processed.
	vector<bool> visited(V, false);
	vector<bool> processed(V, false);

	int mindeg = INT_MAX;
	int startvertex;

	// Store degrees of all vertices
	vector<int> vDegree(V);
	for (int i = 0; i < V; i++)
	{
		vDegree[i] = adj[i].size();

		if (vDegree[i] < mindeg)
		{
			mindeg = vDegree[i];
			startvertex = i;
		}
	}

	DFSUtil(startvertex, visited, vDegree, k);

	// If Graph is disconnected.
	for (int i = 0; i < V; i++)
		if (visited[i] == false)
			DFSUtil(i, visited, vDegree, k);

	// Considering Edge Case (Example 3 in main() function)
	for (int v = 0; v < V; v++) {
		if (vDegree[v] >= k) {
			int cnt = 0;
			vector<int>& links = adj[v];
			for (int i = 0; i < links.size(); ++i)
				if (vDegree[links[i]] >= k)
					cnt++;
			if (cnt < k)
				vDegree[v] = cnt;
		}
	}

	return vDegree;
}

// Prints k cores of an undirected graph
void Graph::printKCores(int k)
{
	vector<int> vDegree = genKCores(k);

	cout << "K-Cores : \n";
	for (int v = 0; v < V; v++)
	{
		// Only considering those vertices which have degree
		// >= K after DFS
		if (vDegree[v] >= k)
		{
			cout << "\n[" << v << "]";

			// Traverse adjacency list of v and print only
			// those adjacent which have vDegree >= k after
			// DFS.
			vector<int>& links = adj[v];
			for (int i = 0; i < links.size(); ++i)
				if (vDegree[links[i]] >= k)
					cout << " -> " << links[i];
		}
	}
}

void convert_network_2_graph(A_Network& x, Graph& g)
{
	// Generate a graph from a network.
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

	g.V = max_no + 1;

	vector<int> links;
	for (int i = 0; i < g.V; ++i)
	{
		g.adj.push_back(links);
	}

	for (int i = 0; i < x.size(); ++i)
	{
		vector <int_double>& listW = x[i].ListW;
		for (int j = 0; j < listW.size(); ++j)
			g.addEdge(i, listW[j].first);
	}
}

void convert2NetworkWithKCores(A_Network& x, int k)
{
	Graph g(x);
	vector<int> vDegree = g.genKCores(k);
	
	for (int i = 0; i < x.size(); ++i)
	{
		if (vDegree[i] < k)
		{
			x[i].ListW.clear();
			x[i].Ops.clear();
			continue;
		}

		vector <int_double>& listW = x[i].ListW;
		vector <int_double> newlistW;

		for (int j = 0; j < listW.size(); j++)
		{
			if (vDegree[listW[j].first] >= k)
			{
				newlistW.push_back(listW[j]);
			}
		}

		x[i].ListW = newlistW;
	}
}


// Driver program to test methods of graph class
int main()
{
	A_Network x;

	//list of all edges of network
	vector<Edge> edges;
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
	edges.push_back(create(6, 8, 1));

	create_Network(&edges, 0, &x, -1);

	print_network(x, "Original Network");

	convert2NetworkWithKCores(x, 3);

	print_network(x, "Network with only K-cores");

	cout << endl << endl;

	return 0;
}

