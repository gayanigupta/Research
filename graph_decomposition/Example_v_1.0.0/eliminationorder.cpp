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
	return max_no;
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

void runStandardEliminationOrdering(A_Network& src, vector <int>& orderList, vector<vector <int> >& tree)
{
	A_Network x;
	convertDirected2UndirectedNetwork(src, x);
	//print_network(x, "Undirected Network");

	int count = x.size();
	// Store degrees of all vertices
	vector<int> degrees(count, 0);
	for (int i = 0; i < count; i++)
		degrees[i] = x[i].ListW.size();

	vector<bool> processed(count, false);

	orderList.clear();
	tree.clear();
	while (orderList.size() != count)
	{
		// Get the vertex with minimum size of its neighbors
		int min_degree = INT_MAX;
		int sel_vertex = -1;
		for (int i = 0; i < count; i++)
		{
			if (processed[i])
				continue;
			if (min_degree > degrees[i])
			{
				min_degree = degrees[i];
				sel_vertex = i;
			}
		}

		if (sel_vertex == -1) // error
			break;

		vector<int> treeVertices;
		treeVertices.push_back(sel_vertex);

		degrees[sel_vertex] = 0;
		for (int i = 0; i < count; ++i)
		{
			if (processed[i])
				continue;

			if (removeUndirectedEdge(x, sel_vertex, i))
			{
				degrees[i]--;
				treeVertices.push_back(i);
			}
		}

		processed[sel_vertex] = true;
		orderList.push_back(sel_vertex);

		tree.push_back(treeVertices);		
		
		//print_network(x, "Undirected Network");

		sel_vertex = 0;
	}
}


// Driver program to test methods of graph class
int main()
{
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

	A_Network x1;
	create_Network(&edges, 0, &x1, -1);

	vector <int> orderList;
	vector<vector <int> > tree;
	runStandardEliminationOrdering(x1, orderList, tree);

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

	cout << endl << endl;

	return 0;
}

