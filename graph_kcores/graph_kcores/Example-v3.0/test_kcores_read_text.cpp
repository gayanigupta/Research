// C++ program to find K-Cores of a graph
#include <iostream>
#include<list>

//INPUT HEADERS
#include "input_to_network.hpp"
#include"structure_defs.hpp"

void print_network(A_Network& x, const char* network_name)
{
	cout << "--" << network_name << "--" << endl;
	for (int i = 0; i < x.size(); ++i)
	{
		ADJ_Bundle& adj = x[i];

		if (adj.ListW.size()>0)
		{
		cout <<"  Core  : "<<adj.ListW.size()<<" - Vertex " << i << "'s - Linked vertices : ";
		for (int j = 0; j < adj.ListW.size(); ++j)
		{
			cout << adj.ListW.at(j).first;
			if (j < adj.ListW.size() - 1)
				cout << ", ";
		}
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

void print_network_kcores(A_Network& x, int count, int k, const char* network_name)
{
	cout << "--" << network_name << "' kcore group--" << endl;

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

	cout << "KCore " << k << ": ";
	for (int i = 1; i < count; ++i)
	{
		if (groups[i].size() == 0)
			continue;

		for (int j = 0; j < groups[i].size(); ++j)
			cout << groups[i][j] << ", ";
	}
	cout << endl;
}

void convert2NetworkWithKCoresEx(A_Network& src, A_Network& x, int k)
{
	convertDirected2UndirectedNetwork(src, x);

	int count = x.size();
	
	cout << endl;
	// Store degrees of all vertices
	vector<int> degrees(count, 0);
	for (int i = 0; i < count; i++)
	{
		degrees[i] = x[i].ListW.size();
		cout << "Vertex " << i << "'s degree: " << degrees[i] << endl;
	}

	vector<bool> processed(count, false);
	while (true)
	{
		int removed_count = 0;
		for (int i = 0; i < count; i++)
		{
			if (processed[i])
				continue;
			
			if (degrees[i] < k)
			{// Remove the vertex
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
			}
		}

		if (removed_count == 0)
			break;
	}

	cout << endl;
	for (int i = 0; i < count; i++)
		cout << "Vertex " << i << "'s degree: " << degrees[i] << endl;
}

// Driver program to test methods of graph class
int main()
{
	fstream fin;
	fin.open("Tests/core.txt", ios::in);
	if (fin.is_open() == false)
		return 0;

	vector<Edge> edges;
	while (!fin.eof())
	{
		int start_no, end_no, weight;
		fin >> start_no >> end_no >> weight;
		edges.push_back(create(start_no, end_no, 1));
	}

	fin.close();

	/*edges.push_back(create(0, 1, 1));
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

	A_Network x1;
	create_Network(&edges, 0, &x1, -1);

	print_network(x1, "Original Network");

	int count = get_max_no_of_network(x1) + 1;

	A_Network kcores;
	convert2NetworkWithKCoresEx(x1, kcores, 3);

	print_network(kcores, "Network with only K-cores");
	
	print_network_kcores(kcores, count, 3, "Network with only K-cores");

	cout << endl << endl;

	return 0;
}

