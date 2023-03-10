// C++ program to find K-Cores of a graph
#include <iostream>
#include<list>
#include <algorithm>

// INPUT HEADERS
#include "translate_from_input.hpp"
#include "input_to_network.hpp"
#include "structure_defs.hpp"
#include "find_neighbors.hpp"

// OUTPUT HEADERS
#include "printout_network.hpp"
#include "printout_others.hpp"


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

void find_paths(A_Network& x, int src, int dst, int& path_count, vector<bool>& visited, vector<vector<int> >& paths)
{
	visited[src] = true;

	if (src == dst) {
		path_count++;
	} else {
		vector<int_double>& neighbours = x[src].ListW;

		int paths_len = paths.size();
		for (int i = 0; i < neighbours.size(); ++i) {
			int neighbour = neighbours[i].first;
			if (!visited[neighbour]) {
				for (int j = 0; j < paths_len; ++j) {
					int path_len = paths[j].size();
					if (path_len > 0 && paths[j][path_len - 1] == src) {
						vector<int>& path = paths[j];
						auto it = std::find(path.begin(), path.end(), neighbour);
						if (it != path.end()) {
							// Prevent for loop path
							continue;
						}
						vector<int> new_path = path;
						new_path.push_back(neighbour);
						paths.push_back(new_path);
					}
				}
			}
		}

		for (int i = 0; i < neighbours.size(); ++i) {
			int neighbour = neighbours[i].first;
			if (!visited[neighbour]) {
				find_paths(x, neighbour, dst, path_count, visited, paths);
			}
		}
	}

	visited[src] = false;
}

int browse_paths(A_Network& x, int src, int dst, int vertices, vector<vector<int> >& paths)
{
	vector<bool> visited(vertices, false);
	int path_count = 0;

	vector<int> init_path;
	init_path.push_back(src);
	paths.push_back(init_path);

	find_paths(x, src, dst, path_count, visited, paths);

	vector<vector<int> > result_paths;
	for (auto path : paths) {
		if (path.size() > 0 && path[path.size() - 1] == dst) {
			result_paths.push_back(path);
		}
	}

	paths = result_paths;

	// Sort the vector
	std::sort(paths.begin(), paths.end());

	// Remove duplicates
	auto it = std::unique(paths.begin(), paths.end());
	paths.erase(it, paths.end());

	return path_count;
}

// Driver program to test methods of graph class
int main(int argc, char *argv[])
{
	/*vector<Edge> edges;

	// add the edges.
	edges.push_back(create(0, 1, 1));
	edges.push_back(create(0, 2, 1));
	edges.push_back(create(0, 4, 1));
	edges.push_back(create(1, 3, 1));
	edges.push_back(create(1, 4, 1));
	edges.push_back(create(2, 3, 1));
	edges.push_back(create(2, 1, 1));
	edges.push_back(create(3, 2, 1));

	// Create a network with edges.
	A_Network x1;
	create_Network(&edges, 0, &x1, -1);
	print_network(x1, "Original Network");
*/

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

	/*

	// Create a network with edges.
	A_Network x1;
	create_Network(&edges, 0, &x1, -1);
	cout << "**** \n";
	*/

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

	// Get the max no of vertex in graph.
	int count = get_max_no_of_network(x1) + 1;

	// Getn all the paths between two nodes and print them
	int start = 0;
	int end = 4;
	vector<vector<int> > paths;

	cout << "------All pathes from " << start << " node to " << end << " node:" << endl;
	cout << "   Count: " << browse_paths(x1, start, end, count, paths) << endl;

	cout << "   Paths:" << endl;
	for (auto path : paths) {
		cout << "      ";
		for (auto node : path) {
			cout << node << " > ";
		}
		cout << endl;
	}

	cout << endl << endl;
	return 0;
}

