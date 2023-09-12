/***

1. Define a simple array of size V + E
2. if the number of v


***/
// define an array of the input file number of lines as the # of edges 
//find the max vertex
//array-size is number of vertices + edges 
// set the E-index as number of vertices
//set current node to -1
// read the input file each line
// if the verex one is not eql to curr_node
//assign current node E-index
//else assign array index the next vertex in the input file read line
// each line of the imnput file ha a start vertexa nd an end vertex
//increment E-index

#include <iostream>
#include <fstream>

struct Edge {
    int start;
    int end;
};

// Function to read the input file and process the data.
void readInputFile(const std::string& filename, int& numEdges, int* graph, int& numVertices) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the input file: " << filename << std::endl;
        return;
    }

   
    int EIndex = numVertices;
    int currentNode = -1;

    // Read each line of the input file.
    while (!inputFile.eof()) {
        int vertex1, vertex2;
        double weight;
        inputFile >> vertex1 >> vertex2 >> weight;
        std::cout<< "vertex1:  "<<vertex1<< std::endl;
         std::cout<< "vertex2: "<<vertex2<< std::endl;
         //std::cout<< vertex2<< std::endl;

        // If vertex1 is not equal to currentNode, update currentNode to vertex1.
        if (vertex1 != currentNode) {
            currentNode = vertex1;
            graph[currentNode] = EIndex;
            std::cout << currentNode <<"::"<< EIndex <<"\n";
        } 
            graph[EIndex] = vertex2;
      

        EIndex++;
    }
    for(int i =0 ;i <= numVertices+numEdges ; i++){
       std::cout << graph[i] << " "<< i << std::endl;
    }
    inputFile.close();
}

int main() {
    std::string filename ="/Users/gayanigupta/Documents/GitHub/Research/graph_decomposition _refined/Example/Tests/core_3.txt";
    int* graph = nullptr;
    int numVertices, numEdges;
    numVertices = 6;
    numEdges = 8;

    // Determine the number of edges by counting the lines in the input file.
    std::ifstream countFile(filename);
    if (countFile.is_open()) {
        numEdges = std::count(std::istreambuf_iterator<char>(countFile), std::istreambuf_iterator<char>(), '\n');
        countFile.close();
    }

    // Create an array to hold the graph data.
    graph = new int[numVertices + numEdges+2];

      for (int i = 0; i <= numVertices+numEdges ; ++i) {
        graph[i] = -1;
    }

    // Call the function to read and process the input file.
    readInputFile(filename, numEdges, graph, numVertices);

    // Print the number of vertices and edges.
    std::cout << "Number of Vertices: " << numVertices << std::endl;
    std::cout << "Number of Edges: " << numEdges << std::endl;

    // Clean up allocated memory.
    delete[] graph;

    return 0;
}
