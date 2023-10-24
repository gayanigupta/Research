#include <iostream>
#include <fstream>


// Function to find chordal edges with elimination order.
void findChordalEdgesWithEliminationOrder(int * graph, int numEdges, int numVertices) {
    std::cout<<"I am in function"<<std::endl;
    

    // Step 6: Initialize the first set in setPartition.
    // Implement Step 6 here.

    // Step 7: While setPartSize is not equal to 0:
    // a. Get the leftmost element from the leftmost set in setPartition.
    // b. Mark u as visited and add to eliminationOrder.
    // c. Check u's neighbors and update sets.
    // Implement Step 7 here.

    // Step 8: Reverse the elements in eliminationOrder.
    // Implement Step 8 here.

    // Step 9: Print the elements of the elimination order.
    // Implement Step 9 here.

    // Step 10: Clean up allocated memory.

    // Step 11: Return the result as a pair.
    
    // Clean up allocated memory.
    delete[] graph;

    // Return the result as a pair.
    return;
}

/*int main() {
    // Example usage:
    Edge edges[] = {{0, 1}, {1, 2}, {2, 3}, {1, 3}, {3, 4}};
    int numVertices = 5;
    int numEdges = sizeof(edges) / sizeof(edges[0]);

    std::pair<Edge*, int*> result = findChordalEdgesWithEliminationOrder(edges, numEdges, numVertices);

    // Step 9: Print elimination order
    std::cout << "Elimination order below:\n";
    for (int i = 0; i < numVertices; i++) {
        std::cout << result.second[i] << " ";
    }
    std::cout << std::endl;

    // Step 10: Clean up allocated memory for chordalEdges and eliminationOrder.
    delete[] result.first;
    delete[] result.second;

    return 0;
}*/

