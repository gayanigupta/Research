#include <iostream>
#include <list>
#include <set>
#include <map>
#include <array>

template<typename T>

/**
 * Making a generic Array Class which will be substitute vector in the program
 * 
*/


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
