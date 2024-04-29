std::pair<std::vector<Edge>, std::vector<int>> findChordalEdgesWithEliminationOrder(std::vector<Edge>& edges, int numVertices)
{
    std::vector<Edge> chordalEdges(edges.size());
    std::vector<int> eliminationOrder(numVertices);
    std::map<int, int> visited;
    std::map<int, std::vector<int>> adjList;
    std::vector<std::vector<int>> setPartition;

    std::vector<int> initSet(numVertices);

    int setPartSize = 0;
    for (int v = 0; v < numVertices; v++)
    {
        getNeighbours(v, edges, adjList, numVertices);
        visited[v] = 0;
        initSet.push_back(v);
    }

    for (int i = 0; i < edges.size(); i++){
        Edge edge = edges[i];
        chordalEdges.push_back(edge);
    }

    setPartition.push_back(initSet);
    setPartSize += 1;

    while (setPartSize != 0){
        int u = getLeftMostElem(setPartition, setPartSize);
        visited[u] = 1;
        eliminationOrder.push_back(u);

        std::vector<int> neighbours;
        for (int i = 0; i < adjList[u].size(); i++){
            int neighbour = adjList[u][i];
            if (visited[neighbour] == 0){
                neighbours.push_back(neighbour);
            }
        }

        updateSet(neighbours, setPartition, setPartSize, numVertices);
    }

    std::reverse(eliminationOrder.begin(), eliminationOrder.end());

    std::cout << "Elimination order below\n";
    for (int i = 0; i < eliminationOrder.size(); i++)
    {
        std::cout << eliminationOrder[i] << std::endl;
    }

    return std::make_pair(chordalEdges, eliminationOrder);
}
