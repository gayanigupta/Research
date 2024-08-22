#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>

using namespace std;

// Function to calculate maximum lateness considering dependencies
int calculateMaxLateness(const vector<int>& completionTimes, const vector<int>& deadlines, const vector<vector<int>>& dependencies) {
    int n = completionTimes.size();

    // Initialize DP and inDegree vectors
    vector<int> dp(n, 0);
    vector<int> inDegree(n, 0);

    // Build the graph and compute in-degrees
    vector<vector<int>> graph(n);
    for (int u = 0; u < n; ++u) {
        for (int v : dependencies[u]) {
            graph[v].push_back(u);
            inDegree[u]++;
        }
    }

    // Perform topological sort
    vector<int> topoOrder;
    vector<int> queue;
    for (int i = 0; i < n; ++i) {
        if (inDegree[i] == 0) {
            queue.push_back(i);
        }
    }

    while (!queue.empty()) {
        int u = queue.back();
        queue.pop_back();
        topoOrder.push_back(u);

        for (int v : graph[u]) {
            inDegree[v]--;
            if (inDegree[v] == 0) {
                queue.push_back(v);
            }
        }
    }

    // Initialize dp array
    for (int i = 0; i < n; ++i) {
        dp[i] = completionTimes[i];
    }

    // Compute the finish times for each task
    for (int u : topoOrder) {
        for (int v : graph[u]) {
            dp[v] = max(dp[v], dp[u] + completionTimes[v]);
        }
    }

    // Calculate maximum lateness
    int maxLateness = 0;
    for (int i = 0; i < n; ++i) {
        int lateness = max(0, dp[i] - deadlines[i]);
        maxLateness = max(maxLateness, lateness);
    }

    return maxLateness;
}

int main() {
    // Example usage
    vector<int> completionTimes = {3, 2, 4, 1};  // Completion times of tasks
    vector<int> deadlines = {4, 5, 6, 3};        // Deadlines of tasks

    // Define dependencies: task[i] depends on tasks in dependencies[i]
    vector<vector<int>> dependencies = {
        {},       // Task 0 (No dependencies)
        {},       // Task 1 (No dependencies)
        {0, 1},   // Task 2 depends on tasks 0 and 1
        {2}        // Task 3 depends on task 2
    };

    int maxLateness = calculateMaxLateness(completionTimes, deadlines, dependencies);
    cout << "The maximum lateness is: " << maxLateness << endl;

    return 0;
}
