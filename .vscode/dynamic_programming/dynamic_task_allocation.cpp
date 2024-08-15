#include <iostream>
#include <climits>

using namespace std;

// Function to find the minimum lateness using dynamic programming
int findMinLateness(const int t[], const int D[], int n) {
    // Calculate maximum possible time
    int max_time = 0;
    for (int i = 0; i < n; ++i) {
        max_time += t[i];
    }

    // Initialize DP table with infinity (or a very large number)
    int dp[n + 1][max_time + 1];
    for (int j = 0; j <= max_time; ++j) {
        dp[0][j] = 0; // Base case: No tasks result in zero lateness
    }
    for (int i = 1; i <= n; ++i) {
        for (int j = 0; j <= max_time; ++j) {
            dp[i][j] = INT_MAX;
        }
    }
    
    // Fill the DP table
    for (int i = 1; i <= n; ++i) {
        for (int j = 0; j <= max_time; ++j) {
            if (j >= t[i - 1]) {
                int lateness_if_scheduled = max(0, j - D[i - 1]);
                dp[i][j] = min(dp[i - 1][j], dp[i - 1][j - t[i - 1]] + lateness_if_scheduled);
            } else {
                dp[i][j] = dp[i - 1][j];
            }
        }
    }
    
    // The optimal solution is the minimum value in the last row of dp table
    int result = INT_MAX;
    for (int j = 0; j <= max_time; ++j) {
        result = min(result, dp[n][j]);
    }
    
    return result;
}

int main() {
    // Example usage
    const int n = 4;
    int completion_times[n] = {3, 2, 1, 4};
    int deadlines[n] = {4, 2, 2, 6};
    
    int optimal_lateness = findMinLateness(completion_times, deadlines, n);
    cout << "The optimal lateness is: " << optimal_lateness << endl;
    
    return 0;
}
