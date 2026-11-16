// Exact all-sources BFS on a connected, unweighted, undirected graph.
// Input: n m, followed by m pairs of zero-based node indices.
// Output: average distance, diameter, ordered distance sum, reachable pair count.
// Integer reductions make results independent of the OpenMP scheduling order.

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>
#include <omp.h>

int main() {
    int n, m;
    if (!(std::cin >> n >> m) || n < 0 || m < 0) {
        return 1;
    }
    std::vector<std::vector<int>> adjacency(n);
    for (int i = 0; i < m; ++i) {
        int a, b;
        if (!(std::cin >> a >> b) || a == b || a < 0 || b < 0 || a >= n || b >= n) {
            return 2;
        }
        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
    }

    long long distance_sum = 0;
    long long reachable_pairs = 0;
    int diameter = 0;

    #pragma omp parallel reduction(+:distance_sum,reachable_pairs) reduction(max:diameter) num_threads(4)
    {
        std::vector<int> distance(n);
        std::vector<int> queue(n);

        #pragma omp for schedule(dynamic,8)
        for (int source = 0; source < n; ++source) {
            std::fill(distance.begin(), distance.end(), -1);
            distance[source] = 0;
            int head = 0;
            int tail = 0;
            queue[tail++] = source;

            while (head < tail) {
                const int u = queue[head++];
                for (int v : adjacency[u]) {
                    if (distance[v] >= 0) {
                        continue;
                    }
                    distance[v] = distance[u] + 1;
                    queue[tail++] = v;
                    distance_sum += distance[v];
                    diameter = std::max(diameter, distance[v]);
                }
            }
            reachable_pairs += tail - 1;
        }
    }

    if (reachable_pairs != 1LL * n * (n - 1)) {
        std::cerr << "Expected connected input\n";
        return 3;
    }
    const double average = reachable_pairs ? double(distance_sum) / reachable_pairs : 0.0;
    std::cout << std::setprecision(15) << average << " " << diameter << " "
              << distance_sum << " " << reachable_pairs << "\n";
}
