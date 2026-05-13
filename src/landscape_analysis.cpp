// Approach 3: Extremal characterization via local operations.
//
// For each sigma in S_n, compute:
//   1. rho(sigma) via DFS
//   2. For each adjacent transposition sigma' (swap sigma[i], sigma[i+1]):
//      compute rho(sigma') and delta = rho(sigma') - rho(sigma)
//   3. Check: is sigma a LOCAL MINIMUM of rho? (no adj transposition decreases rho)
//   4. For local minima: what's their structure? Are they all in the extremal family?
//
// Also: test whether "smoothing" a block boundary (merging two blocks by
// an adjacent transposition) always increases or preserves rho.
//
// Compile: g++ -O3 -march=native -std=c++20 -fopenmp analyze.cpp -o analyze.exe

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <mutex>
#include <atomic>
#include <omp.h>

using u32 = uint32_t;
using u64 = uint64_t;

static int N;

// --- DFS for longest alternating path (vertex count) ---

static thread_local int dfs_best;

static void dfs(int cur, u32 visited, int next_color, int depth,
                const u32* v_nbrs, const u32* p_nbrs, int n) {
    if (depth > dfs_best) dfs_best = depth;
    if (dfs_best >= n) return;

    u32 nbrs = (next_color == 0 ? p_nbrs[cur] : v_nbrs[cur]) & ~visited;
    u32 full = (1u << n) - 1;

    while (nbrs) {
        int nb = __builtin_ctz(nbrs);
        nbrs &= nbrs - 1;
        int can_add = __builtin_popcount(~(visited | (1u << nb)) & full);
        if (depth + 1 + can_add <= dfs_best) continue;
        dfs(nb, visited | (1u << nb), 1 - next_color, depth + 1,
            v_nbrs, p_nbrs, n);
        if (dfs_best >= n) return;
    }
}

static int longest_alt_path(int n, const int* sigma) {
    int inv[20];
    for (int i = 0; i < n; i++) inv[sigma[i]] = i;

    u32 v_nbrs[20] = {}, p_nbrs[20] = {};
    for (int i = 0; i < n; i++) {
        int sv = sigma[i];
        if (sv > 0)   v_nbrs[i] |= 1u << inv[sv - 1];
        if (sv < n-1) v_nbrs[i] |= 1u << inv[sv + 1];
        if (i > 0)    p_nbrs[i] |= 1u << (i - 1);
        if (i < n-1)  p_nbrs[i] |= 1u << (i + 1);
    }

    dfs_best = 1;
    for (int start = 0; start < n; start++) {
        u32 vis = 1u << start;
        dfs(start, vis, 0, 1, v_nbrs, p_nbrs, n);
        if (dfs_best >= n) return dfs_best;
        dfs(start, vis, 1, 1, v_nbrs, p_nbrs, n);
        if (dfs_best >= n) return dfs_best;
    }
    return dfs_best;
}

// --- Check if sigma is a local minimum under adjacent transpositions ---

struct LocalMinInfo {
    int sigma[20];
    int rho;
    int n_neighbors_equal;   // how many adj transpositions preserve rho
    int n_neighbors_higher;  // how many increase rho
    int min_neighbor_rho;    // smallest rho among neighbors
    int max_neighbor_rho;    // largest rho among neighbors
};

static LocalMinInfo check_local_min(int n, const int* sigma) {
    LocalMinInfo info{};
    memcpy(info.sigma, sigma, n * sizeof(int));
    info.rho = longest_alt_path(n, sigma);
    info.min_neighbor_rho = n;
    info.max_neighbor_rho = 0;

    int sigma2[20];
    memcpy(sigma2, sigma, n * sizeof(int));

    for (int i = 0; i < n - 1; i++) {
        // Swap positions i and i+1
        std::swap(sigma2[i], sigma2[i + 1]);
        int rho2 = longest_alt_path(n, sigma2);

        if (rho2 == info.rho) info.n_neighbors_equal++;
        else if (rho2 > info.rho) info.n_neighbors_higher++;
        // rho2 < info.rho means NOT a local minimum

        if (rho2 < info.min_neighbor_rho) info.min_neighbor_rho = rho2;
        if (rho2 > info.max_neighbor_rho) info.max_neighbor_rho = rho2;

        // Swap back
        std::swap(sigma2[i], sigma2[i + 1]);
    }

    return info;
}

// Is this a local minimum? (no neighbor has strictly lower rho)
static bool is_local_min(const LocalMinInfo& info, int n) {
    return info.min_neighbor_rho >= info.rho;
}

// --- Symmetry: canonical form under (rev_P, rev_Q, swap) ---

static void to_canonical(int n, const int* sigma, int* out) {
    int best[20];
    memcpy(best, sigma, n * sizeof(int));

    int tmp[20], tmp2[20];
    for (int mask = 1; mask < 8; mask++) {
        memcpy(tmp, sigma, n * sizeof(int));
        if (mask & 1) { // rev_P
            for (int i = 0; i < n; i++) tmp2[i] = tmp[n-1-i];
            memcpy(tmp, tmp2, n * sizeof(int));
        }
        if (mask & 2) { // rev_Q
            for (int i = 0; i < n; i++) tmp[i] = n - 1 - tmp[i];
        }
        if (mask & 4) { // swap (invert)
            for (int i = 0; i < n; i++) tmp2[tmp[i]] = i;
            memcpy(tmp, tmp2, n * sizeof(int));
        }
        bool smaller = false;
        for (int i = 0; i < n; i++) {
            if (tmp[i] < best[i]) { smaller = true; break; }
            if (tmp[i] > best[i]) break;
        }
        if (smaller) memcpy(best, tmp, n * sizeof(int));
    }
    memcpy(out, best, n * sizeof(int));
}

static std::string sigma_str(int n, const int* sigma) {
    std::string s = "[";
    for (int i = 0; i < n; i++) {
        if (i > 0) s += ",";
        s += std::to_string(sigma[i]);
    }
    s += "]";
    return s;
}

// --- Main ---

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }
    N = atoi(argv[1]);
    if (N < 3 || N > 14) {
        fprintf(stderr, "n must be 3..14\n");
        return 1;
    }

    printf("Local minimum analysis: n=%d\n", N);
    printf("Target rho: %d\n\n", (N + 3 + 1) / 2);

    // Phase 1: Find all sigmas achieving the global minimum rho
    std::atomic<int> global_min_rho{N};
    std::atomic<u64> total{0};

    // Store rho for each sigma (indexed by rank in lex order)
    // For n <= 12 this fits in memory (n! * 1 byte)
    // For n=12: 479M entries. Too much. Let's just find local minima.

    // Phase 1: Find global min rho
    printf("Phase 1: Finding global min rho...\n");

    #pragma omp parallel
    {
        int sigma[20];
        #pragma omp for schedule(dynamic, 1)
        for (int pair = 0; pair < N * (N - 1); pair++) {
            int first = pair / (N - 1);
            int second_raw = pair % (N - 1);
            int second = second_raw >= first ? second_raw + 1 : second_raw;

            sigma[0] = first;
            sigma[1] = second;

            int avail[20], n_avail = 0;
            for (int v = 0; v < N; v++)
                if (v != first && v != second) avail[n_avail++] = v;

            std::sort(avail, avail + n_avail);
            do {
                for (int i = 0; i < n_avail; i++)
                    sigma[i + 2] = avail[i];

                int rho = longest_alt_path(N, sigma);
                total.fetch_add(1, std::memory_order_relaxed);

                int cur_min = global_min_rho.load();
                while (rho < cur_min) {
                    if (global_min_rho.compare_exchange_weak(cur_min, rho))
                        break;
                }
            } while (std::next_permutation(avail, avail + n_avail));
        }
    }

    int min_rho = global_min_rho.load();
    printf("  Total: %llu, global min rho: %d (deficit %d)\n\n",
           (unsigned long long)total.load(), min_rho, N - min_rho);

    // Phase 2: Find all sigmas at or near global min, check local minimum
    printf("Phase 2: Finding local minima (rho <= %d)...\n", min_rho + 2);

    std::mutex results_mutex;
    std::vector<LocalMinInfo> local_minima;
    std::map<int, u64> rho_hist;
    std::set<std::string> seen_canonical;
    int n_checked = 0;

    #pragma omp parallel
    {
        int sigma[20];
        #pragma omp for schedule(dynamic, 1)
        for (int pair = 0; pair < N * (N - 1); pair++) {
            int first = pair / (N - 1);
            int second_raw = pair % (N - 1);
            int second = second_raw >= first ? second_raw + 1 : second_raw;

            sigma[0] = first;
            sigma[1] = second;

            int avail[20], n_avail = 0;
            for (int v = 0; v < N; v++)
                if (v != first && v != second) avail[n_avail++] = v;

            std::sort(avail, avail + n_avail);
            do {
                for (int i = 0; i < n_avail; i++)
                    sigma[i + 2] = avail[i];

                int rho = longest_alt_path(N, sigma);

                if (rho <= min_rho + 2) {
                    // Check if local minimum
                    LocalMinInfo info = check_local_min(N, sigma);
                    bool is_lm = is_local_min(info, N);

                    std::lock_guard<std::mutex> lock(results_mutex);
                    rho_hist[rho]++;

                    if (is_lm && rho <= min_rho + 1) {
                        // Deduplicate by canonical form
                        int canon[20];
                        to_canonical(N, sigma, canon);
                        std::string key = sigma_str(N, canon);
                        if (seen_canonical.find(key) == seen_canonical.end()) {
                            seen_canonical.insert(key);
                            local_minima.push_back(info);
                        }
                    }
                }
            } while (std::next_permutation(avail, avail + n_avail));
        }
    }

    // Sort local minima by rho
    std::sort(local_minima.begin(), local_minima.end(),
              [](const LocalMinInfo& a, const LocalMinInfo& b) {
                  return a.rho < b.rho;
              });

    printf("\nLocal minima (canonical, rho <= %d):\n", min_rho + 1);
    printf("  %4s  %-40s  %4s %6s %6s\n",
           "rho", "sigma", "eq", "higher", "min_nb");

    for (auto& lm : local_minima) {
        printf("  %4d  %-40s  %4d %6d %6d\n",
               lm.rho,
               sigma_str(N, lm.sigma).c_str(),
               lm.n_neighbors_equal,
               lm.n_neighbors_higher,
               lm.min_neighbor_rho);
    }

    printf("\n  Total distinct canonical local minima: %d\n",
           (int)local_minima.size());

    // Summary
    printf("\n  Rho histogram (for rho <= %d):\n", min_rho + 2);
    for (auto& [r, c] : rho_hist) {
        printf("    rho=%d: %llu sigmas\n", r, (unsigned long long)c);
    }

    return 0;
}
