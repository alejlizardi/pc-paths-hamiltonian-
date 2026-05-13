// Approach 2: Monotone block contribution analysis.
//
// For each sigma in S_n, compute:
//   1. The monotone block decomposition (maximal ascending/descending runs)
//   2. Block sizes and types (ascending/descending)
//   3. Boundary types (smooth = V-edge across boundary, break = no V-edge)
//   4. rho(sigma) via DFS
//   5. deficit = n - rho
//
// Then test hypotheses about per-block contribution to rho.
//
// Compile: g++ -O3 -march=native -std=c++20 -fopenmp analyze.cpp -o analyze.exe

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <map>
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

// --- Block decomposition ---

struct BlockInfo {
    int n_blocks;
    int sizes[20];       // block sizes
    int dirs[20];        // +1 = ascending, -1 = descending, 0 = singleton
    int n_breaks;        // boundaries without V-edge
    int n_smooth;        // boundaries with V-edge (parallel edges)
    int max_block;       // largest block size
    int n_desc_blocks;   // number of descending blocks (size >= 2)
    int n_singletons;    // number of singleton blocks
};

static BlockInfo decompose_blocks(int n, const int* sigma) {
    BlockInfo bi{};
    int i = 0;
    while (i < n) {
        int start = i;
        if (i + 1 < n && sigma[i + 1] > sigma[i]) {
            // ascending run
            while (i + 1 < n && sigma[i + 1] > sigma[i]) i++;
            i++;
            bi.sizes[bi.n_blocks] = i - start;
            bi.dirs[bi.n_blocks] = +1;
        } else if (i + 1 < n && sigma[i + 1] < sigma[i]) {
            // descending run
            while (i + 1 < n && sigma[i + 1] < sigma[i]) i++;
            i++;
            bi.sizes[bi.n_blocks] = i - start;
            bi.dirs[bi.n_blocks] = -1;
            bi.n_desc_blocks++;
        } else {
            // singleton
            i++;
            bi.sizes[bi.n_blocks] = 1;
            bi.dirs[bi.n_blocks] = 0;
            bi.n_singletons++;
        }
        if (bi.sizes[bi.n_blocks] > bi.max_block)
            bi.max_block = bi.sizes[bi.n_blocks];
        bi.n_blocks++;
    }

    // Analyze boundaries
    // A boundary between positions j and j+1 is "smooth" if there's a V-edge
    // across it, i.e., |sigma(j) - sigma(j+1)| = 1
    int pos = 0;
    for (int b = 0; b < bi.n_blocks - 1; b++) {
        pos += bi.sizes[b];
        int j = pos - 1; // last position of block b
        // j+1 = first position of block b+1
        if (j + 1 < n && abs(sigma[j] - sigma[j + 1]) == 1) {
            bi.n_smooth++;
        } else {
            bi.n_breaks++;
        }
    }

    return bi;
}

// Encode block composition as a string for grouping
static std::string block_key(const BlockInfo& bi) {
    std::string s;
    for (int b = 0; b < bi.n_blocks; b++) {
        if (b > 0) s += ",";
        s += std::to_string(bi.sizes[b]);
        s += (bi.dirs[b] > 0 ? "+" : (bi.dirs[b] < 0 ? "-" : "s"));
    }
    return s;
}

// Just sizes, ignoring direction
static std::string size_key(const BlockInfo& bi) {
    std::string s;
    for (int b = 0; b < bi.n_blocks; b++) {
        if (b > 0) s += ",";
        s += std::to_string(bi.sizes[b]);
    }
    return s;
}

// --- Main ---

struct StatsKey {
    int n_blocks, n_breaks, max_block;
    bool operator<(const StatsKey& o) const {
        if (n_blocks != o.n_blocks) return n_blocks < o.n_blocks;
        if (n_breaks != o.n_breaks) return n_breaks < o.n_breaks;
        return max_block < o.max_block;
    }
};

struct StatsVal {
    int max_deficit = 0;
    u64 count = 0;
};

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n> [--compositions]\n", argv[0]);
        return 1;
    }
    N = atoi(argv[1]);
    if (N < 3 || N > 16) {
        fprintf(stderr, "n must be 3..16\n");
        return 1;
    }

    bool show_compositions = (argc > 2 && strcmp(argv[2], "--compositions") == 0);

    printf("Block contribution analysis: n=%d\n\n", N);

    std::map<StatsKey, StatsVal> stats;
    std::mutex stats_mutex;
    std::atomic<int> global_max_deficit{0};
    std::atomic<u64> total{0};

    // Track by block SIZE composition
    std::map<std::string, StatsVal> comp_stats;
    std::mutex comp_mutex;

    // Track worst sigmas
    std::mutex worst_mutex;
    struct WorstEntry {
        int sigma[20];
        int deficit;
        BlockInfo bi;
    };
    std::vector<WorstEntry> worst_list;

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

                BlockInfo bi = decompose_blocks(N, sigma);
                int rho = longest_alt_path(N, sigma);
                int deficit = N - rho;

                total.fetch_add(1, std::memory_order_relaxed);

                {
                    StatsKey key{bi.n_blocks, bi.n_breaks, bi.max_block};
                    std::lock_guard<std::mutex> lock(stats_mutex);
                    auto& s = stats[key];
                    s.count++;
                    if (deficit > s.max_deficit) s.max_deficit = deficit;
                }

                if (show_compositions) {
                    std::string sk = size_key(bi);
                    std::lock_guard<std::mutex> lock(comp_mutex);
                    auto& s = comp_stats[sk];
                    s.count++;
                    if (deficit > s.max_deficit) s.max_deficit = deficit;
                }

                int old_max = global_max_deficit.load();
                if (deficit >= old_max) {
                    std::lock_guard<std::mutex> lock(worst_mutex);
                    if (deficit > global_max_deficit.load()) {
                        worst_list.clear();
                        global_max_deficit.store(deficit);
                    }
                    if (deficit >= global_max_deficit.load() &&
                        worst_list.size() < 20) {
                        WorstEntry we;
                        memcpy(we.sigma, sigma, N * sizeof(int));
                        we.deficit = deficit;
                        we.bi = bi;
                        worst_list.push_back(we);
                    }
                }
            } while (std::next_permutation(avail, avail + n_avail));
        }
    }

    printf("Total sigmas: %llu\n", (unsigned long long)total.load());
    printf("Max deficit: %d\n", global_max_deficit.load());
    printf("Target: %d\n\n", N - (N + 3 + 1) / 2);

    printf("Stats by (n_blocks, n_breaks, max_block):\n");
    printf("  %6s %6s %6s %10s %6s\n", "blocks", "breaks", "maxblk", "count", "maxdef");
    for (auto& [key, val] : stats) {
        if (val.max_deficit > 0 || val.count > 1000)
            printf("  %6d %6d %6d %10llu %6d\n",
                   key.n_blocks, key.n_breaks, key.max_block,
                   (unsigned long long)val.count, val.max_deficit);
    }

    if (show_compositions) {
        printf("\nBlock size compositions with max deficit > 0:\n");
        printf("  %-30s %10s %6s\n", "composition", "count", "maxdef");
        for (auto& [key, val] : comp_stats) {
            if (val.max_deficit > 0)
                printf("  %-30s %10llu %6d\n",
                       key.c_str(), (unsigned long long)val.count,
                       val.max_deficit);
        }
    }

    printf("\nWorst sigmas (deficit=%d):\n", global_max_deficit.load());
    for (auto& we : worst_list) {
        printf("  sigma=[");
        for (int i = 0; i < N; i++) printf("%s%d", i ? "," : "", we.sigma[i]);
        printf("] blocks=%d breaks=%d smooth=%d maxblk=%d singletons=%d ",
               we.bi.n_blocks, we.bi.n_breaks, we.bi.n_smooth,
               we.bi.max_block, we.bi.n_singletons);
        // Show block structure
        printf("structure=[");
        for (int b = 0; b < we.bi.n_blocks; b++) {
            if (b > 0) printf("|");
            printf("%d%s", we.bi.sizes[b],
                   we.bi.dirs[b] > 0 ? "+" : (we.bi.dirs[b] < 0 ? "-" : ""));
        }
        printf("]\n");
    }

    return 0;
}
