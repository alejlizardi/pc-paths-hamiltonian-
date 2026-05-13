// C++ verification of Opus 4.7 session 2 counterexamples.
//
// Verifies in multiple independent ways:
//   C1: sigma' at n=18 -> rho=10 (disproves ALT)
//   C1b: sigma_4 at n=19 -> rho=11 (sanity check)
//   C2: stacked family [3,1,3,3,1]^k  k=0..4
//   C3: best-slope module
//   C7: rho >= b for test sigmas (the surviving conjecture)
//
// Uses the same bitmask DFS as the main verifier (well-tested).
//
// Compile: g++ -O3 -march=native -std=c++20 verify_opus_v2_counterexamples.cpp -o verify_opus_v2_counterexamples.exe

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <vector>
#include <chrono>

using u32 = uint32_t;
using u64 = uint64_t;

// DFS finding longest PC path (vertex count)
static int dfs_best;

static void dfs(int cur, u64 visited, int next_color, int depth,
                const u64* v_nbrs, const u64* p_nbrs, int n) {
    if (depth > dfs_best) dfs_best = depth;
    if (dfs_best >= n) return;
    u64 nbrs = (next_color == 0 ? p_nbrs[cur] : v_nbrs[cur]) & ~visited;
    u64 full = (n < 64) ? ((1ULL << n) - 1) : (u64)-1;
    while (nbrs) {
        int nb = __builtin_ctzll(nbrs);
        nbrs &= nbrs - 1;
        int can_add = __builtin_popcountll(~(visited | (1ULL << nb)) & full);
        if (depth + 1 + can_add <= dfs_best) continue;
        dfs(nb, visited | (1ULL << nb), 1 - next_color, depth + 1, v_nbrs, p_nbrs, n);
        if (dfs_best >= n) return;
    }
}

static int rho_f(int n, const int* sigma) {
    if (n > 64) return -1;
    std::vector<int> inv(n);
    for (int i = 0; i < n; i++) inv[sigma[i]] = i;
    std::vector<u64> v_nbrs(n, 0), p_nbrs(n, 0);
    for (int i = 0; i < n; i++) {
        int sv = sigma[i];
        if (sv > 0)   v_nbrs[i] |= 1ULL << inv[sv - 1];
        if (sv < n-1) v_nbrs[i] |= 1ULL << inv[sv + 1];
        if (i > 0)    p_nbrs[i] |= 1ULL << (i - 1);
        if (i < n-1)  p_nbrs[i] |= 1ULL << (i + 1);
    }
    dfs_best = 1;
    for (int start = 0; start < n; start++) {
        u64 vis = 1ULL << start;
        dfs(start, vis, 0, 1, v_nbrs.data(), p_nbrs.data(), n);
        if (dfs_best >= n) return dfs_best;
        dfs(start, vis, 1, 1, v_nbrs.data(), p_nbrs.data(), n);
        if (dfs_best >= n) return dfs_best;
    }
    return dfs_best;
}

// Count +-1-agreement blocks
static int count_blocks(int n, const int* sigma) {
    if (n <= 1) return n;
    int b = 1, cur_dir = 0;
    for (int i = 1; i < n; i++) {
        int diff = sigma[i] - sigma[i-1];
        if (diff == 1 || diff == -1) {
            if (cur_dir == 0) cur_dir = diff;
            else if (diff != cur_dir) { b++; cur_dir = diff; }
        } else {
            b++; cur_dir = 0;
        }
    }
    return b;
}

// Build a block-reversal permutation from composition
static void build_from_comp(const std::vector<int>& comp, std::vector<int>& sigma) {
    sigma.clear();
    int pos = 0;
    for (int s : comp) {
        for (int i = s - 1; i >= 0; i--) sigma.push_back(pos + i);
        pos += s;
    }
}

// Report for one sigma
static void report(const char* label, const std::vector<int>& sigma, int expected_rho = -1) {
    int n = (int)sigma.size();
    auto t0 = std::chrono::steady_clock::now();
    int r = rho_f(n, sigma.data());
    auto t1 = std::chrono::steady_clock::now();
    double dt = std::chrono::duration<double>(t1 - t0).count();

    int b = count_blocks(n, sigma.data());
    int d = n - r;
    int target = (n + 3 + 1) / 2; // ceil((n+3)/2)

    printf("  %s: n=%d rho=%d d=%d b=%d  time=%.3fs\n", label, n, r, d, b, dt);
    printf("    ALT target ceil((n+3)/2) = %d, %s\n",
           target, r >= target ? "ALT OK" : "*** ALT VIOLATED ***");
    printf("    rho >= b? %s (b=%d, rho=%d)\n",
           r >= b ? "YES" : "*** NO ***", b, r);
    printf("    d <= b-1? %s  d <= b? %s  d <= n-b? %s\n",
           d <= b-1 ? "yes" : "**NO**",
           d <= b   ? "yes" : "**NO**",
           d <= n-b ? "yes" : "**NO**");
    if (expected_rho >= 0) {
        printf("    Expected rho=%d: %s\n", expected_rho,
               r == expected_rho ? "MATCHES" : "*** MISMATCH ***");
    }
    printf("    sigma=[");
    for (int i = 0; i < n; i++) printf("%s%d", i?",":"", sigma[i]);
    printf("]\n\n");
}

int main() {
    printf("=== C++ Verification of Opus 4.7 Session 2 Counterexamples ===\n\n");

    // === TASK C1: The primary counterexample ===
    printf("--- Task C1: Primary counterexample (disproves ALT at n=18) ---\n");
    {
        std::vector<int> sigma = {1,0,4,3,2,5,8,7,6,11,10,9,12,15,14,13,17,16};
        report("sigma'", sigma, 10);
    }

    // === TASK C1b: Sanity check sigma_4 ===
    printf("--- Task C1b: Sanity check sigma_4 at n=19 ---\n");
    {
        // sigma_4: composition (2,3,1,3,1,3,1,3,2)
        std::vector<int> comp = {2,3,1,3,1,3,1,3,2};
        std::vector<int> sigma;
        build_from_comp(comp, sigma);
        report("sigma_4", sigma, 11);
    }

    // === Sanity check: sigma_2 and sigma_3 ===
    printf("--- Additional sanity checks ---\n");
    for (int m : {1, 2, 3, 4, 5}) {
        std::vector<int> comp = {2};
        comp.push_back(3);
        for (int t = 0; t < m-1; t++) { comp.push_back(1); comp.push_back(3); }
        comp.push_back(2);
        std::vector<int> sigma;
        build_from_comp(comp, sigma);
        int expected = 2*m + 3;
        char label[32]; snprintf(label, 32, "sigma_%d", m);
        report(label, sigma, expected);
    }

    // === TASK C2: Stacked family with module [3,1,3,3,1] ===
    printf("--- Task C2: Stacked family [2] + [3,1,3,3,1]^k + [3,2] ---\n");
    printf("    Module = (3,1,3,3,1), size 11, outer = (2) + [.] + (3,2) size 7\n\n");
    for (int k = 0; k <= 4; k++) {
        std::vector<int> comp = {2};
        for (int i = 0; i < k; i++) {
            comp.push_back(3); comp.push_back(1); comp.push_back(3);
            comp.push_back(3); comp.push_back(1);
        }
        comp.push_back(3); comp.push_back(2);
        std::vector<int> sigma;
        build_from_comp(comp, sigma);
        int n = (int)sigma.size();
        if (n > 60) continue; // bitmask limit
        char label[64]; snprintf(label, 64, "stack k=%d", k);
        int expected[] = {5, 10, 15, 20, 25, 30};
        report(label, sigma, k <= 5 ? expected[k] : -1);
    }

    // === TASK C2b: The specific n=25 counterexample ===
    printf("--- Task C2b: n=25 counterexample (disproves revised ALT, d<=b) ---\n");
    {
        std::vector<int> comp = {2,3,1,3,3,1,3,3,1,3,2};
        std::vector<int> sigma;
        build_from_comp(comp, sigma);
        report("n=25 bad", sigma, 13);
    }

    // === TASK C3: Best-slope module [3,3,3,3,3,1,3,3,3,3] ===
    printf("--- Task C3: Best module [3,3,3,3,3,1,3,3,3,3] (size 28) ---\n");
    for (int k = 1; k <= 2; k++) {
        std::vector<int> comp = {2};
        for (int i = 0; i < k; i++) {
            comp.push_back(3); comp.push_back(3); comp.push_back(3);
            comp.push_back(3); comp.push_back(3); comp.push_back(1);
            comp.push_back(3); comp.push_back(3); comp.push_back(3);
            comp.push_back(3);
        }
        comp.push_back(3); comp.push_back(2);
        std::vector<int> sigma;
        build_from_comp(comp, sigma);
        int n = (int)sigma.size();
        if (n > 60) { printf("    k=%d skipped (n=%d > 60)\n", k, n); continue; }
        char label[64]; snprintf(label, 64, "best_module k=%d", k);
        report(label, sigma);
    }

    printf("=== DONE ===\n");
    return 0;
}
