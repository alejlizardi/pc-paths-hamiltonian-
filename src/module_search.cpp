// Agenda item 16: Search for exotic worst cases at larger n.
//
// Use simulated annealing to minimize rho(sigma).
// Force searches under structural constraints (e.g., "must contain a 4-block").
//
// Key question: do all minimization runs converge to the known extremal family?
//
// Compile: g++ -O3 -march=native -std=c++20 exotic_search.cpp -o exotic_search.exe

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <random>
#include <chrono>

using u32 = uint32_t;
static thread_local int dfs_best;

static void dfs(int cur, u32 vis, int nc, int dep,
                const u32* v, const u32* p, int n) {
    if (dep > dfs_best) dfs_best = dep;
    if (dfs_best >= n) return;
    u32 nbrs = (nc == 0 ? p[cur] : v[cur]) & ~vis;
    u32 full = (n < 32) ? ((1u << n) - 1) : 0xFFFFFFFFu;
    while (nbrs) {
        int nb = __builtin_ctz(nbrs); nbrs &= nbrs - 1;
        if (dep + 1 + __builtin_popcount(~(vis|(1u<<nb))&full) <= dfs_best) continue;
        dfs(nb, vis|(1u<<nb), 1-nc, dep+1, v, p, n);
        if (dfs_best >= n) return;
    }
}

static int rho_f(int n, const int* sigma) {
    if (n > 30) return -1; // bitmask limit
    int inv[32];
    for (int i = 0; i < n; i++) inv[sigma[i]] = i;
    u32 vn[32]={}, pn[32]={};
    for (int i = 0; i < n; i++) {
        int sv = sigma[i];
        if (sv > 0) vn[i] |= 1u<<inv[sv-1];
        if (sv < n-1) vn[i] |= 1u<<inv[sv+1];
        if (i > 0) pn[i] |= 1u<<(i-1);
        if (i < n-1) pn[i] |= 1u<<(i+1);
    }
    dfs_best = 1;
    for (int s = 0; s < n; s++) {
        u32 vis = 1u<<s;
        dfs(s, vis, 0, 1, vn, pn, n); if (dfs_best >= n) return dfs_best;
        dfs(s, vis, 1, 1, vn, pn, n); if (dfs_best >= n) return dfs_best;
    }
    return dfs_best;
}

static int count_blocks(int n, const int* sigma) {
    if (n <= 1) return n;
    int b = 1, cd = 0;
    for (int i = 1; i < n; i++) {
        int d = sigma[i] - sigma[i-1];
        if (d == 1 || d == -1) { if (cd == 0) cd = d; else if (d != cd) { b++; cd = d; } }
        else { b++; cd = 0; }
    }
    return b;
}

static void build_extremal(int m, int* sigma, int& n) {
    n = 4*m + 3;
    int pos = 0;
    sigma[pos]=pos+1; sigma[pos+1]=pos; pos+=2;
    sigma[pos]=pos+2; sigma[pos+1]=pos+1; sigma[pos+2]=pos; pos+=3;
    for (int t = 0; t < m-1; t++) {
        sigma[pos]=pos; pos+=1;
        sigma[pos]=pos+2; sigma[pos+1]=pos+1; sigma[pos+2]=pos; pos+=3;
    }
    sigma[pos]=pos+1; sigma[pos+1]=pos; pos+=2;
}

// Simulated annealing to minimize rho
static void sa_search(int n, int trials, std::mt19937& rng, const char* constraint_name,
                       bool (*constraint)(int, const int*)) {
    int best_sigma[32], best_rho = n, best_d = 0, best_b = 0;
    int target = (n + 3 + 1) / 2;

    for (int trial = 0; trial < trials; trial++) {
        int sigma[32];
        // Random starting permutation
        for (int i = 0; i < n; i++) sigma[i] = i;
        std::shuffle(sigma, sigma + n, rng);

        // If constraint, rejection-sample start
        if (constraint) {
            int attempts = 0;
            while (!constraint(n, sigma) && attempts++ < 1000)
                std::shuffle(sigma, sigma + n, rng);
            if (!constraint(n, sigma)) continue;
        }

        int cur_rho = rho_f(n, sigma);
        double temp = 2.0;

        for (int step = 0; step < 5000; step++) {
            // Random adjacent transposition
            int i = rng() % (n - 1);
            std::swap(sigma[i], sigma[i+1]);

            if (constraint && !constraint(n, sigma)) {
                std::swap(sigma[i], sigma[i+1]); // undo
                continue;
            }

            int new_rho = rho_f(n, sigma);
            int delta = new_rho - cur_rho;

            if (delta <= 0 || (temp > 0.01 && std::uniform_real_distribution<>(0,1)(rng) < exp(-delta/temp))) {
                cur_rho = new_rho;
            } else {
                std::swap(sigma[i], sigma[i+1]); // undo
            }
            temp *= 0.999;
        }

        if (cur_rho < best_rho || (cur_rho == best_rho && n - cur_rho > best_d)) {
            best_rho = cur_rho;
            best_d = n - cur_rho;
            best_b = count_blocks(n, sigma);
            memcpy(best_sigma, sigma, n * sizeof(int));
        }
    }

    printf("  [%s] n=%d: best_rho=%d d=%d b=%d target=%d  sigma=[",
           constraint_name, n, best_rho, best_d, best_b, target);
    for (int i = 0; i < n; i++) printf("%s%d", i?",":"", best_sigma[i]);
    printf("]  %s\n", best_rho >= target ? "OK" : "BELOW TARGET!");
}

// Constraints
static bool has_4block(int n, const int* s) {
    int run = 1, dir = 0;
    for (int i = 1; i < n; i++) {
        int d = s[i] - s[i-1];
        if ((d == 1 || d == -1) && (dir == 0 || d == dir)) {
            if (dir == 0) dir = d;
            run++;
            if (run >= 4) return true;
        } else { run = 1; dir = 0; if (d == 1 || d == -1) { dir = d; run = 2; } }
    }
    return false;
}

static bool all_blocks_le3(int n, const int* s) { return !has_4block(n, s); }

static bool has_ascending_3block(int n, const int* s) {
    int run = 1, dir = 0;
    for (int i = 1; i < n; i++) {
        int d = s[i] - s[i-1];
        if (d == 1 && (dir == 0 || dir == 1)) { dir = 1; run++; if (run >= 3) return true; }
        else { run = 1; dir = 0; if (d == 1) { dir = 1; run = 2; } else if (d == -1) { dir = -1; run = 2; } }
    }
    return false;
}

int main() {
    printf("=== Exotic Worst-Case Search ===\n\n");

    std::mt19937 rng(42);

    for (int n : {15, 19}) { // n=23,27 too slow for now
        printf("n=%d:\n", n);

        // Unconstrained search
        sa_search(n, 200, rng, "unconstrained", nullptr);

        // Must contain a 4-block
        sa_search(n, 200, rng, "has 4-block", has_4block);

        // All blocks ≤ 3
        sa_search(n, 200, rng, "blocks<=3", all_blocks_le3);

        // Must have ascending 3-block
        sa_search(n, 200, rng, "has asc-3", has_ascending_3block);

        // Compare with known extremal
        int ext[32]; int en;
        build_extremal((n-3)/4, ext, en);
        int er = rho_f(en, ext);
        printf("  [extremal family] n=%d: rho=%d d=%d b=%d\n\n", en, er, en-er, count_blocks(en, ext));
    }

    return 0;
}
