// Agenda item 11: Deletion/contraction experiments on minimizers.
//
// For each sigma with d = max deficit at n:
//   - Delete each value v and standardize: does deficit decrease?
//   - Delete each block and standardize: how does rho change?
//   - Contract singletons: merge adjacent blocks
//
// Key question: do extremizers reduce to smaller extremizers?
//
// Compile: g++ -O3 -march=native -std=c++20 deletion_test.cpp -o deletion_test.exe

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <vector>

using u32 = uint32_t;
static thread_local int dfs_best;

static void dfs(int cur, u32 vis, int nc, int dep,
                const u32* v, const u32* p, int n) {
    if (dep > dfs_best) dfs_best = dep;
    if (dfs_best >= n) return;
    u32 nbrs = (nc == 0 ? p[cur] : v[cur]) & ~vis;
    u32 full = (1u << n) - 1;
    while (nbrs) {
        int nb = __builtin_ctz(nbrs); nbrs &= nbrs - 1;
        if (dep + 1 + __builtin_popcount(~(vis|(1u<<nb))&full) <= dfs_best) continue;
        dfs(nb, vis|(1u<<nb), 1-nc, dep+1, v, p, n);
        if (dfs_best >= n) return;
    }
}

static int rho_f(int n, const int* sigma) {
    int inv[40];
    for (int i = 0; i < n; i++) inv[sigma[i]] = i;
    u32 vn[40]={}, pn[40]={};
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

// Delete value v from sigma and standardize
static void delete_value(const int* sigma, int n, int del_val, int* out, int& out_n) {
    out_n = n - 1;
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (sigma[i] == del_val) continue;
        out[j++] = sigma[i] > del_val ? sigma[i] - 1 : sigma[i];
    }
}

// Delete position pos from sigma and standardize
static void delete_position(const int* sigma, int n, int del_pos, int* out, int& out_n) {
    out_n = n - 1;
    int del_val = sigma[del_pos];
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (i == del_pos) continue;
        out[j++] = sigma[i] > del_val ? sigma[i] - 1 : sigma[i];
    }
}

// Build extremal sigma
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

int main() {
    printf("=== Deletion/Contraction Experiments ===\n\n");

    // Test on extremal family
    for (int m = 1; m <= 6; m++) {
        int sigma[40], n;
        build_extremal(m, sigma, n);
        int r = rho_f(n, sigma);
        int d = n - r;
        int b = count_blocks(n, sigma);

        printf("--- Extremal m=%d, n=%d, rho=%d, d=%d, b=%d ---\n", m, n, r, d, b);

        // Delete each position
        printf("  Delete position i: (n-1, rho', d', b')\n");
        for (int i = 0; i < n; i++) {
            int out[40], on;
            delete_position(sigma, n, i, out, on);
            int r2 = rho_f(on, out);
            int d2 = on - r2;
            int b2 = count_blocks(on, out);
            if (d2 >= d - 1) // only show interesting cases
                printf("    pos %2d (val=%2d): n'=%d rho'=%d d'=%d b'=%d  delta_d=%+d  delta_b=%+d\n",
                       i, sigma[i], on, r2, d2, b2, d2-d, b2-b);
        }

        // Check: does deleting a singleton reduce to a smaller extremal family?
        if (m >= 2) {
            printf("  Singleton deletion test:\n");
            // Singletons are at positions 5, 9, 13, ...
            for (int t = 0; t < m-1; t++) {
                int sing_pos = 5 + 4*t;
                int out[40], on;
                delete_position(sigma, n, sing_pos, out, on);
                int r2 = rho_f(on, out);
                int d2 = on - r2;
                int b2 = count_blocks(on, out);
                printf("    delete singleton at pos %d: n'=%d d'=%d b'=%d\n",
                       sing_pos, on, d2, b2);
                // Is this the extremal family for m-1?
                int ext_sigma[40]; int ext_n;
                if (on == 4*(m-1)+3) {
                    build_extremal(m-1, ext_sigma, ext_n);
                    bool match = (ext_n == on);
                    for (int i = 0; match && i < on; i++)
                        if (out[i] != ext_sigma[i]) match = false;
                    printf("      matches extremal(m-1=%d)? %s\n", m-1, match ? "YES" : "no");
                }
            }
        }
        printf("\n");
    }

    // Also test some non-extremal worst cases
    printf("--- Non-extremal worst cases ---\n");
    // n=12, d=3: sigma=[0,1,2,5,4,3,8,7,6,9,10,11]
    {
        int sigma[] = {0,1,2,5,4,3,8,7,6,9,10,11};
        int n = 12;
        int r = rho_f(n, sigma);
        int d = n - r;
        int b = count_blocks(n, sigma);
        printf("sigma=[0,1,2,5,4,3,8,7,6,9,10,11] n=%d rho=%d d=%d b=%d\n", n, r, d, b);

        printf("  Interesting position deletions:\n");
        for (int i = 0; i < n; i++) {
            int out[40], on;
            delete_position(sigma, n, i, out, on);
            int r2 = rho_f(on, out);
            int d2 = on - r2;
            int b2 = count_blocks(on, out);
            printf("    pos %2d (val=%2d): d'=%d b'=%d  delta_d=%+d\n",
                   i, sigma[i], d2, b2, d2-d);
        }
    }

    return 0;
}
