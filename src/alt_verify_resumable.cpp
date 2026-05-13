// Resumable exhaustive ALT verifier.
//
// The natural checkpoint unit is a "chunk": one fixed (sigma[0], sigma[1])
// pair. For n, there are n*(n-1) chunks.  Each chunk is independent and
// enumerates all (n-2)! sigmas with that prefix.
//
// On startup:
//   If --resume is passed and <chunk_log> exists, read which chunks have
//   already completed and skip them.
//
// During the run:
//   When a chunk completes, ATOMICALLY append a single line to chunk_log:
//     CHUNK <id> <first> <second> <sigma_count> <loss_histogram> <worst_k> <worst_sigma>
//   All threads share the same chunk_log via a mutex.
//
// On completion (normal or interrupted):
//   Read chunk_log and aggregate. Print the final report.
//
// If the program crashes mid-run, just re-run with --resume and it picks up.
//
// Compile:
//   g++ -O3 -march=native -std=c++20 -fopenmp alt_verify_resumable.cpp \
//       -o alt_verify_resumable.exe
//
// Usage:
//   alt_verify_resumable.exe <n> [options]
//
// Options:
//   --log=<file>       Chunk log file (default: n<n>_chunks.log)
//   --resume           Skip chunks already in the log (default: overwrite)
//   --print-thr=<L>    Only echo sigmas with loss>=L to stdout (default: 5)
//   --max-print=<N>    Cap on echoed sigmas (default: 5000)
//   --aggregate-only   Do not run; just read log and print final report
//
// Example restart after a crash:
//   alt_verify_resumable.exe 15 --resume

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <atomic>
#include <mutex>
#include <chrono>
#include <ctime>
#include <thread>
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <omp.h>

using u64 = uint64_t;

// ---------- CLI config ----------
struct Config {
    int n = 0;
    std::string chunk_log = "";
    bool resume = false;
    int print_threshold = 5;
    u64 max_printed = 5000;
    bool aggregate_only = false;
};

static Config cfg;

// ---------- Global state ----------
static std::atomic<u64> total_sigma{0};
static std::atomic<u64> loss_hist[40];
static std::atomic<u64> printed{0};
static std::atomic<u64> chunks_completed{0};
static u64 total_chunks = 0;

// Thread-safe chunk log writing
static std::mutex log_mutex;
static std::mutex print_mutex;
static std::ofstream chunk_log_stream;

// Global worst sigma
static std::mutex worst_mutex;
static int best_global_worst_seen = 999;
static int best_global_sigma[32];

// ---------- Alternating walk ----------
struct AltWalker {
    int n;
    u64 v_nbrs[32];
    u64 p_nbrs[32];
    int best;

    void setup(const int* sig) {
        int inv[32];
        for (int i = 0; i < n; i++) inv[sig[i]] = i;
        for (int i = 0; i < n; i++) {
            v_nbrs[i] = 0;
            p_nbrs[i] = 0;
            int sv = sig[i];
            if (sv > 0)     v_nbrs[i] |= (u64)1 << inv[sv - 1];
            if (sv < n - 1) v_nbrs[i] |= (u64)1 << inv[sv + 1];
            if (i > 0)      p_nbrs[i] |= (u64)1 << (i - 1);
            if (i < n - 1)  p_nbrs[i] |= (u64)1 << (i + 1);
        }
    }

    void dfs(int cur, u64 visited, int next_type, int length) {
        if (length > best) best = length;
        if (best >= n) return;
        u64 nbrs = (next_type == 0 ? v_nbrs[cur] : p_nbrs[cur]) & ~visited;
        while (nbrs) {
            int nb = __builtin_ctzll(nbrs);
            nbrs &= nbrs - 1;
            dfs(nb, visited | ((u64)1 << nb), 1 - next_type, length + 1);
            if (best >= n) return;
        }
    }

    int longest_alt_walk() {
        best = 1;
        for (int start = 0; start < n; start++) {
            u64 vis = (u64)1 << start;
            dfs(start, vis, 0, 1);
            if (best >= n) return best;
            dfs(start, vis, 1, 1);
            if (best >= n) return best;
        }
        return best;
    }
};

// ---------- Symmetry ----------
static void apply_rev_P(const int* s, int* out, int n) {
    for (int i = 0; i < n; i++) out[i] = s[n - 1 - i];
}
static void apply_rev_Q(const int* s, int* out, int n) {
    for (int i = 0; i < n; i++) out[i] = n - 1 - s[i];
}
static void apply_swap(const int* s, int* out, int n) {
    for (int i = 0; i < n; i++) out[s[i]] = i;
}

static bool is_canonical(const int* sigma, int n) {
    int tmp[32], tmp2[32];
    for (int mask = 1; mask < 8; mask++) {
        std::memcpy(tmp, sigma, n * sizeof(int));
        if (mask & 1) { apply_rev_P(tmp, tmp2, n); std::memcpy(tmp, tmp2, n*sizeof(int)); }
        if (mask & 2) { apply_rev_Q(tmp, tmp2, n); std::memcpy(tmp, tmp2, n*sizeof(int)); }
        if (mask & 4) { apply_swap(tmp, tmp2, n); std::memcpy(tmp, tmp2, n*sizeof(int)); }
        for (int i = 0; i < n; i++) {
            if (tmp[i] < sigma[i]) return false;
            if (tmp[i] > sigma[i]) break;
        }
    }
    return true;
}

static int orbit_size(const int* sigma, int n) {
    int orbit[8][32];
    int tmp2[32];
    int count = 0;
    for (int mask = 0; mask < 8; mask++) {
        int tmp[32];
        std::memcpy(tmp, sigma, n * sizeof(int));
        if (mask & 1) { apply_rev_P(tmp, tmp2, n); std::memcpy(tmp, tmp2, n*sizeof(int)); }
        if (mask & 2) { apply_rev_Q(tmp, tmp2, n); std::memcpy(tmp, tmp2, n*sizeof(int)); }
        if (mask & 4) { apply_swap(tmp, tmp2, n); std::memcpy(tmp, tmp2, n*sizeof(int)); }
        bool seen = false;
        for (int k = 0; k < count; k++) {
            bool eq = true;
            for (int i = 0; i < n; i++) if (orbit[k][i] != tmp[i]) { eq = false; break; }
            if (eq) { seen = true; break; }
        }
        if (!seen) {
            std::memcpy(orbit[count], tmp, n * sizeof(int));
            count++;
        }
    }
    return count;
}

// ---------- Per-chunk state ----------
thread_local int sigma_buf[32];
thread_local AltWalker walker;
thread_local u64 chunk_sigma_count = 0;
thread_local u64 chunk_loss_hist[40];
thread_local int chunk_worst_seen;
thread_local int chunk_worst_sigma[32];

static void chunk_process_sigma() {
    if (!is_canonical(sigma_buf, cfg.n)) return;
    walker.n = cfg.n;
    walker.setup(sigma_buf);
    int k = walker.longest_alt_walk();
    int loss = cfg.n - k;
    int mult = orbit_size(sigma_buf, cfg.n);
    chunk_sigma_count += mult;
    chunk_loss_hist[loss] += mult;
    if (k < chunk_worst_seen) {
        chunk_worst_seen = k;
        for (int i = 0; i < cfg.n; i++) chunk_worst_sigma[i] = sigma_buf[i];
    }
    if (loss >= cfg.print_threshold) {
        u64 p = printed.fetch_add(1, std::memory_order_relaxed);
        if (p < cfg.max_printed) {
            std::lock_guard<std::mutex> lock(print_mutex);
            printf("  [orbit %d, loss=%d] sigma =", mult, loss);
            for (int i = 0; i < cfg.n; i++) printf(" %d", sigma_buf[i]);
            printf("\n");
            fflush(stdout);
        }
    }
}

static void enum_perms_rec(int depth, u64 used_mask) {
    if (depth == cfg.n) { chunk_process_sigma(); return; }
    for (int v = 0; v < cfg.n; v++) {
        if (used_mask & ((u64)1 << v)) continue;
        sigma_buf[depth] = v;
        enum_perms_rec(depth + 1, used_mask | ((u64)1 << v));
    }
}

// ---------- Chunk log I/O ----------
static void write_chunk_record(int chunk_id, int first, int second,
                                u64 count, const u64 hist[40],
                                int worst_k, const int* worst_sigma) {
    std::lock_guard<std::mutex> lock(log_mutex);
    chunk_log_stream << "CHUNK " << chunk_id << " " << first << " " << second
                     << " count=" << count << " hist=";
    bool first_h = true;
    for (int L = 0; L < 40; L++) {
        if (hist[L] > 0) {
            if (!first_h) chunk_log_stream << ",";
            chunk_log_stream << L << ":" << hist[L];
            first_h = false;
        }
    }
    chunk_log_stream << " worst_k=" << worst_k << " worst_sigma=";
    for (int i = 0; i < cfg.n; i++) {
        if (i > 0) chunk_log_stream << ",";
        chunk_log_stream << worst_sigma[i];
    }
    chunk_log_stream << "\n";
    chunk_log_stream.flush();
}

struct ChunkRecord {
    int chunk_id;
    int first, second;
    u64 count;
    u64 hist[40];
    int worst_k;
    int worst_sigma[32];
};

static std::vector<ChunkRecord> read_chunk_log(const std::string& path) {
    std::vector<ChunkRecord> out;
    std::ifstream f(path);
    if (!f.is_open()) return out;
    std::string line;
    while (std::getline(f, line)) {
        if (line.substr(0, 6) != "CHUNK ") continue;
        ChunkRecord r{};
        std::istringstream ss(line.substr(6));
        ss >> r.chunk_id >> r.first >> r.second;
        std::string tok;
        while (ss >> tok) {
            if (tok.substr(0, 6) == "count=") {
                r.count = std::stoull(tok.substr(6));
            } else if (tok.substr(0, 5) == "hist=") {
                std::string h = tok.substr(5);
                size_t pos = 0;
                while (pos < h.size()) {
                    size_t colon = h.find(':', pos);
                    size_t comma = h.find(',', colon);
                    if (colon == std::string::npos) break;
                    int L = std::stoi(h.substr(pos, colon - pos));
                    u64 v = std::stoull(h.substr(colon + 1,
                                        (comma == std::string::npos ? h.size() : comma) - colon - 1));
                    if (L < 40) r.hist[L] = v;
                    if (comma == std::string::npos) break;
                    pos = comma + 1;
                }
            } else if (tok.substr(0, 8) == "worst_k=") {
                r.worst_k = std::stoi(tok.substr(8));
            } else if (tok.substr(0, 12) == "worst_sigma=") {
                std::string w = tok.substr(12);
                size_t pos = 0;
                int idx = 0;
                while (pos < w.size() && idx < 32) {
                    size_t comma = w.find(',', pos);
                    r.worst_sigma[idx++] = std::stoi(
                        w.substr(pos, (comma == std::string::npos ? w.size() : comma) - pos));
                    if (comma == std::string::npos) break;
                    pos = comma + 1;
                }
            }
        }
        out.push_back(r);
    }
    return out;
}

// ---------- Progress reporter ----------
static std::atomic<bool> done_flag{false};

static std::string fmt_duration(double seconds) {
    int s = (int)seconds;
    int d = s / 86400;
    s %= 86400;
    int h = s / 3600;
    int m = (s % 3600) / 60;
    int sec = s % 60;
    char buf[64];
    if (d > 0)      snprintf(buf, sizeof(buf), "%dd %02dh %02dm", d, h, m);
    else if (h > 0) snprintf(buf, sizeof(buf), "%dh %02dm %02ds", h, m, sec);
    else if (m > 0) snprintf(buf, sizeof(buf), "%dm %02ds", m, sec);
    else            snprintf(buf, sizeof(buf), "%ds", sec);
    return std::string(buf);
}

static std::string fmt_number(u64 n) {
    char buf[64];
    double d = (double)n;
    if (d >= 1e12)      snprintf(buf, sizeof(buf), "%.2fT", d / 1e12);
    else if (d >= 1e9)  snprintf(buf, sizeof(buf), "%.2fB", d / 1e9);
    else if (d >= 1e6)  snprintf(buf, sizeof(buf), "%.2fM", d / 1e6);
    else if (d >= 1e3)  snprintf(buf, sizeof(buf), "%.2fK", d / 1e3);
    else                snprintf(buf, sizeof(buf), "%llu", (unsigned long long)n);
    return std::string(buf);
}

static void progress_reporter(std::chrono::steady_clock::time_point t0,
                               u64 pre_completed_chunks) {
    using namespace std::chrono;
    {
        std::lock_guard<std::mutex> lock(print_mutex);
        printf("\n[progress] updates every 5m (chunk-based progress)\n");
        printf("[progress]    elapsed  |  chunks    |  sigma seen  | curr k >= | ETA\n");
        fflush(stdout);
    }
    while (!done_flag.load()) {
        for (int s = 0; s < 300; s++) {
            std::this_thread::sleep_for(seconds(1));
            if (done_flag.load()) return;
        }

        auto now = steady_clock::now();
        double elapsed = duration<double>(now - t0).count();
        u64 done = chunks_completed.load();
        u64 total_done = done + pre_completed_chunks;
        u64 seen = total_sigma.load();
        double frac = total_chunks > 0 ? (double)total_done / (double)total_chunks : 0.0;

        // ETA based on THIS session's throughput (exclude pre_completed)
        double eta = -1;
        if (done > 0) {
            double per_chunk = elapsed / (double)done;
            eta = per_chunk * (double)(total_chunks - total_done);
        }

        int max_loss_seen = 0;
        for (int L = 39; L >= 0; L--) {
            if (loss_hist[L].load() > 0) { max_loss_seen = L; break; }
        }
        int curr_min_k = cfg.n - max_loss_seen;

        std::lock_guard<std::mutex> lock(print_mutex);
        printf("[progress] %10s | %4llu/%-4llu | %12s |    %2d     | %s\n",
               fmt_duration(elapsed).c_str(),
               (unsigned long long)total_done, (unsigned long long)total_chunks,
               fmt_number(seen).c_str(),
               curr_min_k,
               eta >= 0 ? fmt_duration(eta).c_str() : "?");
        fflush(stdout);
    }
}

// ---------- Aggregate ----------
static void print_final_report(const std::vector<ChunkRecord>& all_chunks, double dt) {
    u64 grand_total_sigma = 0;
    u64 grand_hist[40] = {0};
    int grand_worst_k = cfg.n + 1;
    int grand_worst_sigma[32] = {0};
    for (const auto& r : all_chunks) {
        grand_total_sigma += r.count;
        for (int L = 0; L < 40; L++) grand_hist[L] += r.hist[L];
        if (r.worst_k < grand_worst_k) {
            grand_worst_k = r.worst_k;
            for (int i = 0; i < cfg.n; i++) grand_worst_sigma[i] = r.worst_sigma[i];
        }
    }

    double n_factorial = 1.0;
    for (int i = 2; i <= cfg.n; i++) n_factorial *= i;

    printf("\n");
    printf("=====================================================================\n");
    printf(" FINAL REPORT (aggregated from %zu chunks)\n", all_chunks.size());
    printf("=====================================================================\n");
    time_t end_time = time(nullptr);
    printf(" finish time           = %s", ctime(&end_time));
    printf(" elapsed (session)     = %s (%.1f s)\n", fmt_duration(dt).c_str(), dt);
    printf(" total sigmas counted  = %llu  (expected %g)\n",
           (unsigned long long)grand_total_sigma, n_factorial);

    if (all_chunks.empty()) {
        printf("\n NO CHUNKS IN LOG YET.\n");
        printf(" The verifier may still be running, but no chunk has completed.\n");
        printf(" Each chunk takes ~%d min for n=%d (estimate).\n",
               cfg.n == 15 ? 60 : cfg.n == 14 ? 5 : 1, cfg.n);
        printf(" Check %s again later, or tail the stdout file.\n", cfg.chunk_log.c_str());
        printf("=====================================================================\n");
        return;
    }

    printf("\n");
    printf(" LOSS HISTOGRAM (counts include symmetry multiplicity)\n");
    printf(" -----------------------------------------------------\n");
    u64 grand_total = 0;
    for (int L = 0; L < 40; L++) {
        u64 v = grand_hist[L];
        if (v > 0) {
            double pct = 100.0 * (double)v / n_factorial;
            printf("   loss=%2d  (k=%2d):  %15llu  (%7.4f%%)\n",
                   L, cfg.n - L, (unsigned long long)v, pct);
            grand_total += v;
        }
    }
    printf("   ----------------------------------------------\n");
    printf("   total              %15llu\n", (unsigned long long)grand_total);

    int worst_loss = cfg.n - grand_worst_k;
    double ratio = (double)grand_worst_k / (double)cfg.n;

    printf("\n");
    printf(" HEADLINES\n");
    printf(" ---------\n");
    printf("   min rho  (longest PC path length) = %d\n", grand_worst_k);
    printf("   max loss (= n - min rho)          = %d\n", worst_loss);
    printf("   rho / n                           = %.4f\n", ratio);
    printf("   ALT-linear conjecture:  rho >= ceil((n+3)/2) = %d\n",
           (cfg.n + 3 + 1) / 2);
    if (cfg.n % 4 == 3) {
        printf("   extremal family (n=4m+3) predicts tight loss = %d\n",
               (cfg.n - 3) / 2);
        printf("   match with extremal family: %s\n",
               worst_loss == (cfg.n - 3) / 2 ? "YES (expected)" :
               worst_loss < (cfg.n - 3) / 2 ? "worst < extremal ??" :
                                               "WORSE than extremal (!!)");
    }
    printf("   ALT-linear SATISFIED for this n: %s\n",
           grand_worst_k >= (cfg.n + 3 + 1) / 2 ? "YES" : "NO (!!)");
    printf("   a worst sigma:");
    for (int i = 0; i < cfg.n; i++) printf(" %d", grand_worst_sigma[i]);
    printf("\n");
    printf("=====================================================================\n");
    fflush(stdout);
}

// ---------- Main ----------
int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n> [options]\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  --log=<file>       chunk log file (default: n<n>_chunks.log)\n");
        fprintf(stderr, "  --resume           skip chunks already in log\n");
        fprintf(stderr, "  --print-thr=<L>    print loss >= L (default 5)\n");
        fprintf(stderr, "  --max-print=<N>    cap echoed sigmas (default 5000)\n");
        fprintf(stderr, "  --aggregate-only   just read log & print report\n");
        return 1;
    }
    cfg.n = atoi(argv[1]);
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--resume") cfg.resume = true;
        else if (a == "--aggregate-only") cfg.aggregate_only = true;
        else if (a.substr(0, 6) == "--log=") cfg.chunk_log = a.substr(6);
        else if (a.substr(0, 12) == "--print-thr=") cfg.print_threshold = atoi(a.c_str() + 12);
        else if (a.substr(0, 12) == "--max-print=") cfg.max_printed = strtoull(a.c_str() + 12, 0, 10);
    }
    if (cfg.chunk_log.empty()) {
        char buf[64];
        snprintf(buf, sizeof(buf), "n%d_chunks.log", cfg.n);
        cfg.chunk_log = buf;
    }
    for (int L = 0; L < 40; L++) loss_hist[L] = 0;

    int nthreads = omp_get_max_threads();
    total_chunks = (u64)cfg.n * (u64)(cfg.n - 1);
    double n_factorial = 1.0;
    for (int i = 2; i <= cfg.n; i++) n_factorial *= i;

    // --- Aggregate-only mode ---
    if (cfg.aggregate_only) {
        printf("Reading %s ...\n", cfg.chunk_log.c_str());
        auto records = read_chunk_log(cfg.chunk_log);
        printf("Found %zu chunk records.\n", records.size());
        print_final_report(records, 0);
        return 0;
    }

    // --- Determine which chunks are already done (if --resume) ---
    std::set<int> done_chunks;
    std::vector<ChunkRecord> existing_records;
    if (cfg.resume) {
        existing_records = read_chunk_log(cfg.chunk_log);
        for (const auto& r : existing_records) done_chunks.insert(r.chunk_id);
        // Seed global totals from existing records
        for (const auto& r : existing_records) {
            total_sigma.fetch_add(r.count, std::memory_order_relaxed);
            for (int L = 0; L < 40; L++) loss_hist[L].fetch_add(r.hist[L]);
            if (r.worst_k < best_global_worst_seen) {
                best_global_worst_seen = r.worst_k;
                for (int i = 0; i < cfg.n; i++) best_global_sigma[i] = r.worst_sigma[i];
            }
        }
    }

    time_t now_time = time(nullptr);
    printf("=====================================================================\n");
    printf(" ALT RESUMABLE VERIFIER\n");
    printf("=====================================================================\n");
    printf(" n                  = %d\n", cfg.n);
    printf(" threads            = %d\n", nthreads);
    printf(" total perms        = %s\n", fmt_number((u64)n_factorial).c_str());
    printf(" total chunks       = %llu\n", (unsigned long long)total_chunks);
    printf(" chunks in log      = %zu  (will skip these if --resume)\n",
           done_chunks.size());
    printf(" chunks to run      = %llu\n", (unsigned long long)(total_chunks - done_chunks.size()));
    printf(" chunk log          = %s\n", cfg.chunk_log.c_str());
    printf(" print threshold    = loss >= %d  (first %llu echoed)\n",
           cfg.print_threshold, (unsigned long long)cfg.max_printed);
    printf(" start time         = %s", ctime(&now_time));
    printf("=====================================================================\n");
    fflush(stdout);

    // Open chunk log for appending
    chunk_log_stream.open(cfg.chunk_log, std::ios::app);
    if (!chunk_log_stream.is_open()) {
        fprintf(stderr, "ERROR: cannot open chunk log: %s\n", cfg.chunk_log.c_str());
        return 1;
    }

    auto t0 = std::chrono::steady_clock::now();
    std::thread reporter_thread(progress_reporter, t0, (u64)done_chunks.size());

    printf("\n--- sigmas with loss >= %d (first %llu echoed) ---\n",
           cfg.print_threshold, (unsigned long long)cfg.max_printed);
    fflush(stdout);

    // Build the list of chunks to run
    std::vector<int> chunks_to_run;
    for (int pair = 0; pair < (int)total_chunks; pair++) {
        if (done_chunks.count(pair) == 0) chunks_to_run.push_back(pair);
    }

    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (int idx = 0; idx < (int)chunks_to_run.size(); idx++) {
            int pair = chunks_to_run[idx];
            int first = pair / (cfg.n - 1);
            int second_raw = pair % (cfg.n - 1);
            int second = (second_raw >= first) ? (second_raw + 1) : second_raw;

            // Reset per-chunk state
            chunk_sigma_count = 0;
            for (int L = 0; L < 40; L++) chunk_loss_hist[L] = 0;
            chunk_worst_seen = cfg.n + 1;
            for (int i = 0; i < cfg.n; i++) chunk_worst_sigma[i] = 0;

            sigma_buf[0] = first;
            sigma_buf[1] = second;
            u64 used = ((u64)1 << first) | ((u64)1 << second);
            enum_perms_rec(2, used);

            // Write chunk record to log
            write_chunk_record(pair, first, second, chunk_sigma_count,
                               chunk_loss_hist, chunk_worst_seen, chunk_worst_sigma);

            // Update globals
            total_sigma.fetch_add(chunk_sigma_count, std::memory_order_relaxed);
            for (int L = 0; L < 40; L++) loss_hist[L].fetch_add(chunk_loss_hist[L]);
            {
                std::lock_guard<std::mutex> lock(worst_mutex);
                if (chunk_worst_seen < best_global_worst_seen) {
                    best_global_worst_seen = chunk_worst_seen;
                    for (int i = 0; i < cfg.n; i++) best_global_sigma[i] = chunk_worst_sigma[i];
                }
            }
            chunks_completed.fetch_add(1, std::memory_order_relaxed);
        }
    }

    done_flag = true;
    reporter_thread.join();
    chunk_log_stream.close();

    auto t1 = std::chrono::steady_clock::now();
    double dt = std::chrono::duration<double>(t1 - t0).count();

    // Re-read log for consistent aggregation
    auto all_records = read_chunk_log(cfg.chunk_log);
    print_final_report(all_records, dt);
    return 0;
}
