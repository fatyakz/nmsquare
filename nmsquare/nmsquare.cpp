#include <iostream>
#include <chrono>
#include <thread>
#include <mutex>
#include <future>
#include <vector>

using namespace std::chrono;
std::mutex mlock;

struct S_matches {
    int matches;
    int n;
    int m;
    int e;
    int root;
};

static S_matches g_best;

bool square(unsigned long long int x)
{
    if (x > 0) { long long sr = sqrt(x);
        return (sr * sr == x); }
    return false;
}

static S_matches thr_SquareGen(unsigned long long int t_start, unsigned long long int t_limit, unsigned long long int t_step) {
    unsigned long long int A, B, C, D, E, F, G, H, I;
    unsigned long long int t_cycles = 0;
    int t_matches = 0;
    S_matches t_best;
    t_best.matches = 0;
    for (unsigned long long int i = t_start; i < t_limit; i += t_step) {
        E = i * i;
        unsigned long long int N_limit = E - 1;
        for (unsigned long long int lN = 1; lN < N_limit; lN++) {
            for (unsigned long long int lM = 1; lM < N_limit; lM++) {
                if (lN == lM) { goto flash; }
                A = E + lN;
                if (square(A) != true) { goto flash; } t_matches++;
                B = E - lN - lM;
                if (square(B) != true) { goto flash; } t_matches++;
                C = E + lM;
                if (square(C) != true) { goto flash; } t_matches++;
                D = E - lN + lM;
                if (square(D) != true) { goto flash; } t_matches++;
                F = E + lN - lM;
                if (square(F) != true) { goto flash; } t_matches++;
                G = E - lM;
                if (square(G) != true) { goto flash; } t_matches++;
                H = E + lN + lM;
                if (square(H) != true) { goto flash; } t_matches++;
                I = E - lN;
                if (square(I) != true) { goto flash; } t_matches++;
                std::cout << "GOTTEM E:" << E << " n:" << lN << " m:" << lM << "\n";
            flash:
                t_cycles++;
                if (t_matches > t_best.matches) { t_best.e = E; t_best.n = lN; t_best.m = lM; t_best.matches = t_matches + 1; t_best.root = i; }
                t_matches = 0;
            }
        }
    }
    mlock.lock();
    if (t_best.matches > g_best.matches) {
        g_best.e = E; g_best.n = t_best.n; g_best.m = t_best.m; g_best.matches = t_best.matches; g_best.root = t_best.root;
    }
    std::cout << "thread " << t_start << ": " << t_cycles << " cycles, E="<< E << 
        " t_best=" << t_best.matches << " e=" << t_best.e << " n=" << t_best.n << " m=" << t_best.m << " r=" << t_best.root << "\n";
    mlock.unlock();
    return t_best;
}

int main()
{
    static const unsigned long long int g_numthreads = 14;
    static unsigned long long int g_start = 3;
    static unsigned long long int g_limit = 16;
    static unsigned long long int g_max = 18446744073709551615;

top:
    g_best.matches = 0;
    std::cout << "Start: "; std::cin >> g_start;
    if (g_start == 0) { return 0; }
    std::cout << "Blocks: "; std::cin >> g_limit; std::cout << "\n\n";

    high_resolution_clock::time_point g1 = high_resolution_clock::now();
    for (int l = g_start; l < g_start + (g_limit * g_numthreads); l += g_numthreads)
    {
        std::cout << "Block: " << l << " to " << l + g_numthreads - 1 << " (lim=" << g_start + (g_limit * g_numthreads) - 1 << ")\n\n";
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        std::thread thr[g_numthreads];
        for (unsigned long long int i = 0; i < g_numthreads; i++) { thr[i] = std::thread(thr_SquareGen, l + i, l + g_numthreads, g_numthreads); }
        for (unsigned long long int i = 0; i < g_numthreads; i++) { thr[i].join(); }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> total_time = duration_cast<duration<double>>(t2 - t1);
        std::cout << "Block time : " << total_time.count() << " s\n\n";
        std::cout << "CURRENT g_best=" << g_best.matches << " e=" << g_best.e << " n=" << g_best.n << " m=" << g_best.m << " r=" << g_best.root <<"\n\n";
    }
    high_resolution_clock::time_point g2 = high_resolution_clock::now();
    duration<double> g_total_time = duration_cast<duration<double>>(g2 - g1);
    std::cout << "Total time : " << g_total_time.count() << " s\n";
    std::cout << "FINAL g_best=" << g_best.matches << " e=" << g_best.e << " n=" << g_best.n << " m=" << g_best.m << " r=" << g_best.root << "\n\n";

    goto top;
    return 0;
}
