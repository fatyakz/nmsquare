#include <iostream>
#include <chrono>
#include <thread>
#include <mutex>
#include <ctime>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>

#pragma warning(disable : 4996)

#define LINE_READ_BUFFER_SIZE 512

struct S_best {
public:
    int matches;
    int n;
    int m;
    int e;
    int r;
};
struct S_thread { 
public:
    uint_fast32_t id;
    std::chrono::duration<double> time;
    std::time_t date;
    unsigned long long int cycles;
    int step; 
    int offset; 
    S_best best = {};
};
struct S_block {
public:
    unsigned long long int id;
    std::chrono::duration<double> time;
    std::time_t date;
    unsigned long long int cycles;
    S_best best = {};
    std::vector<S_thread> thread;
};
struct S_global {
public:
    std::vector<S_block> block;
    std::time_t date;
    std::chrono::duration<double> time;
    unsigned long long int cycles;
    S_best best = {};
    uint_fast32_t G_BLOCK_START;
    uint_fast32_t G_LIMIT;
    
    std::string G_BLOCK_FILE_PATH = "block.dat";
};

S_global g;
S_global data;

static std::mutex mlock;
static const uint_fast32_t G_NUM_THREADS = 14;
static uint_fast32_t G_B_DATA_CELLS = 6;

static uint_fast32_t R_MAX = 4294967295;
static unsigned long long int E_MAX = 18446744073709551615;

static uint_fast32_t T_CYCLES_DIVIDER = 1;
static std::string T_CYCLES_SYMBOL = "";
static uint_fast32_t B_CYCLES_DIVIDER = 1000;
static std::string B_CYCLES_SYMBOL = "K";
static uint_fast32_t G_CYCLES_DIVIDER = 1000;
static std::string G_CYCLES_SYMBOL = "M";

bool square(unsigned long long int x) { 
    if (x > 0) { 
        unsigned long long int sr = sqrt(x);
        return (sr * sr == x); 
    } return false; 
}

std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str) {
    std::vector<std::string> result;
    std::string line;
    std::getline(str, line);

    std::stringstream lineStream(line);
    std::string cell;

    while (std::getline(lineStream, cell, ',')) {
        result.push_back(cell);
    }
    if (!lineStream && cell.empty()) {
        result.push_back("");
    }
    return result;
}

bool B_WriteToFile(std::vector<S_block> block_vector, std::string path) {
    std::ofstream file_out{ path, std::ios_base::app };
    // appending with no new line
    file_out << "test";
    file_out << "test";
    return 1;
}

bool B_ReadFromFile(std::vector<S_block> &block_vector, std::string path) {
    std::ifstream file_out{ path, std::ifstream::binary };

    char line_buffer[LINE_READ_BUFFER_SIZE];
    int line_index = 0;

    S_block block_buffer;

    while (file_out.getline(line_buffer, 256)) {
        std::istringstream csv_buffer(line_buffer);
        std::vector<std::string> cells_buffer = getNextLineAndSplitIntoTokens(csv_buffer);

        block_buffer.best.n = stoi(cells_buffer[0]);
        block_buffer.best.m = stoi(cells_buffer[1]);
        block_buffer.best.e = stoi(cells_buffer[2]);
        block_buffer.best.r = stoi(cells_buffer[3]);
        block_buffer.best.matches = stoi(cells_buffer[4]);
        block_buffer.cycles = stoi(cells_buffer[5]);

        data.block.push_back(block_buffer);

        line_index++;
    }
    return 1;
}




static S_thread thr_SingleE(unsigned long long int t_E, uint_fast32_t t_offset, uint_fast32_t t_step) {
    unsigned long long int A, B, C, D, E, F, G, H, I;
    unsigned long long int t_cycles = 0, t_NM_limit = 0;
    
    uint_fast32_t t_matches = 0;
    
    S_thread t_thread{};
    t_thread.best.matches = 0;
    t_thread.id = t_E;
    t_thread.offset = t_offset;
    t_thread.step = t_step;

    E = t_E * t_E; 
    t_NM_limit = E - 1;
    
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    for (unsigned long long int lN = 1; lN < t_NM_limit; lN++) {
        for (unsigned long long int lM = t_offset + 1; lM < t_NM_limit; lM += t_step) {
            if (lN == lM) { 
                goto end; }
            if (lN + lM >= E) {
                break;
            }

            A = E + lN;
            B = E - lN - lM;
            C = E + lM;
            D = E - lN + lM;
            F = E + lN - lM;
            G = E - lM;
            H = E + lN + lM;
            I = E - lN;

            if (square(A) == true) { t_matches++; }
            if (square(B) == true) { t_matches++; }
            if (square(C) == true) { t_matches++; }
            if (square(D) == true) { t_matches++; }
            if (square(F) == true) { t_matches++; }
            if (square(G) == true) { t_matches++; }
            if (square(H) == true) { t_matches++; }
            if (square(I) == true) { t_matches++; }

        end:
            t_cycles++;
            if (t_matches >= t_thread.best.matches) {
                t_thread.best.e = E;
                t_thread.best.n = lN;
                t_thread.best.m = lM;
                t_thread.best.matches = t_matches + 1;
                t_thread.best.r = t_E; }

            if (t_thread.best.matches >= g.block[t_E].best.matches) {
                g.block[t_E].best = t_thread.best;
            }

            t_matches = 0;
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

    auto t_date = std::chrono::system_clock::now();
    t_thread.date = std::chrono::system_clock::to_time_t(t_date);
    t_cycles /= T_CYCLES_DIVIDER;
    t_thread.cycles = t_cycles;
    t_thread.time = t_time;

    mlock.lock();
        g.block[t_E].cycles += t_cycles / B_CYCLES_DIVIDER;
        g.block[t_E].thread[t_offset] = t_thread;

        std::cout << "THREAD [" << g.block[t_E].thread[t_offset].id << "+" << t_offset << "]" << 
            " cycles = " << g.block[t_E].thread[t_offset].cycles << T_CYCLES_SYMBOL <<
            " time="<< g.block[t_E].thread[t_offset].time.count() << " E=" << E <<
            " matches=" << g.block[t_E].thread[t_offset].best.matches << " n=" << g.block[t_E].thread[t_offset].best.n << 
            " m=" << g.block[t_E].thread[t_offset].best.m <<
            " e=" << g.block[t_E].thread[t_offset].best.e << " r=" << g.block[t_E].thread[t_offset].best.r << "\n";
    mlock.unlock();

    return t_thread;
}

int main()
{
start:
    //data.block.resize(0);
    B_ReadFromFile(g.block, g.G_BLOCK_FILE_PATH);

    std::cout << "Start: "; std::cin >> g.G_BLOCK_START;
    if (g.G_BLOCK_START == 0) { 
        return 0; }
    std::cout << "Blocks: "; std::cin >> g.G_LIMIT; std::cout << "\n";
    std::chrono::high_resolution_clock::time_point g1 = std::chrono::high_resolution_clock::now();

    auto g_date = std::chrono::system_clock::now();
    g.date = std::chrono::system_clock::to_time_t(g_date);

    g.block.resize(g.G_BLOCK_START + g.G_LIMIT);
    
    std::wcout << "BLOCKS " << g.G_BLOCK_START << "->" << g.G_BLOCK_START + g.G_LIMIT - 1 << 
        " date=" << std::ctime(&g.date) << "\n";

    for (int l = g.G_BLOCK_START; l < g.G_BLOCK_START + g.G_LIMIT; l++)
    {
        S_block g_block;
        g_block.id = l;

        g.block[g_block.id].thread.resize(G_NUM_THREADS);

        std::cout << "BLOCK [" << g_block.id << " of " << g.G_BLOCK_START + g.G_LIMIT  - 1 << 
            "] threads=" << G_NUM_THREADS <<"\n";

        auto b_date = std::chrono::system_clock::now(); 
        g.block[g_block.id].date = std::chrono::system_clock::to_time_t(b_date);

        std::chrono::high_resolution_clock::time_point b1 = std::chrono::high_resolution_clock::now();

        std::thread thr[G_NUM_THREADS];
        
        for (uint_fast32_t g_thread_offset = 0; g_thread_offset < G_NUM_THREADS; g_thread_offset++) {
            thr[g_thread_offset] = std::thread(thr_SingleE, g_block.id, g_thread_offset, G_NUM_THREADS);
        }
        for (uint_fast32_t g_thread_id = 0; g_thread_id < G_NUM_THREADS; g_thread_id++) {
            thr[g_thread_id].join();
        }

        std::chrono::high_resolution_clock::time_point b2 = std::chrono::high_resolution_clock::now();
        g.block[g_block.id].time = std::chrono::duration_cast<std::chrono::duration<double>>(b2 - b1);

        g.cycles += g.block[g_block.id].cycles / G_CYCLES_DIVIDER;

        if (g.block[g_block.id].best.matches >= g.best.matches) {
            g.best = g.block[g_block.id].best;
        }

        std::chrono::high_resolution_clock::time_point g2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> g_total_time = std::chrono::duration_cast<std::chrono::duration<double>>(g2 - g1);

        g.time = g_total_time;

        std::cout << "BLOCK [" << g_block.id << "] cycles=" << g.block[g_block.id].cycles << B_CYCLES_SYMBOL 
            << " time=" << g.block[g_block.id].time.count() << "s date=" << std::ctime(&g.block[g_block.id].date);
        std::cout << "BLOCK [" << g_block.id << "] best="<< g.block[g_block.id].best.matches 
            << " n=" << g.block[g_block.id].best.n << " m=" << g.block[g_block.id].best.m
            << " e=" << g.block[g_block.id].best.e << " r=" << g.block[g_block.id].best.r << "\n\n";
        std::cout << "G.CURRENT time=" << g.time.count() << "s cycles=" << g.cycles << G_CYCLES_SYMBOL <<
            " best=" << g.best.matches << " n=" << g.best.n << " m=" << g.best.m <<
            " e=" << g.best.e << " r=" << g.best.r << "\n\n";
    }
    std::chrono::high_resolution_clock::time_point g2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> g_total_time = std::chrono::duration_cast<std::chrono::duration<double>>(g2 - g1);

    g.time = g_total_time;

    std::cout << "END time=" << g.time.count() << "s cycles=" << g.cycles << G_CYCLES_SYMBOL << "\n";
    std::cout << "END best=" << g.best.matches << " n=" << g.best.n << " m=" << g.best.m << 
        " e=" << g.best.e << " r=" << g.best.r << "\n\n";

    goto start;

    return 0;
}
