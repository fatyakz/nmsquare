#define LINE_READ_BUFFER_SIZE 512

#pragma warning(disable : 4996)

#include <cmath>
#include <iostream>
#include <chrono>
#include <thread>
#include <mutex>
#include <ctime>
#include <time.h>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>

struct S_best {
public:
    int         matches;
    int         n;
    int         m;
    int         e;
    int         r;
};
struct S_thread {
public:
    uint_fast32_t                   id;
    std::chrono::duration<double>   time;
    std::time_t                     date;
    unsigned long long int          cycles;
    int                             step;
    int                             offset;
    S_best                          best = {};
};
struct S_block {
public:
    unsigned long long int          id;
    std::chrono::duration<double>   time;
    std::time_t                     date;
    unsigned long long int          cycles;
    S_best                          best = {};
    std::vector<S_thread>           thread;
    uint_fast32_t                   state; //0=incomplete 1=running 2=finished
};
struct S_global {
public:
    std::vector<S_block>            block;
    std::time_t                     date;
    std::chrono::duration<double>   time;
    unsigned long long int          cycles;
    S_best                          best = {};
    uint_fast32_t                   G_BLOCK_START;
    uint_fast32_t                   G_LIMIT;
    std::string                     G_BLOCK_FILE_PATH = "block.dat";
    uint_fast32_t                   G_NUM_THREADS;
};
struct S_cell_index {
    uint_fast32_t _cell_count = 10;

    uint_fast32_t b_state = 0;
    uint_fast32_t b_id = 1;
    uint_fast32_t b_time = 2;
    uint_fast32_t b_date = 3;
    uint_fast32_t b_cycles = 4;
    uint_fast32_t b_n = 5;
    uint_fast32_t b_m = 6;
    uint_fast32_t b_e = 7;
    uint_fast32_t b_r = 8;
    uint_fast32_t b_matches = 9;
};

S_cell_index cell_index;
S_global g;
S_global data;

static std::mutex mlock;

static uint_fast32_t G_NUM_THREADS      = 14;
static uint_fast32_t G_B_DATA_CELLS     = 6;
static uint_fast32_t R_MAX              = 4294967295;
static unsigned long long int E_MAX     = 18446744073709551615;
std::string G_BLOCK_FILE_PATH_DEFAULT   = "block.dat";
static uint_fast32_t T_CYCLES_DIVIDER   = 1;
static std::string T_CYCLES_SYMBOL      = "";
static uint_fast32_t B_CYCLES_DIVIDER   = 1000;
static std::string B_CYCLES_SYMBOL      = "K";
static uint_fast32_t G_CYCLES_DIVIDER   = 1000;
static std::string G_CYCLES_SYMBOL      = "M";

bool square(unsigned long long int x) {
    if (x > 0) {
        unsigned long long int sr = sqrt(x);
        return (sr * sr == x);
    } return false;
}

inline bool file_exists(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

std::vector<std::string> B_ReadLine(std::istream& str) {
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
    std::ofstream file_out{ path, std::ios::trunc };

    if (block_vector.size() > 0) {
        std::vector<std::string> line_buffer;
        
        for (uint_fast32_t block_index = 0; block_index < block_vector.size(); block_index++) {
            line_buffer.resize(cell_index._cell_count);

            line_buffer[cell_index.b_cycles] =      std::to_string(block_vector[block_index].cycles);
            line_buffer[cell_index.b_date] =        std::to_string(block_vector[block_index].date);
            line_buffer[cell_index.b_n] =           std::to_string(block_vector[block_index].best.n);
            line_buffer[cell_index.b_m] =           std::to_string(block_vector[block_index].best.m);
            line_buffer[cell_index.b_e] =           std::to_string(block_vector[block_index].best.e);
            line_buffer[cell_index.b_r] =           std::to_string(block_vector[block_index].best.r);
            line_buffer[cell_index.b_matches] =     std::to_string(block_vector[block_index].best.matches);
            line_buffer[cell_index.b_id] =          std::to_string(block_vector[block_index].id);
            line_buffer[cell_index.b_state] =       std::to_string(block_vector[block_index].state);
            line_buffer[cell_index.b_time] =        std::to_string(block_vector[block_index].time.count());
            
            file_out << line_buffer[0];

            for (uint_fast32_t cell_pointer = 1; cell_pointer < cell_index._cell_count; cell_pointer++) {
                file_out << ",";
                file_out << line_buffer[cell_pointer];
            }

            if (block_index < block_vector.size() - 1) { file_out << "\n"; }

            line_buffer.resize(0);
        }
    }
    return 1;
}

template <typename T>
auto seconds_to_duration(T seconds) {
    return std::chrono::duration<T, std::ratio<1>>(seconds);
}

bool B_ReadFromFile(std::vector<S_block>& block_vector, std::string path) {
    std::ifstream file_out{ path, std::ifstream::binary };

    char line_buffer[LINE_READ_BUFFER_SIZE];
    int line_index = 0;

    cell_index.b_cycles;

    S_block block_buffer;

    std::string s; int file_lines = 0;
    std::ifstream in; in.open(path);
    while (!in.eof()) { getline(in, s); file_lines++; }

    data.block.resize(file_lines);

    while (file_out.getline(line_buffer, LINE_READ_BUFFER_SIZE)) {
        std::istringstream csv_buffer(line_buffer);
        std::vector<std::string> cells_buffer = B_ReadLine(csv_buffer);

        if (stoi(cells_buffer[cell_index.b_state]) != 0) {
            block_buffer.best.n =       stoi(cells_buffer[cell_index.b_n]);
            block_buffer.best.m =       stoi(cells_buffer[cell_index.b_m]);
            block_buffer.best.e =       stoi(cells_buffer[cell_index.b_e]);
            block_buffer.best.r =       stoi(cells_buffer[cell_index.b_r]);
            block_buffer.best.matches = stoi(cells_buffer[cell_index.b_matches]);
            block_buffer.cycles =       stol(cells_buffer[cell_index.b_cycles]);
            block_buffer.date =         stoi(cells_buffer[cell_index.b_date]);
            block_buffer.id =           stoi(cells_buffer[cell_index.b_id]);
            block_buffer.state =        stoi(cells_buffer[cell_index.b_state]);
            block_buffer.time =         seconds_to_duration(stod(cells_buffer[cell_index.b_time]));

            data.block[line_index] = block_buffer;

            if (block_buffer.best.matches >= g.best.matches) { g.best = block_buffer.best; }
        }
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
                goto end;
            }
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
                t_thread.best.matches = t_matches + 1;
                t_thread.best.e = E;
                t_thread.best.n = lN;
                t_thread.best.m = lM;
                t_thread.best.r = t_E;
            }

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
    g.block[t_E].cycles += t_thread.cycles / B_CYCLES_DIVIDER;
    g.block[t_E].thread[t_offset] = t_thread;
    g.block[t_E].id = t_E;
    g.block[t_E].state = 2;

    std::cout << "THR   [" << g.block[t_E].thread[t_offset].id << "+" << t_offset << "]" <<
        " cycles = " << g.block[t_E].thread[t_offset].cycles << T_CYCLES_SYMBOL <<
        " time=" << g.block[t_E].thread[t_offset].time.count() << 
        " E=" << E <<
        " matches=" << g.block[t_E].thread[t_offset].best.matches << 
        " n=" << g.block[t_E].thread[t_offset].best.n <<
        " m=" << g.block[t_E].thread[t_offset].best.m <<
        " e=" << g.block[t_E].thread[t_offset].best.e << 
        " r=" << g.block[t_E].thread[t_offset].best.r << "\n";
    mlock.unlock();

    return t_thread;
}

int main()
{
reset:
    g.time = std::chrono::milliseconds::zero();
    g.cycles = 0;

    std::cout << "[nmSquare]\nINIT  Threads:";
    std::cin >> g.G_NUM_THREADS;
    if (g.G_NUM_THREADS == 0) {
        return 0;
    }
    std::cout << "INIT  (" << g.G_BLOCK_FILE_PATH << "):";
    std::cin >> g.G_BLOCK_FILE_PATH;

    if (!file_exists(g.G_BLOCK_FILE_PATH)) {
        std::cout << "INIT  " << g.G_BLOCK_FILE_PATH << " does not exist, reverting to default: " << G_BLOCK_FILE_PATH_DEFAULT <<"\n";
        g.G_BLOCK_FILE_PATH = G_BLOCK_FILE_PATH_DEFAULT;
    }

    std::cout << "DATA  read=" << g.G_BLOCK_FILE_PATH << "\n";

    B_ReadFromFile(g.block, g.G_BLOCK_FILE_PATH);

    if (data.block.size() > g.block.size()) {
        std::wcout << "DATA  update available, importing " << data.block.size() << " blocks...\n";
        g.block.resize(data.block.size());

        for (uint_fast32_t block_pointer = 1; block_pointer < data.block.size(); block_pointer++) {
            if (g.block[block_pointer].state == 0) {
                g.block[block_pointer] = data.block[block_pointer];
            }
        }

        long block_count = 0;
        for (int i = 0; i < data.block.size(); i++) {
            if (data.block[i].state != 0) { block_count++; }
        }
        std::wcout << "DATA  up to date, " << data.block.size() - 1 << " blocks in file, " << block_count << " completed\n";
    }
    else {
        long block_count = 0;
        for (int i = 0; i < data.block.size(); i++) {
            if (data.block[i].state != 0) { block_count++; }
        }
        std::wcout << "DATA  up to date, " << data.block.size() - 1 << " blocks in file, " << block_count << " completed\n";
    }

start:
    std::cout << "[NMS] g_time=" << g.time.count() <<
        "s g_cycles=" << g.cycles << G_CYCLES_SYMBOL <<
        " best=" << g.best.matches <<
        " n=" << g.best.n <<
        " m=" << g.best.m <<
        " e=" << g.best.e <<
        " r=" << g.best.r << "\n";

    std::cout << "INIT  Start: "; std::cin >> g.G_BLOCK_START;

    if (g.G_BLOCK_START == 0) {
        goto reset;
    }

    std::cout << "INIT  Blocks: "; std::cin >> g.G_LIMIT; std::cout;
    std::chrono::high_resolution_clock::time_point g1 = std::chrono::high_resolution_clock::now();

    auto g_date = std::chrono::system_clock::now();
    g.date = std::chrono::system_clock::to_time_t(g_date);

    if (g.block.size() < (g.G_BLOCK_START + g.G_LIMIT)) { g.block.resize(g.G_BLOCK_START + g.G_LIMIT); }

    std::wcout << "BLOCK " << g.G_BLOCK_START << "->" << g.G_BLOCK_START + g.G_LIMIT - 1 <<
        " date=" << std::ctime(&g.date);

    for (int l = g.G_BLOCK_START; l < g.G_BLOCK_START + g.G_LIMIT; l++) {
        
        
        std::cout << "DATA  read=" << g.G_BLOCK_FILE_PATH << "\n";

        B_ReadFromFile(g.block, g.G_BLOCK_FILE_PATH);

        if (data.block.size() > g.block.size()) {
            std::wcout << "DATA  update available, importing " << data.block.size() << " blocks...\n";
            g.block.resize(data.block.size());

            for (uint_fast32_t block_pointer = 1; block_pointer < data.block.size(); block_pointer++) {
                if (g.block[block_pointer].state == 0) {
                    g.block[block_pointer] = data.block[block_pointer];
                }
            }

            long block_count = 0;
            for (int i = 0; i < data.block.size(); i++) {
                if (data.block[i].state != 0) { block_count++; }
            }
            std::wcout << "DATA  up to date, " << data.block.size() - 1 << " blocks in file, " << block_count << " completed\n";
        }
        else {
            long block_count = 0;
            for (int i = 0; i < data.block.size(); i++) {
                if (data.block[i].state != 0) { block_count++; }
            }
            std::wcout << "DATA  up to date, " << data.block.size() - 1 << " blocks in file, " << block_count << " completed\n";
        }

        if (g.block[l].state == 0) {
            S_block g_block;
            g_block.id = l;
            std::cout << "BLOCK [" << g_block.id << " of " << g.G_BLOCK_START + g.G_LIMIT - 1 <<
                "] threads=" << g.G_NUM_THREADS << " running...\n";

            g.block[g_block.id].thread.resize(g.G_NUM_THREADS);
            auto b_date = std::chrono::system_clock::now();
            g.block[g_block.id].date = std::chrono::system_clock::to_time_t(b_date);

            std::chrono::high_resolution_clock::time_point b1 = std::chrono::high_resolution_clock::now();

            std::vector<std::thread> thr(g.G_NUM_THREADS);

            for (uint_fast32_t g_thread_offset = 0; g_thread_offset < g.G_NUM_THREADS; g_thread_offset++) {
                thr[g_thread_offset] = std::thread(thr_SingleE, g_block.id, g_thread_offset, g.G_NUM_THREADS);
            }
            for (uint_fast32_t g_thread_id = 0; g_thread_id < g.G_NUM_THREADS; g_thread_id++) {
                thr[g_thread_id].join();
            }

            std::chrono::high_resolution_clock::time_point b2 = std::chrono::high_resolution_clock::now();
            g.block[g_block.id].time = std::chrono::duration_cast<std::chrono::duration<double>>(b2 - b1);

            g.cycles += g.block[g_block.id].cycles / G_CYCLES_DIVIDER;

            if (g.block[g_block.id].best.matches >= g.best.matches) {
                g.best = g.block[g_block.id].best;
            }

            std::cout << "BLOCK [" << g_block.id << "] cycles=" << g.block[g_block.id].cycles << B_CYCLES_SYMBOL <<
                " time=" << g.block[g_block.id].time.count() << 
                "s date=" << std::ctime(&g.block[g_block.id].date);

            std::cout << "BLOCK [" << g_block.id << "] best=" << g.block[g_block.id].best.matches <<
                " n=" << g.block[g_block.id].best.n << 
                " m=" << g.block[g_block.id].best.m << 
                " e=" << g.block[g_block.id].best.e << 
                " r=" << g.block[g_block.id].best.r << "\n";

            std::chrono::high_resolution_clock::time_point g2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> g_total_time = std::chrono::duration_cast<std::chrono::duration<double>>(g2 - g1);

            g.time += g_total_time;

            std::cout << "END   time=" << g.time.count() << "s cycles=" << g.cycles << G_CYCLES_SYMBOL << "\n";
            std::cout << "END   best=" << g.best.matches << " n=" << g.best.n << " m=" << g.best.m <<
                " e=" << g.best.e << " r=" << g.best.r << "\n";

            std::cout << "BEST  time=" << g.time.count() << 
                "s cycles=" << g.cycles << G_CYCLES_SYMBOL <<
                " best=" << g.best.matches << 
                " n=" << g.best.n << 
                " m=" << g.best.m <<
                " e=" << g.best.e << 
                " r=" << g.best.r << "\n";

            bool update = 0;

            if (g.block.size() > data.block.size()) {
                data.block.resize(g.block.size());
                update = 1;
            }
            else {
                g.block.resize(data.block.size());
            }

            for (uint_fast32_t block_pointer = 1; block_pointer < g.block.size(); block_pointer++) {
                if (data.block[block_pointer].state == 0 && g.block[block_pointer].state != 0) {
                    data.block[block_pointer] = g.block[block_pointer];
                    update = 1;
                }
        }
        if (update) { 
            std::cout << "DATA  write=" << g.G_BLOCK_FILE_PATH << "\n";
            B_WriteToFile(data.block, g.G_BLOCK_FILE_PATH); 
            std::cout << "FILE  " << g.G_BLOCK_FILE_PATH << " up to date\n";
        }

        }
        else {
            std::cout << "BLOCK r=" << l << " already complete, skipping...\n";
        }
    }

    goto start;

    return 0;
}
