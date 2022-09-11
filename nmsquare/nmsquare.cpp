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
#include <sys/stat.h>
#include <filesystem>
#include <regex>

struct S_best {
public:
    uint_fast32_t                       matches;
    uint_fast32_t                           n;
    uint_fast32_t                           m;
    uint_fast32_t                           e;
    uint_fast32_t                           r;
};
struct S_thread {
public:
    uint_fast32_t                           id;
    std::chrono::duration<double>           time;
    std::time_t                             date;
    unsigned long long int                  cycles;
    int                                     step;
    int                                     offset;
    S_best                                  best = {};
};
struct S_block {
public:
    unsigned long long int                  id;
    std::chrono::duration<double>           time;
    std::time_t                             date;
    unsigned long long int                  cycles;
    S_best                                  best = {};
    std::vector<S_thread>                   thread;
    uint_fast32_t                           state; 
    std::string                             system_name;
};
struct S_rmw {
    uint_fast32_t                           rmw_reads;
    uint_fast32_t                           rmw_writes;
    uint_fast32_t                           rmw_file_size;
    uint_fast32_t                           rmw_file_blocks;
    uint_fast32_t                           rmw_incomplete;
    uint_fast32_t                           rmw_pending;
    uint_fast32_t                           rmw_completed;
};
struct S_var {
    std::string G_BLOCK_FILE_PATH_DEFAULT   = "block.dat";
    uint_fast32_t T_CYCLES_DIVIDER          = 1;
    std::string T_CYCLES_SYMBOL             = "";
    uint_fast32_t B_CYCLES_DIVIDER          = 1000000;
    std::string B_CYCLES_SYMBOL             = "M";
    uint_fast32_t G_CYCLES_DIVIDER          = 1000;
    std::string G_CYCLES_SYMBOL             = "B";
    uint_fast32_t G_TIME_DIVIDER            = 60;
    std::string G_TIME_SYMBOL               = "m";
    uint_fast32_t G_FILESIZE_DIVIDER        = 1024;
    std::string G_FILESIZE_SYMBOL           = "kb";
};
struct S_global {
public:
    std::vector<S_block>                    block;
    std::time_t                             date;
    std::chrono::duration<double>           time;
    unsigned long long int                  cycles;
    S_var                                   var = {};
    S_rmw                                   rmw = {};
    S_best                                  best = {};
    uint_fast32_t                           G_BLOCK_START;
    uint_fast32_t                           G_LIMIT;
    std::string                             G_BLOCK_FILE_PATH = "block.dat";
    uint_fast32_t                           G_NUM_THREADS;
    std::string                             G_SYSTEM_NAME;
    uint_fast32_t                           G_MODE;
};
struct S_cell_index {
    uint_fast32_t                           _cell_count = 11;
    uint_fast32_t                           b_state = 0;
    uint_fast32_t                           b_id = 1;
    uint_fast32_t                           b_time = 2;
    uint_fast32_t                           b_date = 3;
    uint_fast32_t                           b_cycles = 4;
    uint_fast32_t                           b_n = 5;
    uint_fast32_t                           b_m = 6;
    uint_fast32_t                           b_e = 7;
    uint_fast32_t                           b_r = 8;
    uint_fast32_t                           b_matches = 9;
    uint_fast32_t                           b_system_name = 10;
};


    S_cell_index    cell_index;
    S_global        global;
    S_global        file;

class C_rmw {
public:
    uint_fast32_t                           rmw_reads;
    uint_fast32_t                           rmw_writes;
    uint_fast32_t                           rmw_file_size;
    uint_fast32_t                           rmw_file_blocks;
    uint_fast32_t                           rmw_incomplete;
    uint_fast32_t                           rmw_pending;
    uint_fast32_t                           rmw_completed;

    static inline bool FileExists(const std::string& name) {
        struct stat buffer;
        return (stat(name.c_str(), &buffer) == 0);
    }
    static void TimeStamp() {
        char timestamp[50]{ 0 };
        std::time_t time = std::time(nullptr);
        std::strftime(timestamp, 30, "[%H:%M:%S] ", std::localtime(&time));
        std::cout << timestamp;
    }
    static std::vector<std::string> ReadLine(std::istream& str) {
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
    template <typename T>
    static auto seconds_to_duration(T seconds) {
        return std::chrono::duration<T, std::ratio<1>>(seconds);
    }
    static long long GetFileSize(std::string path) {
        std::streampos fsize = 0;
        std::ifstream myfile(path, std::ios::in);
        fsize = myfile.tellg();
        myfile.seekg(0, std::ios::end);
        fsize = myfile.tellg() - fsize;
        myfile.close();
        static_assert(sizeof(fsize) >= sizeof(long long), "Oops.");
        return fsize;
    }
    static void WriteToFile(std::vector<S_block> block_vector, std::string path) {
        std::ofstream file_out{ path, std::ios::trunc };

        if (block_vector.size() > 0) {
            std::vector<std::string> line_buffer;

            for (uint_fast32_t block_index = 0; block_index < block_vector.size(); block_index++) {

                line_buffer.resize(cell_index._cell_count);

                line_buffer[cell_index.b_cycles]        = std::to_string(block_vector[block_index].cycles);
                line_buffer[cell_index.b_date]          = std::to_string(block_vector[block_index].date);
                line_buffer[cell_index.b_n]             = std::to_string(block_vector[block_index].best.n);
                line_buffer[cell_index.b_m]             = std::to_string(block_vector[block_index].best.m);
                line_buffer[cell_index.b_e]             = std::to_string(block_vector[block_index].best.e);
                line_buffer[cell_index.b_r]             = std::to_string(block_vector[block_index].best.r);
                line_buffer[cell_index.b_matches]       = std::to_string(block_vector[block_index].best.matches);
                line_buffer[cell_index.b_id]            = std::to_string(block_vector[block_index].id);
                line_buffer[cell_index.b_state]         = std::to_string(block_vector[block_index].state);
                line_buffer[cell_index.b_time]          = std::to_string(block_vector[block_index].time.count());
                line_buffer[cell_index.b_system_name]   = block_vector[block_index].system_name;

                file_out << line_buffer[0];

                for (uint_fast32_t cell_pointer = 1; cell_pointer < cell_index._cell_count; cell_pointer++) {
                    file_out << ",";
                    file_out << line_buffer[cell_pointer];
                }

                if (block_index < block_vector.size() - 1) { file_out << "\n"; }

                line_buffer.resize(0);
            }
        }
        return;
    }


    static S_rmw ReadFromFile(std::string path, std::vector<S_block>& target_block) {

        S_rmw read_rmw = {};

        char line_buffer[LINE_READ_BUFFER_SIZE];
        int line_index = 0;

        std::ifstream file_out{ path, std::ifstream::binary };
        std::vector<S_block> block_buffer;
        std::string s; int file_lines = 0;
        std::ifstream in; in.open(path);

        while (!in.eof()) { getline(in, s); file_lines++; }

        block_buffer.resize(file_lines);


        while (file_out.getline(line_buffer, LINE_READ_BUFFER_SIZE)) {
            std::istringstream csv_buffer(line_buffer);
            std::vector<std::string> cells_buffer = ReadLine(csv_buffer);

            block_buffer[line_index].best.n         = stoi(cells_buffer[cell_index.b_n]);
            block_buffer[line_index].best.m         = stoi(cells_buffer[cell_index.b_m]);
            block_buffer[line_index].best.e         = stoi(cells_buffer[cell_index.b_e]);
            block_buffer[line_index].best.r         = stoi(cells_buffer[cell_index.b_r]);
            block_buffer[line_index].best.matches   = stoi(cells_buffer[cell_index.b_matches]);
            block_buffer[line_index].cycles         = stol(cells_buffer[cell_index.b_cycles]);
            block_buffer[line_index].date           = stoi(cells_buffer[cell_index.b_date]);
            block_buffer[line_index].id             = stoi(cells_buffer[cell_index.b_id]);
            block_buffer[line_index].state          = stoi(cells_buffer[cell_index.b_state]);
            block_buffer[line_index].time           = seconds_to_duration(stod(cells_buffer[cell_index.b_time]));

            std::string sys_name = cells_buffer[cell_index.b_system_name];

            sys_name = std::regex_replace(sys_name, std::regex("\\r\\n|\\r|\\n"), "");

            block_buffer[line_index].system_name = sys_name;

            line_index++;
        }

        target_block.resize(block_buffer.size());
        target_block = block_buffer;

        read_rmw.rmw_file_size = GetFileSize(global.G_BLOCK_FILE_PATH);
        read_rmw.rmw_file_blocks = file_lines;

        return read_rmw;
    }
    static void ReadMergeWrite(std::string path) {

        global.rmw = ReadFromFile(path, file.block);

        if (file.block.size() > global.block.size()) {
            global.block.resize(file.block.size());
        }
        else {
            file.block.resize(global.block.size());
        }

        global.rmw.rmw_reads = 0;
        global.rmw.rmw_writes = 0;
        global.rmw.rmw_incomplete = 0;
        global.rmw.rmw_pending = 0;
        global.rmw.rmw_completed = 0;

        for (unsigned long long int i = 0; i < file.block.size(); i++) {


            if (file.block[i].state != global.block[i].state) {
                if (file.block[i].state > global.block[i].state) {
                    global.block[i] = file.block[i];
                    global.rmw.rmw_reads += 1;
                }
                else {
                    file.block[i] = global.block[i];
                    global.rmw.rmw_writes += 1;
                }
            }

            if (global.block[i].best.matches >= global.best.matches) { global.best = global.block[i].best; }

            if (file.block[i].state == 0) { global.rmw.rmw_incomplete += 1; }
            if (file.block[i].state == 1) { global.rmw.rmw_pending += 1; }
            if (file.block[i].state == 2) { global.rmw.rmw_completed += 1; }
        }
        if (global.rmw.rmw_reads > 0) {

            TimeStamp();
            std::cout <<
                "READ " <<
                " reads=" << global.rmw.rmw_reads <<
                " blocks=" << global.rmw.rmw_file_blocks - 1 <<
                " size=" << global.rmw.rmw_file_size / global.var.G_FILESIZE_DIVIDER << global.var.G_FILESIZE_SYMBOL <<
                " complete=" << global.rmw.rmw_completed <<
                " pending=" << global.rmw.rmw_pending <<
                " incomplete=" << global.rmw.rmw_incomplete - 1 << "\n";
        }

        if (global.rmw.rmw_writes > 0) {

            TimeStamp();
            std::cout << "WRITE writes=" << global.rmw.rmw_writes << " -> " << path << " as [" << global.G_SYSTEM_NAME << "]\n";
            WriteToFile(file.block, path);
        }

        
    }
    void Stat() {
        if (global.rmw.rmw_writes > 0 || global.rmw.rmw_reads > 0) {

            double cps = 0;
            if (((double)global.cycles * global.var.G_CYCLES_DIVIDER) / global.time.count() > 0) {
                cps = ((double)global.cycles * global.var.G_CYCLES_DIVIDER) / global.time.count();
            }
            else {
                cps = 0;
            }

            TimeStamp();
            std::cout <<
                "STAT " <<
                " time=" << global.time.count() / global.var.G_TIME_DIVIDER << global.var.G_TIME_SYMBOL <<
                " cycles=" << global.cycles << global.var.G_CYCLES_SYMBOL <<
                " cps=" << cps << global.var.B_CYCLES_SYMBOL <<
                " best=" << global.best.matches <<
                " n=" << global.best.n <<
                " m=" << global.best.m <<
                " e=" << global.best.e <<
                " r=" << global.best.r << "\n";
        }
    }
};

C_rmw rmw;

static std::mutex mlock;

static S_thread thr_Single(unsigned long long int t_E, uint_fast32_t t_offset, uint_fast32_t t_step) {

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

    auto square = [](unsigned long long int x) {
        if (x > 0) {
            unsigned long long int sr = sqrt(x);
            return (sr * sr == x);
        } return false;
    };

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

            if (t_thread.best.matches >= global.block[t_E].best.matches) {
                global.block[t_E].best = t_thread.best;
            }

            t_matches = 0;
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

    auto t_date = std::chrono::system_clock::now();
    t_thread.date = std::chrono::system_clock::to_time_t(t_date);
    t_cycles /= global.var.T_CYCLES_DIVIDER;
    t_thread.cycles = t_cycles;
    t_thread.time = t_time;

    mlock.lock();
    global.block[t_E].cycles += t_thread.cycles / global.var.B_CYCLES_DIVIDER;
    global.block[t_E].thread[t_offset] = t_thread;
    global.block[t_E].id = t_E;
    global.block[t_E].state = 2;

    rmw.TimeStamp();
    double cps = (double)t_thread.cycles / t_thread.time.count();

    std::cout << "PROC  [" << global.block[t_E].thread[t_offset].id << "+" << t_offset << "]" <<
        " cycles=" << global.block[t_E].thread[t_offset].cycles << global.var.T_CYCLES_SYMBOL <<
        " t=" << global.block[t_E].thread[t_offset].time.count() <<
        " best=" << global.block[t_E].thread[t_offset].best.matches <<
        " n=" << global.block[t_E].thread[t_offset].best.n <<
        " m=" << global.block[t_E].thread[t_offset].best.m <<
        " e=" << global.block[t_E].thread[t_offset].best.e <<
        " r=" << global.block[t_E].thread[t_offset].best.r <<
        " cps=" << cps / global.var.B_CYCLES_DIVIDER << global.var.B_CYCLES_SYMBOL <<
        "\n";
    mlock.unlock();

    return t_thread;
}

static S_thread thr_Nines(unsigned long long int t_E, uint_fast32_t t_offset, uint_fast32_t t_step) {

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

    auto square = [](unsigned long long int x) {
        if (x > 0) {
            unsigned long long int sr = sqrt(x);
            return (sr * sr == x);
        } return false;
    };

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
            if (square(A) != true) { goto end; }
            
            t_matches++;
            B = E - lN - lM;
            if (square(B) != true) { goto end; }

            t_matches++;
            C = E + lM;
            if (square(C) != true) { goto end; }

            t_matches++;
            D = E - lN + lM;
            if (square(D) != true) { goto end; }

            t_matches++;
            F = E + lN - lM;
            if (square(F) != true) { goto end; }

            t_matches++;
            G = E - lM;
            if (square(G) != true) { goto end; }

            t_matches++;
            H = E + lN + lM;
            if (square(H) != true) { goto end; }

            t_matches++;
            I = E - lN;
            if (square(I) != true) { goto end; }

            t_matches++;

end:
            t_cycles++;
            if (t_matches >= t_thread.best.matches) {
                t_thread.best.matches = t_matches + 1;
                t_thread.best.e = E;
                t_thread.best.n = lN;
                t_thread.best.m = lM;
                t_thread.best.r = t_E;
            }

            if (t_thread.best.matches >= global.block[t_E].best.matches) {
                global.block[t_E].best = t_thread.best;
            }

            t_matches = 0;
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

    auto t_date = std::chrono::system_clock::now();
    t_thread.date = std::chrono::system_clock::to_time_t(t_date);
    t_cycles /= global.var.T_CYCLES_DIVIDER;
    t_thread.cycles = t_cycles;
    t_thread.time = t_time;

    mlock.lock();
    global.block[t_E].cycles += t_thread.cycles / global.var.B_CYCLES_DIVIDER;
    global.block[t_E].thread[t_offset] = t_thread;
    global.block[t_E].id = t_E;
    global.block[t_E].state = 2;

    rmw.TimeStamp();
    double cps = (double)t_thread.cycles / t_thread.time.count();

    std::cout << "PROC  [" << global.block[t_E].thread[t_offset].id << "+" << t_offset << "]" <<
        " cycles=" << global.block[t_E].thread[t_offset].cycles << global.var.T_CYCLES_SYMBOL <<
        " t=" << global.block[t_E].thread[t_offset].time.count() <<
        " best=" << global.block[t_E].thread[t_offset].best.matches <<
        " n=" << global.block[t_E].thread[t_offset].best.n <<
        " m=" << global.block[t_E].thread[t_offset].best.m <<
        " e=" << global.block[t_E].thread[t_offset].best.e <<
        " r=" << global.block[t_E].thread[t_offset].best.r <<
        " cps=" << cps / global.var.B_CYCLES_DIVIDER << global.var.B_CYCLES_SYMBOL <<
        "\n";
    mlock.unlock();

    return t_thread;
}

int main()
{
    const auto processor_count = std::thread::hardware_concurrency();

reset:
    global.time = std::chrono::milliseconds::zero();
    global.cycles = 0;

    std::cout << "[nmSquare]\n";
    rmw.TimeStamp();
    std::cout << "INIT  Threads (" << processor_count << " cores):";
    std::cin >> global.G_NUM_THREADS;
    if (global.G_NUM_THREADS == 0) {
        return 0;
    }

    rmw.TimeStamp();
    std::cout << "INIT  Mode (0=default, 1=nines):";
    std::cin >> global.G_MODE;

    rmw.TimeStamp();
    std::cout << "INIT  System name:";
    std::cin >> global.G_SYSTEM_NAME;

    rmw.TimeStamp();
    std::cout << "INIT  (" << global.var.G_BLOCK_FILE_PATH_DEFAULT << "):";
    std::cin >> global.G_BLOCK_FILE_PATH;
    
    if (!rmw.FileExists(global.G_BLOCK_FILE_PATH)) {
        rmw.TimeStamp();
        std::cout << "INIT  " << global.G_BLOCK_FILE_PATH << " does not exist, reverting to default: " << global.var.G_BLOCK_FILE_PATH_DEFAULT <<"\n";
        global.G_BLOCK_FILE_PATH = global.var.G_BLOCK_FILE_PATH_DEFAULT;
    }

    if (!rmw.FileExists(global.G_BLOCK_FILE_PATH)) {
        rmw.TimeStamp();
        std::cout << "ERROR !" << global.G_BLOCK_FILE_PATH << " does not exist. Restarting...";
        goto reset;
    }

    rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);

    if (global.rmw.rmw_pending > 0) {
        std::string clear_pending = "n";

        rmw.TimeStamp();
        std::cout << "INIT  Clear pending? (y/n): "; std::cin >> clear_pending;

        if (clear_pending == "y") {
            for (uint_fast32_t i = 0; i < file.block.size(); i++) {
                if (file.block[i].state == 1) {
                    file.block[i].state = 0;
                    global.block[i].state = 0;
                }
            }
            rmw.WriteToFile(file.block, global.G_BLOCK_FILE_PATH);

            rmw.TimeStamp();
            std::cout << "INIT  All pending blocks reset to incomplete\n";
            rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);
            goto reset;
        }
    }

start:
    global.time = std::chrono::milliseconds::zero();
    global.cycles = 0;

    rmw.TimeStamp();
    std::cout << "INIT  Start:"; std::cin >> global.G_BLOCK_START;

    if (global.G_BLOCK_START == 0) {
        goto reset;
    }
    rmw.TimeStamp();
    std::cout << "INIT  Blocks:"; std::cin >> global.G_LIMIT;
    std::chrono::high_resolution_clock::time_point g1 = std::chrono::high_resolution_clock::now();

    auto g_date = std::chrono::system_clock::now();
    global.date = std::chrono::system_clock::to_time_t(g_date);

    if (global.block.size() < (global.G_BLOCK_START + global.G_LIMIT)) { global.block.resize(global.G_BLOCK_START + global.G_LIMIT); }

    rmw.TimeStamp();
    std::wcout << "BLOCK [" << global.G_BLOCK_START << "] -> [" << global.G_BLOCK_START + global.G_LIMIT - 1 << "]\n";

    for (uint_fast32_t id = global.G_BLOCK_START; id < global.G_BLOCK_START + global.G_LIMIT; id++) {

        rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);
        
        rmw.Stat();

        if (global.block[id].state == 0) {
            S_block g_block;
            g_block.id = id;

            global.block[g_block.id].state = 1;
            global.block[g_block.id].system_name = global.G_SYSTEM_NAME;

            rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);

            rmw.TimeStamp();
            std::cout << "BLOCK [" << g_block.id << " of " << global.G_BLOCK_START + global.G_LIMIT - 1 <<
                "] threads=" << global.G_NUM_THREADS << " flagged as PENDING...\n";

            global.block[g_block.id].thread.resize(global.G_NUM_THREADS);

            auto b_date = std::chrono::system_clock::now();
            global.block[g_block.id].date = std::chrono::system_clock::to_time_t(b_date);

            std::chrono::high_resolution_clock::time_point b1 = std::chrono::high_resolution_clock::now();

            std::vector<std::thread> thr(global.G_NUM_THREADS);

            if (global.G_MODE == 0) {
                for (uint_fast32_t g_thread_offset = 0; g_thread_offset < global.G_NUM_THREADS; g_thread_offset++) {
                    thr[g_thread_offset] = std::thread(thr_Single, g_block.id, g_thread_offset, global.G_NUM_THREADS);
                }
                for (uint_fast32_t g_thread_id = 0; g_thread_id < global.G_NUM_THREADS; g_thread_id++) {
                    thr[g_thread_id].join();
                }
            }
            if (global.G_MODE == 1) {
                for (uint_fast32_t g_thread_offset = 0; g_thread_offset < global.G_NUM_THREADS; g_thread_offset++) {
                    thr[g_thread_offset] = std::thread(thr_Nines, g_block.id, g_thread_offset, global.G_NUM_THREADS);
                }
                for (uint_fast32_t g_thread_id = 0; g_thread_id < global.G_NUM_THREADS; g_thread_id++) {
                    thr[g_thread_id].join();
                }
            }

            std::chrono::high_resolution_clock::time_point b2 = std::chrono::high_resolution_clock::now();
            global.block[g_block.id].time = std::chrono::duration_cast<std::chrono::duration<double>>(b2 - b1);

            global.cycles += global.block[g_block.id].cycles / global.var.G_CYCLES_DIVIDER;

            if (global.block[g_block.id].best.matches >= global.best.matches) {
                global.best = global.block[g_block.id].best;
            }

            rmw.TimeStamp();
            double cps = (double)global.block[g_block.id].cycles / global.block[g_block.id].time.count(); 
            std::cout << "BLOCK [" << g_block.id << "] flagged as COMPLETE" << 
                " cycles="  << global.block[g_block.id].cycles << global.var.B_CYCLES_SYMBOL <<
                " time="    << global.block[g_block.id].time.count() / global.var.G_TIME_DIVIDER << global.var.G_TIME_SYMBOL <<
                " cps="    << cps << global.var.B_CYCLES_SYMBOL << "\n";

            rmw.Stat();

            global.block[g_block.id].state = 2;

            std::chrono::high_resolution_clock::time_point g2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> g_total_time = std::chrono::duration_cast<std::chrono::duration<double>>(g2 - g1);

            global.time += g_total_time;

        } 
        else {
            if (global.block[id].state == 1) {
                rmw.TimeStamp();
                std::cout << "BLOCK [" << id << "] pending (" << global.block[id].system_name << ") skipping...\n";
            }
            if (global.block[id].state == 2) {
                rmw.TimeStamp();
                std::cout << "BLOCK [" << id << "] completed (" << global.block[id].system_name << ") skipping...\n";
            }
        }
    }
    rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);

    goto start;

    return 0;
}
