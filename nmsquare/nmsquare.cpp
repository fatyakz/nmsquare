#define LINE_READ_BUFFER_SIZE 512
//#define checklimits
//#define printchecksquares
//#define countvalids

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
	uint_fast64_t                           matches;
	uint_fast64_t                           n;
	uint_fast64_t                           m;
	uint_fast64_t                           e;
	uint_fast64_t                           r;
};

struct S_thread {
public:
	uint_fast64_t                           id;
	std::chrono::duration<double>           time;
	std::time_t                             date;
	uint_fast64_t                           cycles;
	int                                     step;
	int                                     offset;
	S_best                                  best = {};
};

struct S_block {
public:
	uint_fast64_t                           id;
	std::chrono::duration<double>           time;
	std::time_t                             date;
	uint_fast64_t                           cycles;
	S_best                                  best = {};
	std::vector<S_thread>                   thread;
	uint_fast64_t                           state; 
	std::string                             system_name;
};

struct S_rmw {
	uint_fast64_t                           rmw_reads;
	uint_fast64_t                           rmw_writes;
	uint_fast64_t                           rmw_file_size;
	uint_fast64_t                           rmw_file_blocks;
	uint_fast64_t                           rmw_incomplete;
	uint_fast64_t                           rmw_pending;
	uint_fast64_t                           rmw_completed;
	std::string                             rmw_writemode;
	
};

struct S_var {
	std::string     G_BLOCK_PATH_DEFAULT    = "block.dat";
	uint_fast64_t   T_CYCLES_DIVIDER        = 1;
	std::string     T_CYCLES_SYMBOL         = "";

	uint_fast64_t   B_CYCLES_DIVIDER        = 1000000000;
	std::string     B_CYCLES_SYMBOL         = "B";
	uint_fast64_t   M_CYCLES_DIVIDER        = 1000000;
	std::string     M_CYCLES_SYMBOL         = "M";

	uint_fast64_t   G_TIME_DIVIDER          = 60;
	std::string     G_TIME_SYMBOL           = "m"; 

	uint_fast64_t   H_TIME_DIVIDER          = 360;
	std::string     H_TIME_SYMBOL           = "h";
	uint_fast64_t   M_TIME_DIVIDER          = 60;
	std::string     M_TIME_SYMBOL           = "m";

	uint_fast64_t   G_FILESIZE_DIVIDER      = 1024;
	std::string     G_FILESIZE_SYMBOL       = "kb";
	std::string     T_TIME_SYMBOL           = "s";
	uint_fast64_t   G_PRECISION             = 2;
	uint_fast64_t   G_MIN_WIDTH             = 2;
	uint_fast64_t   G_MIN_WIDTH_NM          = 5;
	std::string     G_COL_SPACE             = "  ";
	uint_fast64_t   G_AVG_CPS_RANGE         = 30;
};

struct S_sizelimits {
	static const uint_fast64_t              ss_uint64 = 18446744073709551615ULL;
	const double                            sl_double = 9007199254740991;
};

struct S_global {
public:
	std::vector<S_block>                    block;
	std::time_t                             date;
	std::chrono::duration<double>           time;
	uint_fast64_t                           cycles;
	S_var                                   var = {};
	S_rmw                                   rmw = {};
	S_best                                  best = {};
	S_sizelimits                            limits;
	uint_fast64_t                           G_BLOCK_START;
	uint_fast64_t                           G_LIMIT;
	std::string                             G_BLOCK_FILE_PATH = "block.dat";
	uint_fast64_t                           G_NUM_THREADS;
	std::string                             G_SYSTEM_NAME;
	uint_fast64_t                           G_MODE;
	bool                                    G_QUIT = 0;
	bool                                    G_CLEAR = 0;
	bool                                    G_BELL = 0;
};

struct S_cell_index {
	uint_fast64_t                           _cell_count = 11;
	uint_fast64_t                           b_state = 0;
	uint_fast64_t                           b_id = 1;
	uint_fast64_t                           b_time = 2;
	uint_fast64_t                           b_date = 3;
	uint_fast64_t                           b_cycles = 4;
	uint_fast64_t                           b_n = 5;
	uint_fast64_t                           b_m = 6;
	uint_fast64_t                           b_e = 7;
	uint_fast64_t                           b_r = 8;
	uint_fast64_t                           b_matches = 9;
	uint_fast64_t                           b_system_name = 10;
};

struct format {
	long double num;
	std::string symbol;
	std::string string;
};

	S_cell_index        cell_index;
	S_global            global;
	S_global            file;

long long check_square(long long n, long long m, long long e) {
	long long a, b, c, d, f, g, h, i;

#ifdef printcheckedsquares 
	bool as{}, bs{}, cs{}, ds{}, es{}, fs{}, gs{}, hs{}, is{};
#endif

	a = e + n;
	b = e - n - m;
	c = e + m;
	d = e - n + m;
	f = e + n - m;
	g = e - m;
	h = e + n + m;
	i = e - n;

	long long matches = 0;

	auto square = [](long double x) {
		if (x > 0) {
			long long sr = sqrt(x);
			return (sr * sr == x);
		} return false;
	};

#ifndef printcheckedsquares
	if (square(a)) { matches++; }
	if (square(b)) { matches++; }
	if (square(c)) { matches++; }
	if (square(d)) { matches++; }
	if (square(e)) { matches++; }
	if (square(f)) { matches++; }
	if (square(g)) { matches++; }
	if (square(h)) { matches++; }
	if (square(i)) { matches++; }
#endif

#ifdef printcheckedsquares
	if (square(a)) { matches++; as = 1; }
	if (square(b)) { matches++; bs = 1; }
	if (square(c)) { matches++; cs = 1; }
	if (square(d)) { matches++; ds = 1; }
	if (square(e)) { matches++; es = 1; }
	if (square(f)) { matches++; fs = 1; }
	if (square(g)) { matches++; gs = 1; }
	if (square(h)) { matches++; hs = 1; }
	if (square(i)) { matches++; is = 1; }

	if (matches > threshold) {
		if (as) { std::cout << "*"; } std::cout << a << ":";
		if (bs) { std::cout << "*"; } std::cout << b << ":";
		if (cs) { std::cout << "*"; } std::cout << c << "\n";
		if (ds) { std::cout << "*"; } std::cout << d << ":";
		if (es) { std::cout << "*"; } std::cout << e << ":";
		if (fs) { std::cout << "*"; } std::cout << f << "\n";
		if (gs) { std::cout << "*"; } std::cout << g << ":";
		if (hs) { std::cout << "*"; } std::cout << h << ":";
		if (is) { std::cout << "*"; } std::cout << i << "\n";
		std::cout << "n:" << n << " m:" << m << " matches:" << matches << "\n\n";
	}
#endif
	return matches;
}

namespace color {
	enum Code {
		FG_RED = 31,
		FG_GREEN = 32,
		FG_BLUE = 34,
		FG_DEFAULT = 39,
		BG_RED = 41,
		BG_GREEN = 42,
		BG_BLUE = 44,
		FG_BLACK = 30,
		FG_YELLOW = 33,
		FG_MAGENTA = 35,
		FG_CYAN = 36,
		FG_LIGHT_GRAY = 37,
		FG_DARK_GRAY = 90,
		FG_LIGHT_RED = 91,
		FG_LIGHT_GREEN = 92,
		FG_LIGHT_YELLOW = 93,
		FG_LIGHT_BLUE = 94,
		FG_LIGHT_MAGENTA = 95,
		FG_LIGHT_CYAN = 96,
		FG_WHITE = 97,
		BG_DEFAULT = 49
	};
	class set {
		Code code;
	public:
		set(Code pCode) : code(pCode) {}
		friend std::ostream&
			operator<<(std::ostream& os, const set& mod) {
			return os << "\033[" << mod.code << "m";
		}
	};
}

struct S_tag {
	std::string S_TAG;
	std::string TAG_COLOR;
	std::string LINE_COLOR;
};

format format_seconds(long double num, S_tag tag) {
	std::ostringstream oss;

	std::string NUM_COLOR = "\033[" + std::to_string(color::FG_DEFAULT) + "m";
	std::string SYM_COLOR = "\033[" + std::to_string(color::FG_DARK_GRAY) + "m";

	format f;

	if (num < 60) { f.num = num; f.symbol = "s"; }
	if (num > 60 && num < 3600) { f.num = num / 60; f.symbol = "m"; }
	if (num > 3600) { f.num = num / 3600; f.symbol = "h"; }

	oss.setf(std::ios::fixed, std::ios::floatfield);
	oss.precision(global.var.G_PRECISION);
	oss << NUM_COLOR << f.num << SYM_COLOR << f.symbol << tag.LINE_COLOR;

	f.string = oss.str();

	return f;
}

format format_long(unsigned long long num, S_tag tag) {
	std::ostringstream oss;

	std::string NUM_COLOR = "\033[" + std::to_string(color::FG_DEFAULT) + "m";
	std::string SYM_COLOR = "\033[" + std::to_string(color::FG_DARK_GRAY) + "m";

	format f;

	if (num < 1000) { f.num = num; f.symbol = ""; }
	if (num > 1000 && num < 1000000) { f.num = num / 1000.0f; f.symbol = "k"; }
	if (num > 1000000 && num < 1000000000) { f.num = num / 1000000.0f; f.symbol = "m"; }
	if (num > 1000000000 && num < 1000000000000) { f.num = num / 1000000000.0f; f.symbol = "b"; }
	if (num > 1000000000000 && num < 1000000000000000) { f.num = num / 1000000000000.0f; f.symbol = "t"; }
	if (num > 1000000000000000) { f.num = num / 1000000000000000.0f; f.symbol = "q"; }

	oss.setf(std::ios::fixed, std::ios::floatfield);
	oss.precision(global.var.G_PRECISION);
	oss << NUM_COLOR << f.num << SYM_COLOR << f.symbol << tag.LINE_COLOR;

	f.string = oss.str();

	return f;
}

format format_filesize(unsigned long long num, S_tag tag) {
	std::ostringstream oss;

	std::string NUM_COLOR = "\033[" + std::to_string(color::FG_DEFAULT) + "m";
	std::string SYM_COLOR = "\033[" + std::to_string(color::FG_DARK_GRAY) + "m";

	format f;

	if (num < 1024) { f.num = num; f.symbol = "b"; }
	if (num > 1024 && num < 1048576) { f.num = num / 1024.0f; f.symbol = "kb"; }
	if (num > 1048576) { f.num = num / 1048576.0f; f.symbol = "mb"; }

	oss.setf(std::ios::fixed, std::ios::floatfield);
	oss.precision(global.var.G_PRECISION);
	oss << NUM_COLOR << f.num << SYM_COLOR << f.symbol << tag.LINE_COLOR;

	f.string = oss.str();

	return f;
}

format format_commas(unsigned long long num, S_tag tag) {
	std::ostringstream oss;
	std::vector<std::string> digits;

	std::string NUM_COLOR = "\033[" + std::to_string(color::FG_DEFAULT) + "m";
	std::string SYM_COLOR = "\033[" + std::to_string(color::FG_DARK_GRAY) + "m";

	format f;

	uint_fast64_t n = num, count = 0;

	while (n > 0) {
		int digit = n % 10;
		n /= 10;
		
		digits.insert(digits.begin(), std::to_string(digit));

		count++;

		if (count == 3 && n > 0) {
			digits.insert(digits.begin(), NUM_COLOR);
			digits.insert(digits.begin(), ",");
			digits.insert(digits.begin(), SYM_COLOR);
			count = 0;
		}
	}

	oss << NUM_COLOR;

	for (auto &i : digits) {
		oss << i;
	}

	oss << tag.LINE_COLOR;

	f.string = oss.str();

	return f;
}

static S_tag tINIT;
static S_tag tSTART;
static S_tag tBLOCK;
static S_tag tPROC;
static S_tag tNMS;
static S_tag tSTATS;
static S_tag tSKIP;
static S_tag tREAD;
static S_tag tWRITE;
static S_tag tERROR;

color::set ly(color::FG_LIGHT_YELLOW);
color::set lr(color::FG_LIGHT_RED);
color::set lg(color::FG_LIGHT_GRAY);
color::set lgrn(color::FG_LIGHT_GREEN);
color::set dg(color::FG_DARK_GRAY);
color::set lcy(color::FG_LIGHT_CYAN);
color::set def(color::FG_DEFAULT);

static void TimeStamp() {
	char            timestamp[50]{ 0 };
	std::time_t     time = std::time(nullptr);

	std::strftime(timestamp, 30, "[%H:%M:%S]", std::localtime(&time));
	std::cout << dg << timestamp << def << global.var.G_COL_SPACE;
}

class C_tag {
public:
	void INIT() {
		tINIT.TAG_COLOR = "\033[" + std::to_string(color::FG_LIGHT_YELLOW) + "m";
		tINIT.LINE_COLOR = "\033[" + std::to_string(color::FG_DEFAULT) + "m";
		tINIT.S_TAG = "INIT ";

		TimeStamp();

		std::cout << tINIT.TAG_COLOR << tINIT.S_TAG << tINIT.LINE_COLOR << global.var.G_COL_SPACE;
	}
	void START() {
		tSTART.TAG_COLOR = "\033[" + std::to_string(color::FG_LIGHT_YELLOW) + "m";
		tSTART.LINE_COLOR = "\033[" + std::to_string(color::FG_DEFAULT) + "m";
		tSTART.S_TAG = "START";

		TimeStamp();

		std::cout << tSTART.TAG_COLOR << tSTART.S_TAG << tSTART.LINE_COLOR << global.var.G_COL_SPACE;
	}
	void BLOCK() {
		tBLOCK.TAG_COLOR = "\033[" + std::to_string(color::FG_LIGHT_CYAN) + "m";
		tBLOCK.LINE_COLOR = "\033[" + std::to_string(color::FG_LIGHT_CYAN) + "m";
		tBLOCK.S_TAG = "BLOCK";

		TimeStamp();

		std::cout << tBLOCK.TAG_COLOR << tBLOCK.S_TAG << tBLOCK.LINE_COLOR << global.var.G_COL_SPACE;
	}
	void PROC() {
		tPROC.TAG_COLOR = "\033[" + std::to_string(color::FG_LIGHT_GRAY) + "m";
		tPROC.LINE_COLOR = "\033[" + std::to_string(color::FG_DARK_GRAY) + "m";
		tPROC.S_TAG = "PROC ";

		TimeStamp();

		std::cout << tPROC.TAG_COLOR << tPROC.S_TAG << tPROC.LINE_COLOR << global.var.G_COL_SPACE;
	}
	static void NMS() {
		tNMS.TAG_COLOR = "\033[" + std::to_string(color::FG_LIGHT_CYAN) + "m";
		tNMS.LINE_COLOR = "\033[" + std::to_string(color::FG_DEFAULT) + "m";
		tNMS.S_TAG = "[nmSquare]";

		std::cout << tNMS.TAG_COLOR << tNMS.S_TAG << tNMS.LINE_COLOR << global.var.G_COL_SPACE;
	}
	void STATS() {
		tSTATS.TAG_COLOR = "\033[" + std::to_string(color::FG_LIGHT_YELLOW) + "m";
		tSTATS.LINE_COLOR = "\033[" + std::to_string(color::FG_LIGHT_YELLOW) + "m";
		tSTATS.S_TAG = "STATS";

		TimeStamp();
		std::cout << tSTATS.TAG_COLOR << tSTATS.S_TAG << tSTATS.LINE_COLOR << global.var.G_COL_SPACE;
	}
	void SKIP() {
		tSKIP.TAG_COLOR = "\033[" + std::to_string(color::FG_LIGHT_CYAN) + "m";
		tSKIP.LINE_COLOR = "\033[" + std::to_string(color::FG_LIGHT_CYAN) + "m";
		tSKIP.S_TAG = "SKIP ";

		TimeStamp();

		std::cout << tSKIP.TAG_COLOR << tSKIP.S_TAG << tSKIP.LINE_COLOR << global.var.G_COL_SPACE;
	}
	void READ() {
		tREAD.TAG_COLOR = "\033[" + std::to_string(color::FG_YELLOW) + "m";
		tREAD.LINE_COLOR = "\033[" + std::to_string(color::FG_DARK_GRAY) + "m";
		tREAD.S_TAG = "READ ";

		TimeStamp();

		std::cout << tREAD.TAG_COLOR << tREAD.S_TAG << tREAD.LINE_COLOR << global.var.G_COL_SPACE;
	}
	void WRITE() {
		tWRITE.TAG_COLOR = "\033[" + std::to_string(color::FG_YELLOW) + "m";
		tWRITE.LINE_COLOR = "\033[" + std::to_string(color::FG_DARK_GRAY) + "m";
		tWRITE.S_TAG = "WRITE";

		TimeStamp();

		std::cout << tWRITE.TAG_COLOR << tWRITE.S_TAG << tWRITE.LINE_COLOR << global.var.G_COL_SPACE;
	}
	void ERROR() {
		tERROR.TAG_COLOR = "\033[" + std::to_string(color::FG_LIGHT_RED) + "m";
		tERROR.LINE_COLOR = "\033[" + std::to_string(color::FG_DEFAULT) + "m";
		tERROR.S_TAG = "ERROR";

		TimeStamp();

		std::cout << tERROR.TAG_COLOR << tERROR.S_TAG << tERROR.LINE_COLOR << global.var.G_COL_SPACE;
	}
	void LIMIT() {
		tWRITE.TAG_COLOR = "\033[" + std::to_string(color::FG_LIGHT_RED) + "m";
		tWRITE.LINE_COLOR = "\033[" + std::to_string(color::FG_DEFAULT) + "m";
		tWRITE.S_TAG = "LIMIT";

		TimeStamp();

		std::cout << tWRITE.TAG_COLOR << tWRITE.S_TAG << tWRITE.LINE_COLOR << global.var.G_COL_SPACE;
	}
};

C_tag tag;

class C_rmw {
public:
	uint_fast64_t                           rmw_reads;
	uint_fast64_t                           rmw_writes;
	uint_fast64_t                           rmw_file_size;
	uint_fast64_t                           rmw_file_blocks;
	uint_fast64_t                           rmw_incomplete;
	uint_fast64_t                           rmw_pending;
	uint_fast64_t                           rmw_completed;

	static void CheckLimits(uint_fast64_t id) {

		double percent = 0;

		tag.LIMIT();

		percent = (global.cycles / global.limits.ss_uint64) * 100;
		std::cout << "[g.c:" << global.cycles << "=" << percent << " ";

		percent = (global.block[id].cycles / global.limits.ss_uint64) * 100;
		std::cout << "b.c:" << global.block[id].cycles << "=" << percent << " ";

		percent = (global.block[id].thread[0].cycles / global.limits.ss_uint64) * 100;
		std::cout << "t.c:" << global.block[id].thread[0].cycles << "=" << percent << " ";

		percent = ((id * id) / global.limits.ss_uint64) * 100;
		std::cout << "e:" << id * id << "=" << percent << "]\n";
  
	}

	static double GetAverageCPS(std::string system, uint_fast64_t range) {
		
		uint_fast64_t                   scanned = 0;
		uint_fast64_t                   t_cycles = 0;
		double                          t_cps = 0;
		std::chrono::duration<double>   t_time = std::chrono::milliseconds::zero();

		for (uint_fast64_t id = global.rmw.rmw_completed; id > 0; id--) {

			if (global.block[id].system_name == system) {

				t_cycles += global.block[id].cycles;
				t_time += global.block[id].time;
				scanned++;
			}

			if (scanned >= range) { 

				break; 
			}
		}

		if (t_cycles > 0 && t_time.count() > 0) {

			t_cps = ((double)t_cycles / t_time.count()) / global.G_NUM_THREADS;
		}

		return t_cps;
	}

	static inline bool FileExists(const std::string& name) {

		struct stat buffer;
		return (stat(name.c_str(), &buffer) == 0);
	}

	static std::vector<std::string> ReadLine(std::istream& str) {

		std::vector<std::string>    result;
		std::string                 line;

		std::getline(str, line);

		std::stringstream           lineStream(line);
		std::string                 cell;

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

		std::streampos  fsize = 0;
		std::ifstream   myfile(path, std::ios::in);

		fsize = myfile.tellg();
		myfile.seekg(0, std::ios::end);
		fsize = myfile.tellg() - fsize;
		myfile.close();

		static_assert(sizeof(fsize) >= sizeof(long long), "Unable to read file size in [GetFileSize()]");

		return fsize;
	}

	uint_fast64_t DigitCount(uint_fast64_t number) {

		return uint_fast64_t(log10(number) + 1);
	}

	static void WriteToFile(std::vector<S_block> block_vector, std::string path) {

		std::ofstream file_out{path, std::ios::trunc};

		if (block_vector.size() > 0) {
			std::vector<std::string> line_buffer;

			for (uint_fast64_t block_index = 0; block_index < block_vector.size(); block_index++) {

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

				for (uint_fast64_t cell_pointer = 1; cell_pointer < cell_index._cell_count; cell_pointer++) {
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

		read_rmw.rmw_writemode = global.rmw.rmw_writemode;

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

			block_buffer[line_index].best.n         = stoll(cells_buffer[cell_index.b_n]);
			block_buffer[line_index].best.m         = stoll(cells_buffer[cell_index.b_m]);
			block_buffer[line_index].best.e         = stoll(cells_buffer[cell_index.b_e]);
			block_buffer[line_index].best.r         = stoll(cells_buffer[cell_index.b_r]);
			block_buffer[line_index].best.matches   = stoll(cells_buffer[cell_index.b_matches]);
			block_buffer[line_index].cycles         = stoll(cells_buffer[cell_index.b_cycles]);
			block_buffer[line_index].date           = stoll(cells_buffer[cell_index.b_date]);
			block_buffer[line_index].id             = stoll(cells_buffer[cell_index.b_id]);
			block_buffer[line_index].state          = stoll(cells_buffer[cell_index.b_state]);
			block_buffer[line_index].time           = seconds_to_duration(stold(cells_buffer[cell_index.b_time]));

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

		for (uint_fast64_t i = 0; i < file.block.size(); i++) {

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

			tag.READ();
			std::cout <<
				"[<-" << global.rmw.rmw_reads << "]" <<
				" blocks:" << format_commas(global.rmw.rmw_file_blocks - 1, tREAD).string <<
				" size:" << format_filesize(global.rmw.rmw_file_size, tREAD).string << 
				" complete:" << format_commas(global.rmw.rmw_completed, tREAD).string <<
				" pending:" << format_commas(global.rmw.rmw_pending, tREAD).string <<
				" incomplete:" << format_commas(global.rmw.rmw_incomplete - 1, tREAD).string << "\n";
		}

		if (global.rmw.rmw_writes > 0) {

			tag.WRITE();
			std::cout <<
				"[->" << global.rmw.rmw_writes << "] " << global.rmw.rmw_writemode << " -> " << path << " (" << 
				format_filesize(global.rmw.rmw_file_size, tREAD).string << ")\n";
			WriteToFile(file.block, path);
			global.rmw.rmw_writemode = "";
		}  
	}
	void Stat() {
		if (global.rmw.rmw_writes > 0 || global.rmw.rmw_reads > 0) {

			double cps = 0;
			if ((double)global.cycles / global.time.count() > 0) {
				cps = (((double)global.cycles) / global.time.count()) / global.G_NUM_THREADS;
			}
			else {
				cps = 0;
			}

			tag.STATS();
			std::cout <<
				"[t:" << format_seconds(global.time.count(), tSTATS).string <<
				" c:" << format_long(global.cycles, tSTATS).string <<
				" cps:" << format_long(cps, tSTATS).string <<
				"] [b:" << format_commas(global.best.matches, tSTATS).string <<
				" n:" << format_commas(global.best.n, tSTATS).string <<
				" m:" << format_commas(global.best.m, tSTATS).string <<
				" e:" << format_commas(global.best.e, tSTATS).string <<
				" r:" << format_commas(global.best.r, tSTATS).string << "]\n";
		}
	}
};

C_rmw rmw;

static std::mutex mlock;
std::vector<std::string> cmd;

void print_square(long long n, long long m, long long e) {
	long long a, b, c, d, f, g, h, i;
	bool as{}, bs{}, cs{}, ds{}, es{}, fs{}, gs{}, hs{}, is{};

	a = e + n;
	b = e - n - m;
	c = e + m;
	d = e - n + m;
	f = e + n - m;
	g = e - m;
	h = e + n + m;
	i = e - n;

	auto square = [](long double x) {
		if (x > 0) {
			long long sr = sqrt(x);
			return (sr * sr == x);
		} return false;
	};

	if (square(a)) { as = 1; }
	if (square(b)) { bs = 1; }
	if (square(c)) { cs = 1; }
	if (square(d)) { ds = 1; }
	if (square(e)) { es = 1; }
	if (square(f)) { fs = 1; }
	if (square(g)) { gs = 1; }
	if (square(h)) { hs = 1; }
	if (square(i)) { is = 1; }

	if (as) { std::cout << "*"; } std::cout << a << ":";
	if (bs) { std::cout << "*"; } std::cout << b << ":";
	if (cs) { std::cout << "*"; } std::cout << c << "\n";
	if (ds) { std::cout << "*"; } std::cout << d << ":";
	if (es) { std::cout << "*"; } std::cout << e << ":";
	if (fs) { std::cout << "*"; } std::cout << f << "\n";
	if (gs) { std::cout << "*"; } std::cout << g << ":";
	if (hs) { std::cout << "*"; } std::cout << h << ":";
	if (is) { std::cout << "*"; } std::cout << i << "\n";

	std::cout << "n:" << n << " m:" << m << "\n\n";

}

static S_thread thr_Single(uint_fast64_t t_E, uint_fast64_t t_offset, uint_fast64_t t_step) {

	uint_fast64_t A, B, C, D, E, F, G, H, I;
	uint_fast64_t t_cycles = 0, t_NM_limit = 0;

	uint_fast64_t t_matches = 0;

	S_thread t_thread{};
	t_thread.best.matches = 0;
	t_thread.id = t_E;
	t_thread.offset = t_offset;
	t_thread.step = t_step;

	E = t_E * t_E;
	t_NM_limit = E - 1;

	auto square = [](uint_fast64_t x) {
		if (x > 0) {
			uint_fast64_t sr = sqrt(x);
			return (sr * sr == x);
		} return false;
	};

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	for (uint_fast64_t lN = 1; lN < t_NM_limit; lN++) {
		for (uint_fast64_t lM = t_offset + 1; lM < t_NM_limit; lM += t_step) {
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

		global.block[t_E].cycles += t_thread.cycles;
		global.block[t_E].thread[t_offset] = t_thread;
		global.block[t_E].id = t_E;
		global.block[t_E].state = 2;

		global.cycles += t_thread.cycles / global.var.B_CYCLES_DIVIDER;

		double cps = (double)t_thread.cycles / t_thread.time.count();

		std::string t_offset_spacer = "";
		if (t_offset < 10) {
			t_offset_spacer = " ";
		}

		if (rmw.DigitCount(global.block[t_E].thread[t_offset].best.n) > global.var.G_MIN_WIDTH_NM) {
			global.var.G_MIN_WIDTH_NM = rmw.DigitCount(global.block[t_E].thread[t_offset].best.n);
		}

		if (rmw.DigitCount(global.block[t_E].thread[t_offset].best.m) > global.var.G_MIN_WIDTH_NM) {
			global.var.G_MIN_WIDTH_NM = rmw.DigitCount(global.block[t_E].thread[t_offset].best.m);
		}

		tag.PROC();
		std::cout << 
			"[r:" << global.block[t_E].thread[t_offset].id << t_offset_spacer << "+" << t_offset << "]" <<
			" cps:" << std::setw(global.var.G_MIN_WIDTH) << cps / global.var.B_CYCLES_DIVIDER << global.var.B_CYCLES_SYMBOL <<
			" c:" << std::setw(global.var.G_MIN_WIDTH) << global.block[t_E].thread[t_offset].cycles << global.var.T_CYCLES_SYMBOL <<
			" t:" << std::setw(global.var.G_MIN_WIDTH) << global.block[t_E].thread[t_offset].time.count() << global.var.T_TIME_SYMBOL <<
			" b:" << global.block[t_E].thread[t_offset].best.matches <<
			" n:" << std::setw(global.var.G_MIN_WIDTH_NM) << global.block[t_E].thread[t_offset].best.n <<
			" m:" << std::setw(global.var.G_MIN_WIDTH_NM) << global.block[t_E].thread[t_offset].best.m <<
			" e:" << global.block[t_E].thread[t_offset].best.e <<
			"\n";

	mlock.unlock();

	return t_thread;
}

static S_thread thr_Nines(uint_fast64_t t_E, uint_fast64_t t_offset, uint_fast64_t t_step) {

	uint_fast64_t A, B, C, D, E, F, G, H, I;
	uint_fast64_t t_cycles = 0, t_NM_limit = 0;

	uint_fast64_t t_matches = 0;

	S_thread t_thread{};
	t_thread.best.matches = 0;
	t_thread.id = t_E;
	t_thread.offset = t_offset;
	t_thread.step = t_step;

	E = t_E * t_E;
	t_NM_limit = E - 1;

	auto square = [](uint_fast64_t x) {
		if (x > 0) {
			uint_fast64_t sr = sqrt(x);
			return (sr * sr == x);
		} return false;
	};

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	for (uint_fast64_t lN = 1; lN < t_NM_limit; lN++) {
		for (uint_fast64_t lM = t_offset + 1; lM < t_NM_limit; lM += t_step) {
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

		global.block[t_E].cycles += t_thread.cycles;
		global.block[t_E].thread[t_offset] = t_thread;
		global.block[t_E].id = t_E;
		global.block[t_E].state = 2;

		global.cycles += t_thread.cycles / global.var.B_CYCLES_DIVIDER;

		double cps = (double)t_thread.cycles / t_thread.time.count();

		std::string t_offset_spacer = "";
		if (t_offset < 10) {
			t_offset_spacer = " ";
		}

		if (rmw.DigitCount(global.block[t_E].thread[t_offset].best.n) > global.var.G_MIN_WIDTH_NM) {
			global.var.G_MIN_WIDTH_NM = rmw.DigitCount(global.block[t_E].thread[t_offset].best.n);
		}

		if (rmw.DigitCount(global.block[t_E].thread[t_offset].best.m) > global.var.G_MIN_WIDTH_NM) {
			global.var.G_MIN_WIDTH_NM = rmw.DigitCount(global.block[t_E].thread[t_offset].best.m);
		}

		TimeStamp();
		std::cout << "PROC " << global.var.G_COL_SPACE << 
			"[r:" << global.block[t_E].thread[t_offset].id << t_offset_spacer << "+" << t_offset << "]" <<
			" cps:" << std::setw(global.var.G_MIN_WIDTH) << cps / global.var.B_CYCLES_DIVIDER << global.var.B_CYCLES_SYMBOL <<
			" c:" << std::setw(global.var.G_MIN_WIDTH) << global.block[t_E].thread[t_offset].cycles << global.var.T_CYCLES_SYMBOL <<
			" t:" << std::setw(global.var.G_MIN_WIDTH) << global.block[t_E].thread[t_offset].time.count() << global.var.T_TIME_SYMBOL <<
			" b:" << global.block[t_E].thread[t_offset].best.matches <<
			" n:" << std::setw(global.var.G_MIN_WIDTH_NM) << global.block[t_E].thread[t_offset].best.n <<
			" m:" << std::setw(global.var.G_MIN_WIDTH_NM) << global.block[t_E].thread[t_offset].best.m <<
			" e:" << global.block[t_E].thread[t_offset].best.e <<
			"\n";

	mlock.unlock();

	return t_thread;
}

static S_thread  thr_nms2(uint_fast32_t start, uint_fast32_t offset, uint_fast32_t threadcount) {
	S_thread t_thread;

	std::mutex tex;
	uint_fast32_t a, b, c, d, e, f, g, h, i;
	uint_fast64_t cycles = 0; uint_fast32_t matches = 0; uint_fast32_t best = 0;
	e = start * start;
	uint_fast32_t nmlimit = e;
	uint_fast32_t bestn = 0, bestm = 0;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	auto isPerfectSquare = [](uint_fast64_t x) {
		if (x > 0) {
			uint_fast64_t sr = sqrt(x);
			return (sr * sr == x);
		} return false;
	};

	for (uint_fast32_t n = 1; n < nmlimit; n++) {
		for (uint_fast32_t m = 1 + offset; m < nmlimit; m += threadcount) {

			if (n + m >= e) { break; }

			a = e + n;
			if (!isPerfectSquare(a)) { goto end; }
			matches++;
			b = e - n - m;
			if (!isPerfectSquare(b)) { goto end; }
			matches++;
			c = e + n;
			if (!isPerfectSquare(c)) { goto end; }
			matches++;
			d = e - n + m;
			if (!isPerfectSquare(d)) { goto end; }
			matches++;
			f = e + n - m;
			if (!isPerfectSquare(f)) { goto end; }
			matches++;
			g = e - n;
			if (!isPerfectSquare(g)) { goto end; }
			matches++;
			h = e + n + m;
			if (!isPerfectSquare(h)) { goto end; }
			matches++;
			i = e - n;
			if (!isPerfectSquare(i)) { goto end; }
			matches++;

		end:

			if (matches > best) {
				best = matches;
				bestn = n; bestm = m;
			}
			matches = 1;
		}
	}

	tex.lock();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> t_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

	cycles = ((e * e * e * e) / 2) / threadcount;

	cycles /= global.var.T_CYCLES_DIVIDER;

	double cps = (double)cycles / t_time.count();

	t_thread.offset = offset;
	t_thread.step = threadcount;
	t_thread.best.matches = best;
	t_thread.best.m = bestm;
	t_thread.best.n = bestn;
	t_thread.best.e = e;
	t_thread.best.r = start;
	t_thread.id = start;
	t_thread.time = t_time;
	auto t_date = std::chrono::system_clock::now();
	t_thread.date = std::chrono::system_clock::to_time_t(t_date);
	t_thread.cycles = cycles;

	global.cycles += t_thread.cycles / global.var.B_CYCLES_DIVIDER;

	global.block[start].cycles += t_thread.cycles / global.var.B_CYCLES_DIVIDER;
	global.block[start].thread[offset] = t_thread;
	global.block[start].id = start;
	global.block[start].state = 2;

	if (t_thread.best.matches >= global.block[start].best.matches) {
		global.block[start].best = t_thread.best;
	}

	std::string t_offset_spacer = "  ";
	if (offset > 9) {
		t_offset_spacer = " ";
	}
	if (offset > 99) {
		t_offset_spacer = "";
	}

	if (rmw.DigitCount(global.block[start].thread[offset].best.n) > global.var.G_MIN_WIDTH_NM) {
		global.var.G_MIN_WIDTH_NM = rmw.DigitCount(global.block[start].thread[offset].best.n);
	}

	if (rmw.DigitCount(global.block[start].thread[offset].best.m) > global.var.G_MIN_WIDTH_NM) {
		global.var.G_MIN_WIDTH_NM = rmw.DigitCount(global.block[start].thread[offset].best.m);
	}

	TimeStamp();
	std::cout << "PROC " << global.var.G_COL_SPACE <<
		"[r:" << global.block[start].thread[offset].id << t_offset_spacer << "+" << offset << "]" <<
		" cps:" << std::setw(global.var.G_MIN_WIDTH) << format_long(cps, tPROC).string <<
		" c:" << std::setw(global.var.G_MIN_WIDTH) << format_long(global.block[start].thread[offset].cycles, tPROC).string <<
		" t:" << std::setw(global.var.G_MIN_WIDTH) << format_seconds(global.block[start].thread[offset].time.count(), tPROC).string <<
		" b:" << global.block[start].thread[offset].best.matches <<
		" n:" << std::setw(global.var.G_MIN_WIDTH_NM) << global.block[start].thread[offset].best.n <<
		" m:" << std::setw(global.var.G_MIN_WIDTH_NM) << global.block[start].thread[offset].best.m <<
		" e:" << global.block[start].thread[offset].best.e <<
		"\n";

	tex.unlock();

	return t_thread;
}

static int thr_find_from_r(long long r, long long offset, long long step) {
	S_best ffr{};
	format form;
	S_thread t_thread;

	long long a, b, c, d, e, f, g, h, i;
	long long n, m;

	e = r * r;
	long long rlimit = sqrt(e) * 2;
#ifdef countvalids 
	long long validnm = 0;
	long long validai = 0;
#endif
	long long matches = 0;
	long long best = 0;
	long long bestn = 0;
	long long bestm = 0;
	long long beste = 0;

	std::chrono::high_resolution_clock::time_point ffr1 = std::chrono::high_resolution_clock::now();

	for (long long fn = 1; fn < rlimit; fn += step) {
		for (long long fm = 1; fm < rlimit; fm++) {

			a = fn * fn;
			c = fm * fm;

			n = a - e;
			m = c - e;

			if (n > 0 && m > 0 && n != m) {

#ifdef countvalids 
				validnm++; 
#endif

				a = e + n;
				b = e - n - m;
				c = e + m;
				d = e - n + m;
				f = e + n - m;
				g = e - m;
				h = e + n + m;
				i = e - n;

				if (a > 0 && b > 0 && c > 0 && d > 0 && f > 0 && g > 0 && h > 0 && i > 0) {

					matches = check_square(n, m, e);
					if (matches >= best) {
						best = matches;
						bestn = n;
						bestm = m;
						beste = e;
					}

#ifdef countvalids 
					validai++; 
#endif

				}
			}
		}
	}

	std::chrono::high_resolution_clock::time_point ffr2 = std::chrono::high_resolution_clock::now();
	auto t_time = std::chrono::duration_cast<std::chrono::duration<double>>(ffr2 - ffr1);

	mlock.lock();

	t_thread.offset = offset;
	t_thread.step = step;
	t_thread.best.matches = best;
	t_thread.best.m = bestm;
	t_thread.best.n = bestn;
	t_thread.best.e = e;
	t_thread.best.r = r;
	t_thread.id = r;
	t_thread.time = t_time;
	auto t_date = std::chrono::system_clock::now();
	t_thread.date = std::chrono::system_clock::to_time_t(t_date);

	t_thread.cycles = ((r + r) - 1) * ((r + r) - 1);

	global.block[r].cycles += t_thread.cycles;
	global.block[r].thread[offset] = t_thread;
	global.block[r].id = r;
	global.block[r].state = 2;

	if (t_thread.best.matches >= global.block[r].best.matches) {
		global.block[r].best = t_thread.best;
	}

	ffr.matches = best;
	ffr.n = bestn;
	ffr.m = bestm;
	ffr.e = beste;
	ffr.r = sqrt(ffr.e);

	double cps = ((double)t_thread.cycles / t_time.count()) / global.G_NUM_THREADS;

	std::string t_offset_spacer = "  ";
	if (offset > 9) {
		t_offset_spacer = " ";
	}
	if (offset > 99) {
		t_offset_spacer = "";
	}

	if (rmw.DigitCount(global.block[r].thread[offset].best.n) > global.var.G_MIN_WIDTH_NM) {
		global.var.G_MIN_WIDTH_NM = rmw.DigitCount(global.block[r].thread[offset].best.n);
	}

	if (rmw.DigitCount(global.block[r].thread[offset].best.m) > global.var.G_MIN_WIDTH_NM) {
		global.var.G_MIN_WIDTH_NM = rmw.DigitCount(global.block[r].thread[offset].best.m);
	}

	tag.PROC();
	std::cout <<
		"[r:" << format_commas(global.block[r].thread[offset].id, tPROC).string << t_offset_spacer << "+" << offset << "]" <<
		" cps:" << format_long(cps, tPROC).string <<
		" t:" << format_seconds(global.block[r].thread[offset].time.count(), tPROC).string <<
		" b:" << format_commas(global.block[r].thread[offset].best.matches, tPROC).string <<
		" n:" << std::setw(global.var.G_MIN_WIDTH_NM) << format_commas(global.block[r].thread[offset].best.n, tPROC).string <<
		" m:" << std::setw(global.var.G_MIN_WIDTH_NM) << format_commas(global.block[r].thread[offset].best.m, tPROC).string <<
		" e:" << format_commas(global.block[r].thread[offset].best.e, tPROC).string <<
		" (" << std::setw(global.var.G_MIN_WIDTH) << format_long(global.block[r].thread[offset].best.e, tPROC).string << ")" <<
		"\n";

	mlock.unlock();
   
#ifdef countvalids
	std::cout << "r:" << r << " cycles:" << cycles << "\n";
	std::cout << "time:" << t_time.count() << "s cps:" << cps / 1000000 << "m \n";
	std::cout << "validnm:" << validnm << " validai:" << validai << "\n";
	std::cout << "best:" << best << " bestn:" << bestn << " bestm:" << bestm << "\n\n";
#endif

	return 1;
}

int read_args(std::vector<std::string> args) {
	bool f{}, s{}, t{}, r{}, b{}, m{};

	for (long unsigned int i = 0; i < args.size(); i++) {
		if (args[i] == "-f") {
			global.G_BLOCK_FILE_PATH = args[i + 1];
			if (!rmw.FileExists(global.G_BLOCK_FILE_PATH)) {
				std::ofstream output(global.G_BLOCK_FILE_PATH);
			}
			f = 1;
		}

		if (args[i] == "-s") {
			global.G_SYSTEM_NAME = args[i + 1];
			s = 1;
		}

		if (args[i] == "-t") {
			global.G_NUM_THREADS = stoll(args[i + 1]);
			t = 1;
		}

		if (args[i] == "-r") {
			global.G_BLOCK_START = stoll(args[i + 1]);
			r = 1;
		}

		if (args[i] == "-b") {
			global.G_LIMIT = stoll(args[i + 1]);
			b = 1;
		}

		if (args[i] == "-m") {
			global.G_MODE = stoll(args[i + 1]);
			m = 1;
		}

		if (args[i] == "-q") {
			global.G_QUIT = 1;
		}

		if (args[i] == "-clear") {
			global.G_CLEAR = 1;
		}

		if (args[i] == "-bell") {
			global.G_BELL = 1;
		}
	}

	if (f && s && t && m) {
		if (r && b) {
			return 0;
		}
		else if (!r && !b) {
			return 1;
		}
		
	}
	else {
		return 2;
	}

	return 3;
}

int main(int argc, char** argv)
{
	const auto processor_count = std::thread::hardware_concurrency();
	std::chrono::high_resolution_clock::time_point g1;
	auto g_date = std::chrono::system_clock::now();

	std::cout << std::setprecision(global.var.G_PRECISION);
	std::cout << std::setw(global.var.G_MIN_WIDTH);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	
	for (int i = 0; i < argc; ++i) {
		cmd.push_back(argv[i]);
	}

reset:

	global.time = std::chrono::milliseconds::zero();
	global.cycles = 0;

	tag.NMS();

	if (cmd.size() > 1) {
		std::cout << dg << "( ";
		for (uint_fast64_t i = 1; i < cmd.size(); ++i) {
			std::cout << cmd[i] << " ";
		}
		std::cout << ")" << def;
	}
	std::cout << "\n";

	if (read_args(cmd) == 0) {
		global.date = std::chrono::system_clock::to_time_t(g_date);
		global.cycles = 0;

		rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);
		cmd.clear();

		goto loophead;
	}
	else if (read_args(cmd) == 1) {
		cmd.clear();
		rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);

		goto start;
	}
	else if (read_args(cmd) == 3) {
		std::cout << "Command line argument error\n";

		return 3;
	}

	tag.INIT();
	std::cout << "Threads" << dg <<" (" << processor_count << " cores):" << def;
	std::cin >> global.G_NUM_THREADS;
	if (global.G_NUM_THREADS == 0) {
		return 0;
	}

	tag.INIT();
	std::cout << "Mode (0:default, 1:nines, 2:nms2, 3:nms3):";
	std::cin >> global.G_MODE;

	tag.INIT();
	std::cout << "INIT " << "System name:";
	std::cin >> global.G_SYSTEM_NAME;

	tag.INIT();
	std::cout << "INIT " << "(" << global.var.G_BLOCK_PATH_DEFAULT << "):";
	std::cin >> global.G_BLOCK_FILE_PATH;
	
	if (!rmw.FileExists(global.G_BLOCK_FILE_PATH)) {
		tag.INIT();
		std::cout << "INIT " << global.G_BLOCK_FILE_PATH << " does not exist, reverting to default: " <<
			global.var.G_BLOCK_PATH_DEFAULT <<"\n";
		global.G_BLOCK_FILE_PATH = global.var.G_BLOCK_PATH_DEFAULT;
	}

	if (!rmw.FileExists(global.G_BLOCK_FILE_PATH)) {
		tag.ERROR();
		std::cout << "!" << global.G_BLOCK_FILE_PATH << " does not exist. Restarting...\n";
		goto reset;
	}

	rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);

pending:

	if (global.rmw.rmw_pending > 0) {
		std::string clear_pending = "n";

		tag.INIT();
		std::cout << "Clear pending? (y/n): "; std::cin >> clear_pending;

		if (clear_pending == "y") {
			for (uint_fast64_t i = 0; i < file.block.size(); i++) {
				if (file.block[i].state == 1) {
					file.block[i].state = 0;
					global.block[i].state = 0;
				}
			}
			rmw.WriteToFile(file.block, global.G_BLOCK_FILE_PATH);

			tag.INIT();
			std::cout << "All pending blocks reset to incomplete\n";
			rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);

			rmw.Stat();
			goto pending;
		}
	}

start:

	rmw.Stat();

	global.time = std::chrono::milliseconds::zero();
	global.cycles = 0;

	tag.INIT();
	std::cout << "Start:"; std::cin >> global.G_BLOCK_START;

	if (global.G_BLOCK_START == 0) {
		goto reset;
	}

	tag.INIT();
	std::cout << "Blocks:"; std::cin >> global.G_LIMIT;
	
	g_date = std::chrono::system_clock::now();
	global.date = std::chrono::system_clock::to_time_t(g_date);

	tag.START();
	std::cout << "[" << global.G_BLOCK_START << "] -> ["
		<< global.G_BLOCK_START + global.G_LIMIT - 1 << "]\n";

loophead:

	cmd.clear();

	for (uint_fast64_t id = global.G_BLOCK_START; id < global.G_BLOCK_START + global.G_LIMIT; id++) {

		g1 = std::chrono::high_resolution_clock::now();

		if (global.block.size() < (id + 1)) { global.block.resize(id+1); }

		rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);
		
		rmw.Stat();

		if (file.block[id].state == 0) {
			S_block g_block;
			g_block.id = id;

			global.block[g_block.id].state = 1;
			global.block[g_block.id].system_name = global.G_SYSTEM_NAME;
			global.rmw.rmw_writemode =
				"[s:" + std::to_string(global.block[g_block.id].state) + " id:" + std::to_string(g_block.id) + " (" + global.G_SYSTEM_NAME + ")]";

			rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);

			uint_fast64_t   predicted_cycles        = (((g_block.id + g_block.id) - 1) * ((g_block.id + g_block.id) - 1)) * global.G_NUM_THREADS;
			double          predicted_cps           = rmw.GetAverageCPS(global.G_SYSTEM_NAME, global.var.G_AVG_CPS_RANGE);
			double          predicted_seconds       = (predicted_cycles / predicted_cps) / global.G_NUM_THREADS;
			uint_fast64_t   predicted_total_cycles  = 0;
			double          predicted_total_seconds = 0.0f;

			for (uint_fast64_t i = g_block.id; i < global.G_BLOCK_START + global.G_LIMIT; i++) {
				predicted_total_cycles += (((i + i) - 1) * ((i + i) - 1)) * global.G_NUM_THREADS;
				predicted_total_seconds += (predicted_total_cycles / rmw.GetAverageCPS(global.G_SYSTEM_NAME, global.var.G_AVG_CPS_RANGE))
					/ global.G_NUM_THREADS;
			}

			tag.BLOCK();
			std::cout << "[" << format_commas(g_block.id, tBLOCK).string << "/" << 
				format_commas(global.G_BLOCK_START + global.G_LIMIT - 1, tBLOCK).string <<
				" [avg cps:" << format_long(predicted_cps, tBLOCK).string << 
				"(" << global.var.G_AVG_CPS_RANGE << ")]" <<
				" [est c:" << format_long(predicted_cycles, tBLOCK).string <<
				"/" << format_long(predicted_total_cycles, tBLOCK).string <<
				" t:" << format_seconds(predicted_seconds, tBLOCK).string <<
				"/" << format_seconds(predicted_total_seconds, tBLOCK).string <<
				"] PENDING...\n";

			global.block[g_block.id].thread.resize(global.G_NUM_THREADS);

			auto b_date = std::chrono::system_clock::now();
			global.block[g_block.id].date = std::chrono::system_clock::to_time_t(b_date);

			std::chrono::high_resolution_clock::time_point b1 = std::chrono::high_resolution_clock::now();

			std::vector<std::thread> thr(global.G_NUM_THREADS);

			if (global.G_MODE == 0) {
				for (uint_fast64_t g_thread_offset = 0; g_thread_offset < global.G_NUM_THREADS; g_thread_offset++) {
					thr[g_thread_offset] = std::thread(thr_Single, g_block.id, g_thread_offset, global.G_NUM_THREADS);
				}
				for (uint_fast64_t g_thread_id = 0; g_thread_id < global.G_NUM_THREADS; g_thread_id++) {
					thr[g_thread_id].join();
				}
			}
			if (global.G_MODE == 1) {
				for (uint_fast64_t g_thread_offset = 0; g_thread_offset < global.G_NUM_THREADS; g_thread_offset++) {
					thr[g_thread_offset] = std::thread(thr_Nines, g_block.id, g_thread_offset, global.G_NUM_THREADS);
				}
				for (uint_fast64_t g_thread_id = 0; g_thread_id < global.G_NUM_THREADS; g_thread_id++) {
					thr[g_thread_id].join();
				}
			}
			if (global.G_MODE == 2) {
				for (uint_fast64_t g_thread_offset = 0; g_thread_offset < global.G_NUM_THREADS; g_thread_offset++) {
					thr[g_thread_offset] = std::thread(thr_nms2, g_block.id, g_thread_offset, global.G_NUM_THREADS);
				}
				for (uint_fast64_t g_thread_id = 0; g_thread_id < global.G_NUM_THREADS; g_thread_id++) {
					thr[g_thread_id].join();
				}
			}
			if (global.G_MODE == 3) {
				for (uint_fast64_t g_thread_offset = 0; g_thread_offset < global.G_NUM_THREADS; g_thread_offset++) {
					thr[g_thread_offset] = std::thread(thr_find_from_r, g_block.id, g_thread_offset, global.G_NUM_THREADS);
				}
				for (uint_fast64_t g_thread_id = 0; g_thread_id < global.G_NUM_THREADS; g_thread_id++) {
					thr[g_thread_id].join();
				}
			}

			if (global.G_CLEAR) {
				if (!std::system("clear")) {
					std::cout << "";
				}
			}

			std::chrono::high_resolution_clock::time_point b2 = std::chrono::high_resolution_clock::now();
			global.block[g_block.id].time = std::chrono::duration_cast<std::chrono::duration<double>>(b2 - b1);

			if (global.block[g_block.id].best.matches >= global.best.matches) {
				global.best = global.block[g_block.id].best;
			}

			double cps = ((double)global.block[g_block.id].cycles / global.block[g_block.id].time.count()) / global.G_NUM_THREADS;

			global.cycles += global.block[g_block.id].cycles;

			tag.BLOCK();
			std::cout << "[" << format_commas(g_block.id, tBLOCK).string << "] COMPLETE" <<
				" [c:"  << format_long(global.block[g_block.id].cycles, tBLOCK).string <<
				" t:"    << format_seconds(global.block[g_block.id].time.count(), tBLOCK).string <<
				" cps:"     << format_long(cps, tBLOCK).string << "]\n";

			rmw.Stat();

			global.block[g_block.id].state = 2;
			global.rmw.rmw_writemode = 
				"[s:" + std::to_string(global.block[g_block.id].state) + " id:" + std::to_string(g_block.id) + 
				" (" + global.G_SYSTEM_NAME + ")]";

			std::chrono::high_resolution_clock::time_point g2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> g_total_time = std::chrono::duration_cast<std::chrono::duration<double>>(g2 - g1);

			global.time += g_total_time;
			

#ifdef checklimits
			rmw.CheckLimits(g_block.id);
#endif
			
		} 
		else {
			if (global.block[id].state == 1) {
				tag.SKIP();
				std::cout << "[" << id << "] PENDING  ("
					<< global.block[id].system_name << "), skipping...\n" << def;
			}
			if (global.block[id].state == 2) {
				tag.SKIP();
				std::cout << "[" << id << "] COMPLETE ("
					<< global.block[id].system_name << "), skipping...\n" << def;
			}
		}
	}
	rmw.ReadMergeWrite(global.G_BLOCK_FILE_PATH);

	if (global.G_BELL) {
		std::cout << "\a";
	}

	if (global.G_QUIT) {
		return 0;
	}

	goto start;

	return 0;
}