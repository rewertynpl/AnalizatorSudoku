
// author copyright Marcin Matysek (rewertynPL)

//g++ sudoku_analyzer-v1.1.cpp -o sudoku_analyzer-v1.1.exe -O2 -std=c++17 -static -lole32 -lshell32 -luuid

#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <shobjidl.h>
#include <shlobj.h>
#include <algorithm>
#include <array>
#include <cctype>
#include <climits>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#pragma comment(lib, "ole32.lib")
#pragma comment(lib, "shell32.lib")
#pragma comment(lib, "uuid.lib")

namespace fs = std::filesystem;

struct Cell { int value = 0; bool revealed = false; };
struct SudokuBoard {
    long long seed = 0;
    int block_rows = 0, block_cols = 0, side_size = 0, total_cells = 0;
    std::vector<Cell> cells;
    bool valid = false;
    std::string error;
};

enum class Strategy {
    NakedSingle, HiddenSingle,
    NakedPair, HiddenPair,
    PointingPairsTriples, BoxLineReduction,
    NakedTriple, HiddenTriple,
    NakedQuad, HiddenQuad,
    XWing, YWing, XYZWing, WXYZWing, Swordfish, Jellyfish, FrankenMutantFish, KrakenFish,
    Skyscraper, TwoStringKite, SimpleColoring, ThreeDMedusa, FinnedXWingSashimi, FinnedSwordfish, FinnedJellyfish, EmptyRectangle,
    UniqueRectangle, UniqueLoop, BivalueOddagon, AvoidableRectangle, BUGPlus1,
    RemotePairs, WWing, GroupedXCycle, XChain, XYChain, GroupedAIC, AIC, ContinuousNiceLoop,
    ALSXZ, ALSXYWing, ALSChain, DeathBlossom, SueDeCoq, MSLS, Exocet, SeniorExocet, SKLoop, PatternOverlayMethod,
    ForcingChains, Backtracking
};

struct AnalysisReport {
    bool contradiction = false;
    bool solved_logically = false;
    bool requires_guessing = false;
    bool solved_with_backtracking = false;
    bool unique_solution = false;
    int solution_count = 0;
    int initial_clues = 0;
    int hardest_rank = 0;
    long long backtracking_nodes = 0;
    long long backtracking_decisions = 0;
    long long backtracking_backtracks = 0;
    std::map<Strategy, int> strategy_usage;
    std::string hardest_strategy;
    std::vector<std::string> debug_logic_logs;
};

struct PuzzleReportEntry {
    std::string source_file;
    int line_no = 0;
    bool valid = false;
    std::string sudoku_type = "Nieznany";
    std::string board_type = "Nieznany";
    std::string parse_error;
    int initial_clues = 0;
    int difficulty_level = 0;
    bool solved_logically = false;
    bool requires_guessing = false;
    bool solved_with_backtracking = false;
    bool contradiction = false;
    int solution_count = 0;
    long long backtracking_nodes = 0;
    long long backtracking_decisions = 0;
    long long backtracking_backtracks = 0;
    std::map<Strategy, int> strategy_usage;
    std::string hardest_strategy = "Brak";
    std::vector<std::string> debug_logic_logs;
};

struct FolderStats {
    fs::path relative_folder;
    long long non_empty_lines = 0;
    long long invalid_lines = 0;
    long long analyzed_puzzles = 0;
    long long contradictions = 0;
    long long solved_logically = 0;
    long long requires_guessing = 0;
    long long solved_with_backtracking = 0;
    long long unique_solutions = 0;
    long long multiple_solutions = 0;
    long long no_solution = 0;
    long long backtracking_nodes_sum = 0;
    long long backtracking_decisions_sum = 0;
    long long backtracking_backtracks_sum = 0;
    long long clues_sum = 0;
    long long difficulty_sum = 0;
    long long difficulty_count = 0;
    int max_difficulty = 0;
    int hardest_rank_seen = 0;
    std::string hardest_name_seen = "Brak";
    std::map<Strategy, long long> strategy_usage;
    std::map<std::string, long long> hardest_histogram;
    std::vector<PuzzleReportEntry> puzzle_reports;
};

static std::string trim(const std::string& s) {
    const std::string ws = " \t\r\n";
    const std::size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const std::size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

static bool parseIntStrict(const std::string& t, int& out) {
    try {
        std::size_t p = 0;
        const int v = std::stoi(t, &p);
        if (p != t.size()) return false;
        out = v;
        return true;
    } catch (...) { return false; }
}

static bool parseLLStrict(const std::string& t, long long& out) {
    try {
        std::size_t p = 0;
        const long long v = std::stoll(t, &p);
        if (p != t.size()) return false;
        out = v;
        return true;
    } catch (...) { return false; }
}

static int bits(unsigned int m) { int c = 0; while (m) { m &= (m - 1U); ++c; } return c; }
static int firstDigit(unsigned int m) {
    int d = 1;
    while (m) { if (m & 1U) return d; m >>= 1U; ++d; }
    return 0;
}
static std::vector<int> digitsFromMask(unsigned int m) {
    std::vector<int> out;
    int d = 1;
    while (m) {
        if (m & 1U) out.push_back(d);
        m >>= 1U;
        ++d;
    }
    return out;
}
static int maskMinDigit(unsigned int m) {
    int d = 1;
    while (m) {
        if (m & 1U) return d;
        m >>= 1U;
        ++d;
    }
    return 0;
}

static std::string strategyName(Strategy s) {
    switch (s) {
        case Strategy::NakedSingle: return "Naked Single";
        case Strategy::HiddenSingle: return "Hidden Single";
        case Strategy::NakedPair: return "Naked Pair";
        case Strategy::HiddenPair: return "Hidden Pair";
        case Strategy::PointingPairsTriples: return "Intersection Removal (Pointing Pairs/Triples)";
        case Strategy::BoxLineReduction: return "Intersection Removal (Box/Line Reduction)";
        case Strategy::NakedTriple: return "Naked Triple";
        case Strategy::HiddenTriple: return "Hidden Triple";
        case Strategy::NakedQuad: return "Naked Quad";
        case Strategy::HiddenQuad: return "Hidden Quad";
        case Strategy::XWing: return "X-Wing";
        case Strategy::YWing: return "Y-Wing (XY-Wing)";
        case Strategy::XYZWing: return "XYZ-Wing";
        case Strategy::WXYZWing: return "WXYZ-Wing";
        case Strategy::Swordfish: return "Swordfish";
        case Strategy::Jellyfish: return "Jellyfish";
        case Strategy::FrankenMutantFish: return "Franken / Mutant Fish";
        case Strategy::KrakenFish: return "Kraken Fish";
        case Strategy::Skyscraper: return "Skyscraper";
        case Strategy::TwoStringKite: return "2-String Kite";
        case Strategy::SimpleColoring: return "Simple Coloring";
        case Strategy::ThreeDMedusa: return "3D Medusa";
        case Strategy::FinnedXWingSashimi: return "Finned X-Wing / Sashimi X-Wing";
        case Strategy::FinnedSwordfish: return "Finned Swordfish";
        case Strategy::FinnedJellyfish: return "Finned Jellyfish";
        case Strategy::EmptyRectangle: return "Empty Rectangle";
        case Strategy::UniqueRectangle: return "Unique Rectangle";
        case Strategy::UniqueLoop: return "Unique Loop (6+)";
        case Strategy::BivalueOddagon: return "Bivalue Oddagon";
        case Strategy::AvoidableRectangle: return "Avoidable Rectangle";
        case Strategy::BUGPlus1: return "BUG+1";
        case Strategy::RemotePairs: return "Remote Pairs";
        case Strategy::WWing: return "W-Wing";
        case Strategy::GroupedXCycle: return "Grouped X-Cycle";
        case Strategy::XChain: return "X-Chain";
        case Strategy::XYChain: return "XY-Chain";
        case Strategy::GroupedAIC: return "Grouped AIC";
        case Strategy::AIC: return "AIC (Alternating Inference Chains)";
        case Strategy::ContinuousNiceLoop: return "Continuous Nice Loop";
        case Strategy::ALSXZ: return "ALS (Rule XZ)";
        case Strategy::ALSXYWing: return "ALS-XY-Wing";
        case Strategy::ALSChain: return "ALS-Chain";
        case Strategy::DeathBlossom: return "Death Blossom";
        case Strategy::SueDeCoq: return "Sue de Coq";
        case Strategy::MSLS: return "MSLS (Multi-Sector Locked Sets)";
        case Strategy::Exocet: return "Exocet";
        case Strategy::SeniorExocet: return "Senior Exocet";
        case Strategy::SKLoop: return "SK Loop";
        case Strategy::PatternOverlayMethod: return "Pattern Overlay Method";
        case Strategy::ForcingChains: return "Forcing Chains";
        case Strategy::Backtracking: return "Backtracking";
    }
    return "Unknown";
}

static int strategyRank(Strategy s) {
    switch (s) {
        case Strategy::NakedSingle: return 1;
        case Strategy::HiddenSingle: return 2;
        case Strategy::NakedPair: return 2;
        case Strategy::HiddenPair: return 2;
        case Strategy::PointingPairsTriples: return 3;
        case Strategy::BoxLineReduction: return 3;
        case Strategy::NakedTriple: return 3;
        case Strategy::HiddenTriple: return 4;
        case Strategy::NakedQuad: return 4;
        case Strategy::HiddenQuad: return 5;
        case Strategy::XWing: return 6;
        case Strategy::YWing: return 6;
        case Strategy::XYZWing: return 7;
        case Strategy::WXYZWing: return 8;
        case Strategy::Swordfish: return 7;
        case Strategy::KrakenFish: return 11;
        case Strategy::Skyscraper: return 7;
        case Strategy::TwoStringKite: return 7;
        case Strategy::SimpleColoring: return 7;
        case Strategy::ThreeDMedusa: return 7;
        case Strategy::Jellyfish: return 8;
        case Strategy::FrankenMutantFish: return 10;
        case Strategy::FinnedXWingSashimi: return 8;
        case Strategy::FinnedSwordfish: return 9;
        case Strategy::FinnedJellyfish: return 10;
        case Strategy::EmptyRectangle: return 8;
        case Strategy::UniqueRectangle: return 8;
        case Strategy::UniqueLoop: return 9;
        case Strategy::BivalueOddagon: return 10;
        case Strategy::AvoidableRectangle: return 8;
        case Strategy::BUGPlus1: return 8;
        case Strategy::RemotePairs: return 8;
        case Strategy::WWing: return 9;
        case Strategy::GroupedXCycle: return 9;
        case Strategy::XChain: return 9;
        case Strategy::XYChain: return 9;
        case Strategy::GroupedAIC: return 10;
        case Strategy::AIC: return 10;
        case Strategy::ContinuousNiceLoop: return 10;
        case Strategy::ALSXZ: return 10;
        case Strategy::ALSXYWing: return 10;
        case Strategy::ALSChain: return 11;
        case Strategy::DeathBlossom: return 11;
        case Strategy::SueDeCoq: return 10;
        case Strategy::MSLS: return 11;
        case Strategy::Exocet: return 11;
        case Strategy::SeniorExocet: return 11;
        case Strategy::SKLoop: return 11;
        case Strategy::PatternOverlayMethod: return 11;
        case Strategy::ForcingChains: return 11;
        case Strategy::Backtracking: return 11;
    }
    return 0;
}

static std::string strategyImplementationStatus(Strategy s) {
    switch (s) {
        case Strategy::NakedSingle:
        case Strategy::HiddenSingle:
        case Strategy::NakedPair:
        case Strategy::HiddenPair:
        case Strategy::PointingPairsTriples:
        case Strategy::BoxLineReduction:
        case Strategy::NakedTriple:
        case Strategy::HiddenTriple:
        case Strategy::NakedQuad:
        case Strategy::HiddenQuad:
        case Strategy::XWing:
        case Strategy::YWing:
        case Strategy::XYZWing:
        case Strategy::WXYZWing:
        case Strategy::Swordfish:
        case Strategy::Jellyfish:
        case Strategy::FrankenMutantFish:
        case Strategy::KrakenFish:
        case Strategy::Skyscraper:
        case Strategy::TwoStringKite:
        case Strategy::SimpleColoring:
        case Strategy::ThreeDMedusa:
        case Strategy::EmptyRectangle:
        case Strategy::BUGPlus1:
        case Strategy::RemotePairs:
        case Strategy::WWing:
        case Strategy::AvoidableRectangle:
        case Strategy::GroupedXCycle:
        case Strategy::XChain:
        case Strategy::XYChain:
        case Strategy::AIC:
        case Strategy::ContinuousNiceLoop:
        case Strategy::ALSXZ:
        case Strategy::ALSXYWing:
        case Strategy::ALSChain:
        case Strategy::DeathBlossom:
        case Strategy::SueDeCoq:
        case Strategy::MSLS:
        case Strategy::SeniorExocet:
        case Strategy::SKLoop:
        case Strategy::PatternOverlayMethod:
        case Strategy::Backtracking:
            return "ZAIMPLEMENTOWANE";
        case Strategy::UniqueRectangle:
            return "ZAIMPLEMENTOWANE (Type 1/2/3/4/5/6 + Hidden UR, konserwatywne)";
        case Strategy::UniqueLoop:
        case Strategy::BivalueOddagon:
            return "ZAIMPLEMENTOWANE (konserwatywne, single-extra)";
        case Strategy::FinnedXWingSashimi:
            return "CZESCIOWO (finned core)";
        case Strategy::FinnedSwordfish:
        case Strategy::FinnedJellyfish:
            return "CZESCIOWO (finned fish, konserwatywne)";
        case Strategy::GroupedAIC:
            return "ZAIMPLEMENTOWANE (implikacyjne, bez backtrackingu)";
        case Strategy::Exocet:
            return "CZESCIOWO (target-check, bez backtrackingu)";
        case Strategy::ForcingChains:
            return "ZAIMPLEMENTOWANE (implikacyjne, bez backtrackingu)";
        default:
            return "NIEZAIMPLEMENTOWANE";
    }
}

class SudokuAnalyzer {
public:
    explicit SudokuAnalyzer(const SudokuBoard& b);
    AnalysisReport run();

private:
    const SudokuBoard& b_;
    int N_ = 0, BR_ = 0, BC_ = 0, NN_ = 0;
    unsigned int all_ = 0U;
    bool contradiction_ = false;
    int hardest_rank_ = 0;
    std::string hardest_name_;
    std::map<Strategy, int> usage_;
    bool debug_logic_enabled_ = true;
    std::size_t debug_logic_limit_ = 300;
    bool debug_logic_truncated_ = false;
    std::vector<std::string> debug_logic_logs_;
    std::vector<int> grid_;
    std::vector<unsigned int> cand_;
    std::vector<std::vector<int>> houses_;
    std::vector<std::array<int, 3>> cell_houses_;
    std::vector<std::vector<int>> peers_;

    unsigned int bit(int d) const { return 1U << (d - 1); }
    int row(int i) const { return i / N_; }
    int col(int i) const { return i % N_; }
    int box(int r, int c) const {
        const int bpr = N_ / BC_;
        return (r / BR_) * bpr + (c / BC_);
    }
    bool solved() const;
    void buildTopology();
    bool removeCandidate(int i, int d, bool& changed);
    bool assignValue(int i, int d);
    void initCandidates();
    void use(Strategy s, int amount);
    bool applyNakedSingles(int& n);
    bool applyHiddenSingles(int& n);
    void combosRec(const std::vector<int>& src, int k, int start, std::vector<int>& cur,
                   const std::function<void(const std::vector<int>&)>& cb);
    void forEachCombo(const std::vector<int>& src, int k, const std::function<void(const std::vector<int>&)>& cb);
    bool applyNakedSubset(int k, int& n);
    bool applyHiddenSubset(int k, int& n);
    bool applyPointingPairsTriples(int& n);
    bool applyBoxLineReduction(int& n);
    bool applyFish(int size, int& n);
    bool applyKrakenFish(int& n);
    bool applyYWing(int& n);
    bool applyXYZWing(int& n);
    bool applyWXYZWing(int& n);
    bool applySkyscraper(int& n);
    bool applyTwoStringKite(int& n);
    bool applySimpleColoring(int& n);
    bool applyFinnedXWingSashimi(int& n);
    bool applyFinnedFish(int size, int& n);
    bool applyFrankenMutantFish(int size, int& n);
    bool applyEmptyRectangle(int& n);
    bool applyUniqueRectangleType1(int& n);
    bool applyUniqueRectangleType2to6(int& n);
    bool applyUniqueLoop(int& n);
    bool applyBivalueOddagon(int& n);
    bool applyAvoidableRectangle(int& n);
    bool applyBUGPlus1(int& n);
    bool applyRemotePairs(int& n);
    bool applyWWing(int& n);
    bool applyXChain(int& n);
    bool applyXYChain(int& n);
    bool applyAIC(int& n);
    bool applyContinuousNiceLoop(int& n);
    bool applyALSXZ(int& n);
    bool applyALSXYWing(int& n);
    bool applyALSChain(int& n);
    bool applyDeathBlossom(int& n);
    bool applySueDeCoq(int& n);
    bool applyMSLS(int& n);
    bool applyExocet(int& n);
    bool applySeniorExocet(int& n);
    bool applySKLoop(int& n);
    bool applyPatternOverlayMethod(int& n);
    bool applyForcingChains(int& n);
    bool applyGroupedXCycle(int& n);
    bool applyGroupedAIC(int& n);
    bool applyThreeDMedusa(int& n);
    void logicalSolve();
    bool isPeerCell(int a, int b) const;
    std::string cellName(int idx) const;
    void pushDebugLog(const std::string& line);
    bool hasLogicalSupportWithAssignments(const std::vector<std::pair<int, int>>& assignments) const;
    int cluesCount() const;
};

class BacktrackingCounter {
public:
    BacktrackingCounter(int br, int bc, int n, std::vector<int> grid);
    int countSolutions(int limit);

private:
    int BR_, BC_, N_, NN_, limit_ = 2, solutions_ = 0;
    unsigned int all_ = 0U;
    std::vector<int> grid_;
    unsigned int bit(int d) const { return 1U << (d - 1); }
    int row(int i) const { return i / N_; }
    int col(int i) const { return i % N_; }
    int box(int r, int c) const {
        const int bpr = N_ / BC_;
        return (r / BR_) * bpr + (c / BC_);
    }
    unsigned int allowed(int idx) const;
    bool validState() const;
    void search();
};

struct BacktrackingSolveStats {
    bool solved = false;
    long long nodes = 0;
    long long decisions = 0;
    long long backtracks = 0;
};

class BacktrackingSolver {
public:
    BacktrackingSolver(int br, int bc, int n, std::vector<int> grid);
    BacktrackingSolveStats solve();

private:
    int BR_, BC_, N_, NN_;
    unsigned int all_ = 0U;
    std::vector<int> grid_;
    BacktrackingSolveStats stats_;
    unsigned int bit(int d) const { return 1U << (d - 1); }
    int row(int i) const { return i / N_; }
    int col(int i) const { return i % N_; }
    int box(int r, int c) const {
        const int bpr = N_ / BC_;
        return (r / BR_) * bpr + (c / BC_);
    }
    unsigned int allowed(int idx) const;
    bool validState() const;
    bool search();
};

static SudokuBoard parseSudokuLine(const std::string& line);
static int countSolutionsWithBacktracking(const SudokuBoard& b, int limit = 2);
static BacktrackingSolveStats solveWithBacktracking(const SudokuBoard& b, const std::vector<int>& initialGrid);
static std::string selectFolderModern();
static bool isTxtFile(const fs::path& p);
static bool isPathWithin(const fs::path& path, const fs::path& parent);
static std::vector<fs::path> collectTxtFilesRecursive(const fs::path& root, const fs::path& excludedRoot);
static long long countNonEmptyLinesInTxtFiles(const std::vector<fs::path>& files);
static std::string folderKeyFromRelativePath(const fs::path& rel);
static std::string sanitizeFileName(const std::string& name);
static std::string csvEscape(const std::string& field);
static int difficultyLevelFromReport(const AnalysisReport& report);
static std::string difficultyTypeFromReport(const AnalysisReport& report);
static std::string boardTypeFromBoard(const SudokuBoard& board);
static void appendInvalidPuzzleReport(FolderStats& stats, const std::string& sourceFile, int lineNo,
                                      const std::string& parseError);
static void appendValidPuzzleReport(FolderStats& stats, const std::string& sourceFile, int lineNo,
                                    const SudokuBoard& board, const AnalysisReport& report);
static void updateFolderStats(FolderStats& stats, const AnalysisReport& report);
static void writeFolderReport(const fs::path& outDir, const std::string& folderKey, const FolderStats& stats);
static void writeGlobalSummary(const fs::path& outDir, const std::map<std::string, FolderStats>& allStats,
                               long long txtFilesScanned);
static void writeFolderCsv(const fs::path& outDir, const std::map<std::string, FolderStats>& allStats);

SudokuAnalyzer::SudokuAnalyzer(const SudokuBoard& b) : b_(b) {
    N_ = b_.side_size;
    BR_ = b_.block_rows;
    BC_ = b_.block_cols;
    NN_ = N_ * N_;
    all_ = (N_ >= 31) ? 0U : ((1U << N_) - 1U);
    buildTopology();
    initCandidates();
}

int SudokuAnalyzer::cluesCount() const {
    int c = 0;
    for (const Cell& x : b_.cells) if (x.revealed) ++c;
    return c;
}

bool SudokuAnalyzer::solved() const {
    if (contradiction_) return false;
    for (int v : grid_) if (v == 0) return false;
    return true;
}

bool SudokuAnalyzer::isPeerCell(int a, int b) const {
    if (a == b) return false;
    const int ra = row(a), ca = col(a);
    const int rb = row(b), cb = col(b);
    if (ra == rb || ca == cb) return true;
    return box(ra, ca) == box(rb, cb);
}

std::string SudokuAnalyzer::cellName(int idx) const {
    std::ostringstream ss;
    ss << "r" << (row(idx) + 1) << "c" << (col(idx) + 1);
    return ss.str();
}

void SudokuAnalyzer::pushDebugLog(const std::string& line) {
    if (!debug_logic_enabled_) return;
    if (debug_logic_logs_.size() >= debug_logic_limit_) {
        if (!debug_logic_truncated_) {
            debug_logic_logs_.push_back("... log obciety (osiagnieto limit wpisow) ...");
            debug_logic_truncated_ = true;
        }
        return;
    }
    debug_logic_logs_.push_back(line);
}

bool SudokuAnalyzer::hasLogicalSupportWithAssignments(const std::vector<std::pair<int, int>>& assignments) const {
    std::vector<int> g = grid_;
    std::vector<unsigned int> c = cand_;
    std::vector<int> queue;
    queue.reserve(NN_);

    auto assignLocal = [&](int cell, int digit) -> bool {
        if (cell < 0 || cell >= NN_ || digit < 1 || digit > N_) return false;
        const unsigned int b = bit(digit);
        if (g[cell] == digit) return true;
        if (g[cell] != 0 && g[cell] != digit) return false;
        if ((c[cell] & b) == 0U) return false;
        g[cell] = digit;
        c[cell] = b;
        queue.push_back(cell);
        return true;
    };

    auto removeLocal = [&](int cell, int digit) -> bool {
        if (cell < 0 || cell >= NN_ || digit < 1 || digit > N_) return false;
        if (g[cell] != 0) return g[cell] != digit;
        const unsigned int b = bit(digit);
        if ((c[cell] & b) == 0U) return true;
        c[cell] &= ~b;
        if (c[cell] == 0U) return false;
        if (bits(c[cell]) == 1) {
            if (!assignLocal(cell, firstDigit(c[cell]))) return false;
        }
        return true;
    };

    for (const auto& asg : assignments) {
        if (!assignLocal(asg.first, asg.second)) return false;
    }

    std::size_t qHead = 0;
    while (true) {
        while (qHead < queue.size()) {
            const int cell = queue[qHead++];
            const int d = g[cell];
            if (d <= 0) return false;
            for (int p : peers_[cell]) {
                if (!removeLocal(p, d)) return false;
            }
        }

        bool pushed = false;
        for (int h = 0; h < 3 * N_; ++h) {
            for (int d = 1; d <= N_; ++d) {
                int solvedCnt = 0;
                int lastPos = -1;
                int places = 0;
                for (int idx : houses_[h]) {
                    if (g[idx] == d) {
                        ++solvedCnt;
                        if (solvedCnt > 1) return false;
                        continue;
                    }
                    if (g[idx] == 0 && (c[idx] & bit(d))) {
                        ++places;
                        lastPos = idx;
                    }
                }
                if (solvedCnt == 0 && places == 0) return false;
                if (solvedCnt == 0 && places == 1) {
                    const std::size_t before = queue.size();
                    if (!assignLocal(lastPos, d)) return false;
                    if (queue.size() > before) pushed = true;
                }
            }
        }

        if (qHead < queue.size()) continue;
        if (!pushed) break;
    }

    return true;
}

void SudokuAnalyzer::buildTopology() {
    houses_.assign(3 * N_, {});
    cell_houses_.assign(NN_, {0, 0, 0});
    peers_.assign(NN_, {});

    for (int r = 0; r < N_; ++r) {
        for (int c = 0; c < N_; ++c) {
            const int idx = r * N_ + c;
            houses_[r].push_back(idx);
            houses_[N_ + c].push_back(idx);
            houses_[2 * N_ + box(r, c)].push_back(idx);
        }
    }
    for (int i = 0; i < NN_; ++i) {
        const int r = row(i), c = col(i), b = box(r, c);
        cell_houses_[i] = {r, N_ + c, 2 * N_ + b};
    }
    for (int i = 0; i < NN_; ++i) {
        std::vector<char> seen(NN_, 0);
        for (int h : cell_houses_[i]) {
            for (int p : houses_[h]) {
                if (p == i || seen[p]) continue;
                seen[p] = 1;
                peers_[i].push_back(p);
            }
        }
    }
}

bool SudokuAnalyzer::removeCandidate(int i, int d, bool& changed) {
    if (grid_[i] != 0) {
        if (grid_[i] == d) {
            contradiction_ = true;
            return false;
        }
        return true;
    }
    const unsigned int b = bit(d);
    if ((cand_[i] & b) == 0U) return true;
    const unsigned int next = cand_[i] & (~b);
    if (next == 0U) {
        contradiction_ = true;
        return false;
    }
    cand_[i] = next;
    changed = true;
    return true;
}

bool SudokuAnalyzer::assignValue(int i, int d) {
    if (d < 1 || d > N_) {
        contradiction_ = true;
        return false;
    }
    if (grid_[i] == d) return true;
    if (grid_[i] != 0 && grid_[i] != d) {
        contradiction_ = true;
        return false;
    }
    if ((cand_[i] & bit(d)) == 0U) {
        contradiction_ = true;
        return false;
    }

    grid_[i] = d;
    cand_[i] = bit(d);
    for (int p : peers_[i]) {
        if (grid_[p] == d) {
            contradiction_ = true;
            return false;
        }
        bool changed = false;
        if (!removeCandidate(p, d, changed)) return false;
    }
    return true;
}

void SudokuAnalyzer::initCandidates() {
    grid_.assign(NN_, 0);
    cand_.assign(NN_, all_);
    for (int i = 0; i < NN_; ++i) {
        if (!b_.cells[i].revealed) continue;
        if (!assignValue(i, b_.cells[i].value)) {
            contradiction_ = true;
            return;
        }
    }
}

void SudokuAnalyzer::use(Strategy s, int amount) {
    if (amount <= 0) return;
    usage_[s] += amount;
    const int r = strategyRank(s);
    if (r > hardest_rank_) {
        hardest_rank_ = r;
        hardest_name_ = strategyName(s);
    }
}

bool SudokuAnalyzer::applyNakedSingles(int& n) {
    n = 0;
    bool changed = false;
    for (int i = 0; i < NN_; ++i) {
        if (grid_[i] != 0) continue;
        if (bits(cand_[i]) != 1) continue;
        if (!assignValue(i, firstDigit(cand_[i]))) return changed;
        changed = true;
        ++n;
    }
    return changed;
}

bool SudokuAnalyzer::applyHiddenSingles(int& n) {
    n = 0;
    for (const auto& h : houses_) {
        for (int d = 1; d <= N_; ++d) {
            bool already = false;
            std::vector<int> where;
            for (int i : h) {
                if (grid_[i] == d) { already = true; break; }
                if (grid_[i] == 0 && (cand_[i] & bit(d))) where.push_back(i);
            }
            if (already) continue;
            if (where.empty()) {
                contradiction_ = true;
                return false;
            }
            if (where.size() == 1U) {
                if (!assignValue(where[0], d)) return false;
                n = 1;
                return true;
            }
        }
    }
    return false;
}

void SudokuAnalyzer::combosRec(const std::vector<int>& src, int k, int start, std::vector<int>& cur,
                               const std::function<void(const std::vector<int>&)>& cb) {
    if (static_cast<int>(cur.size()) == k) {
        cb(cur);
        return;
    }
    const int need = k - static_cast<int>(cur.size());
    const int lim = static_cast<int>(src.size()) - need;
    for (int i = start; i <= lim; ++i) {
        cur.push_back(src[i]);
        combosRec(src, k, i + 1, cur, cb);
        cur.pop_back();
    }
}

void SudokuAnalyzer::forEachCombo(const std::vector<int>& src, int k,
                                  const std::function<void(const std::vector<int>&)>& cb) {
    std::vector<int> cur;
    combosRec(src, k, 0, cur, cb);
}

bool SudokuAnalyzer::applyNakedSubset(int k, int& n) {
    n = 0;
    for (const auto& h : houses_) {
        std::vector<int> pool;
        for (int i : h) {
            if (grid_[i] != 0) continue;
            const int bc = bits(cand_[i]);
            if (bc >= 2 && bc <= k) pool.push_back(i);
        }
        if (static_cast<int>(pool.size()) < k) continue;

        bool found = false;
        // Naked subset: k komorek tworzy lacznie dokladnie k cyfr,
        // wiec te cyfry mozna usunac z pozostalych komorek domu.
        forEachCombo(pool, k, [&](const std::vector<int>& cset) {
            if (found || contradiction_) return;
            unsigned int uni = 0U;
            for (int i : cset) uni |= cand_[i];
            if (bits(uni) != k) return;

            std::vector<int> subset;
            for (int i : h) if (grid_[i] == 0 && (cand_[i] & ~uni) == 0U) subset.push_back(i);
            if (static_cast<int>(subset.size()) != k) return;

            std::set<int> subset_set(subset.begin(), subset.end());
            int local = 0;
            for (int i : h) {
                if (grid_[i] != 0 || subset_set.count(i)) continue;
                unsigned int m = uni;
                while (m) {
                    const unsigned int one = m & (~m + 1U);
                    bool changed = false;
                    if (!removeCandidate(i, firstDigit(one), changed)) return;
                    if (changed) ++local;
                    m &= (m - 1U);
                }
            }
            if (local > 0) {
                n = local;
                found = true;
            }
        });
        if (found) return true;
    }
    return false;
}

bool SudokuAnalyzer::applyHiddenSubset(int k, int& n) {
    n = 0;
    std::vector<int> digits;
    for (int d = 1; d <= N_; ++d) digits.push_back(d);

    for (const auto& h : houses_) {
        bool found = false;
        // Hidden subset: k cyfr wystepuje tylko w k komorkach,
        // wiec w tych komorkach usuwamy wszystkich innych kandydatow.
        forEachCombo(digits, k, [&](const std::vector<int>& dset) {
            if (found || contradiction_) return;
            unsigned int dm = 0U;
            for (int d : dset) dm |= bit(d);

            std::vector<int> union_cells;
            for (int i : h) if (grid_[i] == 0 && (cand_[i] & dm)) union_cells.push_back(i);
            if (static_cast<int>(union_cells.size()) != k) return;

            for (int d : dset) {
                bool present = false;
                for (int i : union_cells) if (cand_[i] & bit(d)) { present = true; break; }
                if (!present) return;
            }

            int local = 0;
            for (int i : union_cells) {
                unsigned int extra = cand_[i] & (~dm);
                while (extra) {
                    const unsigned int one = extra & (~extra + 1U);
                    bool changed = false;
                    if (!removeCandidate(i, firstDigit(one), changed)) return;
                    if (changed) ++local;
                    extra &= (extra - 1U);
                }
            }
            if (local > 0) {
                n = local;
                found = true;
            }
        });
        if (found) return true;
    }
    return false;
}

bool SudokuAnalyzer::applyPointingPairsTriples(int& n) {
    n = 0;
    for (int bh = 2 * N_; bh < 3 * N_; ++bh) {
        const auto& box_cells = houses_[bh];
        for (int d = 1; d <= N_; ++d) {
            std::vector<int> w;
            for (int i : box_cells) if (grid_[i] == 0 && (cand_[i] & bit(d))) w.push_back(i);
            if (w.size() < 2U) continue;

            bool same_row = true;
            const int r0 = row(w[0]);
            for (int i : w) if (row(i) != r0) { same_row = false; break; }
            if (same_row) {
                int local = 0;
                for (int i : houses_[r0]) {
                    if (grid_[i] != 0 || box(row(i), col(i)) == (bh - 2 * N_)) continue;
                    bool changed = false;
                    if (!removeCandidate(i, d, changed)) return false;
                    if (changed) ++local;
                }
                if (local > 0) { n = local; return true; }
            }

            bool same_col = true;
            const int c0 = col(w[0]);
            for (int i : w) if (col(i) != c0) { same_col = false; break; }
            if (same_col) {
                int local = 0;
                for (int i : houses_[N_ + c0]) {
                    if (grid_[i] != 0 || box(row(i), col(i)) == (bh - 2 * N_)) continue;
                    bool changed = false;
                    if (!removeCandidate(i, d, changed)) return false;
                    if (changed) ++local;
                }
                if (local > 0) { n = local; return true; }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyBoxLineReduction(int& n) {
    n = 0;
    for (int r = 0; r < N_; ++r) {
        for (int d = 1; d <= N_; ++d) {
            std::vector<int> w;
            for (int i : houses_[r]) if (grid_[i] == 0 && (cand_[i] & bit(d))) w.push_back(i);
            if (w.size() < 2U) continue;
            const int b0 = box(row(w[0]), col(w[0]));
            bool same_box = true;
            for (int i : w) if (box(row(i), col(i)) != b0) { same_box = false; break; }
            if (!same_box) continue;

            int local = 0;
            for (int i : houses_[2 * N_ + b0]) {
                if (grid_[i] != 0 || row(i) == r) continue;
                bool changed = false;
                if (!removeCandidate(i, d, changed)) return false;
                if (changed) ++local;
            }
            if (local > 0) { n = local; return true; }
        }
    }

    for (int c = 0; c < N_; ++c) {
        for (int d = 1; d <= N_; ++d) {
            std::vector<int> w;
            for (int i : houses_[N_ + c]) if (grid_[i] == 0 && (cand_[i] & bit(d))) w.push_back(i);
            if (w.size() < 2U) continue;
            const int b0 = box(row(w[0]), col(w[0]));
            bool same_box = true;
            for (int i : w) if (box(row(i), col(i)) != b0) { same_box = false; break; }
            if (!same_box) continue;

            int local = 0;
            for (int i : houses_[2 * N_ + b0]) {
                if (grid_[i] != 0 || col(i) == c) continue;
                bool changed = false;
                if (!removeCandidate(i, d, changed)) return false;
                if (changed) ++local;
            }
            if (local > 0) { n = local; return true; }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyFish(int size, int& n) {
    n = 0;
    if (size < 2 || size > N_) return false;

    for (int d = 1; d <= N_; ++d) {
        for (int mode = 0; mode < 2; ++mode) {
            std::vector<std::vector<int>> positions(N_);
            std::vector<int> eligible;
            for (int line = 0; line < N_; ++line) {
                std::vector<int> pos;
                for (int p = 0; p < N_; ++p) {
                    const int idx = (mode == 0) ? (line * N_ + p) : (p * N_ + line);
                    if (grid_[idx] == 0 && (cand_[idx] & bit(d))) {
                        pos.push_back(p);
                    }
                }
                positions[line] = std::move(pos);
                const int cnt = static_cast<int>(positions[line].size());
                if (cnt >= 2 && cnt <= size) {
                    eligible.push_back(line);
                }
            }
            if (static_cast<int>(eligible.size()) < size) continue;

            bool found = false;
            forEachCombo(eligible, size, [&](const std::vector<int>& lines) {
                if (found || contradiction_) return;

                std::vector<char> inLine(N_, 0), inPos(N_, 0);
                for (int line : lines) {
                    inLine[line] = 1;
                    for (int p : positions[line]) inPos[p] = 1;
                }
                int unionCount = 0;
                for (char f : inPos) if (f) ++unionCount;
                if (unionCount != size) return;

                int local = 0;
                for (int p = 0; p < N_; ++p) {
                    if (!inPos[p]) continue;
                    for (int line = 0; line < N_; ++line) {
                        if (inLine[line]) continue;
                        const int idx = (mode == 0) ? (line * N_ + p) : (p * N_ + line);
                        bool changed = false;
                        if (!removeCandidate(idx, d, changed)) return;
                        if (changed) ++local;
                    }
                }
                if (local > 0) {
                    n = local;
                    found = true;
                }
            });
            if (found) return true;
        }
    }
    return false;
}

bool SudokuAnalyzer::applyKrakenFish(int& n) {
    n = 0;
    if (N_ != 9) return false;

    auto buildStrongAdj = [&](int d) {
        std::vector<std::vector<int>> adj(NN_);
        for (int h = 0; h < 3 * N_; ++h) {
            std::vector<int> where;
            for (int idx : houses_[h]) {
                if (grid_[idx] == 0 && (cand_[idx] & bit(d))) where.push_back(idx);
            }
            if (where.size() == 2U) {
                adj[where[0]].push_back(where[1]);
                adj[where[1]].push_back(where[0]);
            }
        }
        return adj;
    };

    for (int d = 1; d <= N_; ++d) {
        const std::vector<std::vector<int>> strongAdj = buildStrongAdj(d);
        auto victimSeesFinByChain = [&](int victim, int fin) {
            if (isPeerCell(victim, fin)) return true;
            for (int s1 : strongAdj[fin]) {
                if (isPeerCell(victim, s1)) return true;
                for (int s2 : strongAdj[s1]) {
                    if (s2 == fin) continue;
                    if (isPeerCell(victim, s2)) return true;
                }
            }
            return false;
        };

        auto process = [&](bool rowBased, int fishSize) -> bool {
            std::vector<std::vector<int>> positions(N_);
            std::vector<int> eligible;
            for (int line = 0; line < N_; ++line) {
                std::vector<int> pos;
                for (int p = 0; p < N_; ++p) {
                    const int idx = rowBased ? (line * N_ + p) : (p * N_ + line);
                    if (grid_[idx] == 0 && (cand_[idx] & bit(d))) pos.push_back(p);
                }
                positions[line] = std::move(pos);
                const int cnt = static_cast<int>(positions[line].size());
                if (cnt >= 2 && cnt <= fishSize + 2) eligible.push_back(line);
            }
            if (static_cast<int>(eligible.size()) < fishSize) return false;

            bool found = false;
            forEachCombo(eligible, fishSize, [&](const std::vector<int>& lines) {
                if (found || contradiction_) return;

                std::vector<char> inLine(N_, 0), inUnionPos(N_, 0);
                std::vector<int> unionPos;
                for (int line : lines) {
                    inLine[line] = 1;
                    for (int p : positions[line]) {
                        if (!inUnionPos[p]) {
                            inUnionPos[p] = 1;
                            unionPos.push_back(p);
                        }
                    }
                }
                const int uc = static_cast<int>(unionPos.size());
                if (uc <= fishSize || uc > fishSize + 2) return;

                forEachCombo(unionPos, fishSize, [&](const std::vector<int>& coverPos) {
                    if (found || contradiction_) return;
                    std::vector<char> inCover(N_, 0);
                    for (int p : coverPos) inCover[p] = 1;

                    std::vector<int> finCells;
                    for (int line : lines) {
                        for (int p : positions[line]) {
                            if (inCover[p]) continue;
                            finCells.push_back(rowBased ? (line * N_ + p) : (p * N_ + line));
                        }
                    }
                    if (finCells.empty()) return;

                    int finBox = -1;
                    bool sameBox = true;
                    for (int f : finCells) {
                        const int b = box(row(f), col(f));
                        if (finBox < 0) finBox = b;
                        else if (finBox != b) { sameBox = false; break; }
                    }
                    if (!sameBox || finBox < 0) return;

                    int local = 0;
                    for (int p : coverPos) {
                        const std::vector<int>& coverHouse = rowBased ? houses_[N_ + p] : houses_[p];
                        for (int idx : coverHouse) {
                            const int line = rowBased ? row(idx) : col(idx);
                            if (inLine[line]) continue;
                            if (grid_[idx] != 0 || (cand_[idx] & bit(d)) == 0U) continue;
                            if (box(row(idx), col(idx)) != finBox) continue;

                            bool seesAll = true;
                            for (int f : finCells) {
                                if (!victimSeesFinByChain(idx, f)) { seesAll = false; break; }
                            }
                            if (!seesAll) continue;

                            bool changed = false;
                            if (!removeCandidate(idx, d, changed)) return;
                            if (changed) ++local;
                        }
                    }
                    if (local > 0) {
                        std::ostringstream ss;
                        ss << "KrakenFish(" << fishSize << "): remove " << d << " from " << local
                           << " cell(s) via fin-chain support";
                        pushDebugLog(ss.str());
                        n = local;
                        found = true;
                    }
                });
            });
            return found;
        };

        if (process(true, 2) || process(false, 2) || process(true, 3) || process(false, 3)) return true;
    }
    return false;
}

bool SudokuAnalyzer::applyFrankenMutantFish(int size, int& n) {
    n = 0;
    if (N_ != 9 || size < 2 || size > 3) return false;

    auto tryMode = [&](bool rowBase) -> bool {
        std::vector<int> basePool;
        std::vector<int> coverPool;
        for (int r = 0; r < N_; ++r) basePool.push_back(r);
        for (int c = 0; c < N_; ++c) coverPool.push_back(N_ + c);
        for (int b = 0; b < N_; ++b) {
            const int boxHouse = 2 * N_ + b;
            basePool.push_back(boxHouse);
            coverPool.push_back(boxHouse);
        }
        if (!rowBase) std::swap(basePool, coverPool);

        for (int d = 1; d <= N_; ++d) {
            std::vector<int> baseEligible;
            std::vector<int> coverEligible;
            for (int h : basePool) {
                int cnt = 0;
                for (int idx : houses_[h]) if (grid_[idx] == 0 && (cand_[idx] & bit(d))) ++cnt;
                if (cnt >= 2 && cnt <= 6) baseEligible.push_back(h);
            }
            for (int h : coverPool) {
                int cnt = 0;
                for (int idx : houses_[h]) if (grid_[idx] == 0 && (cand_[idx] & bit(d))) ++cnt;
                if (cnt >= 2 && cnt <= 6) coverEligible.push_back(h);
            }
            if (static_cast<int>(baseEligible.size()) < size || static_cast<int>(coverEligible.size()) < size) continue;

            bool found = false;
            forEachCombo(baseEligible, size, [&](const std::vector<int>& baseSets) {
                if (found || contradiction_) return;
                std::vector<char> inBase(NN_, 0);
                std::vector<int> baseCandCells;
                for (int h : baseSets) {
                    for (int idx : houses_[h]) {
                        if (grid_[idx] != 0 || (cand_[idx] & bit(d)) == 0U) continue;
                        if (!inBase[idx]) {
                            inBase[idx] = 1;
                            baseCandCells.push_back(idx);
                        }
                    }
                }
                if (baseCandCells.size() < static_cast<std::size_t>(size)) return;

                forEachCombo(coverEligible, size, [&](const std::vector<int>& coverSets) {
                    if (found || contradiction_) return;
                    std::vector<char> inCover(NN_, 0);
                    for (int h : coverSets) {
                        for (int idx : houses_[h]) {
                            if (grid_[idx] == 0 && (cand_[idx] & bit(d))) inCover[idx] = 1;
                        }
                    }

                    for (int idx : baseCandCells) {
                        if (!inCover[idx]) return;
                    }

                    int local = 0;
                    for (int idx = 0; idx < NN_; ++idx) {
                        if (!inCover[idx] || inBase[idx]) continue;
                        if (grid_[idx] != 0 || (cand_[idx] & bit(d)) == 0U) continue;
                        bool changed = false;
                        if (!removeCandidate(idx, d, changed)) return;
                        if (changed) ++local;
                    }
                    if (local > 0) {
                        std::ostringstream ss;
                        ss << "FrankenMutantFish(" << size << "): remove " << d
                           << " from " << local << " cell(s)";
                        pushDebugLog(ss.str());
                        n = local;
                        found = true;
                    }
                });
            });
            if (found) return true;
        }
        return false;
    };

    if (tryMode(true)) return true;
    return tryMode(false);
}

bool SudokuAnalyzer::applyYWing(int& n) {
    n = 0;
    std::vector<int> bivalueCells;
    for (int i = 0; i < NN_; ++i) {
        if (grid_[i] == 0 && bits(cand_[i]) == 2) bivalueCells.push_back(i);
    }

    for (int pivot : bivalueCells) {
        const unsigned int pivotMask = cand_[pivot];
        const std::vector<int> pd = digitsFromMask(pivotMask);
        if (pd.size() != 2U) continue;
        const int a = pd[0], b = pd[1];

        for (int wing1 : peers_[pivot]) {
            if (grid_[wing1] != 0 || bits(cand_[wing1]) != 2) continue;
            const unsigned int m1 = cand_[wing1];
            const unsigned int common1 = m1 & pivotMask;
            if (bits(common1) != 1) continue;
            const int shared1 = firstDigit(common1);
            const int z1 = firstDigit(m1 & ~common1);
            if (z1 == 0) continue;

            const int neededShared = (shared1 == a) ? b : ((shared1 == b) ? a : 0);
            if (neededShared == 0) continue;

            for (int wing2 : peers_[pivot]) {
                if (wing2 == wing1) continue;
                if (grid_[wing2] != 0 || bits(cand_[wing2]) != 2) continue;
                const unsigned int m2 = cand_[wing2];
                const unsigned int common2 = m2 & pivotMask;
                if (bits(common2) != 1) continue;
                const int shared2 = firstDigit(common2);
                const int z2 = firstDigit(m2 & ~common2);
                if (shared2 != neededShared || z2 != z1) continue;

                int local = 0;
                for (int i = 0; i < NN_; ++i) {
                    if (i == pivot || i == wing1 || i == wing2) continue;
                    if (grid_[i] != 0 || (cand_[i] & bit(z1)) == 0U) continue;
                    if (!isPeerCell(i, wing1) || !isPeerCell(i, wing2)) continue;
                    bool changed = false;
                    if (!removeCandidate(i, z1, changed)) return false;
                    if (changed) ++local;
                }
                if (local > 0) { n = local; return true; }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyXYZWing(int& n) {
    n = 0;
    for (int pivot = 0; pivot < NN_; ++pivot) {
        if (grid_[pivot] != 0 || bits(cand_[pivot]) != 3) continue;
        const unsigned int pm = cand_[pivot];

        std::vector<int> wingCells;
        for (int p : peers_[pivot]) {
            if (grid_[p] != 0 || bits(cand_[p]) != 2) continue;
            if ((cand_[p] & ~pm) == 0U) wingCells.push_back(p);
        }
        if (wingCells.size() < 2U) continue;

        for (std::size_t i = 0; i < wingCells.size(); ++i) {
            const int w1 = wingCells[i];
            const unsigned int m1 = cand_[w1];
            for (std::size_t j = i + 1; j < wingCells.size(); ++j) {
                const int w2 = wingCells[j];
                const unsigned int m2 = cand_[w2];
                if ((m1 | m2) != pm) continue;
                const unsigned int zMask = m1 & m2;
                if (bits(zMask) != 1) continue;
                const int z = firstDigit(zMask);

                int local = 0;
                for (int c = 0; c < NN_; ++c) {
                    if (c == pivot || c == w1 || c == w2) continue;
                    if (grid_[c] != 0 || (cand_[c] & bit(z)) == 0U) continue;
                    if (!isPeerCell(c, pivot) || !isPeerCell(c, w1) || !isPeerCell(c, w2)) continue;
                    bool changed = false;
                    if (!removeCandidate(c, z, changed)) return false;
                    if (changed) ++local;
                }
                if (local > 0) { n = local; return true; }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyWXYZWing(int& n) {
    n = 0;

    for (int pivot = 0; pivot < NN_; ++pivot) {
        if (grid_[pivot] != 0 || bits(cand_[pivot]) != 4) continue;
        const unsigned int pm = cand_[pivot];

        std::vector<int> wings;
        for (int p : peers_[pivot]) {
            if (grid_[p] != 0) continue;
            const int bc = bits(cand_[p]);
            if (bc < 2 || bc > 3) continue;
            if ((cand_[p] & ~pm) == 0U) wings.push_back(p);
        }
        if (wings.size() < 3U) continue;

        for (std::size_t i = 0; i < wings.size(); ++i) {
            const int w1 = wings[i];
            for (std::size_t j = i + 1; j < wings.size(); ++j) {
                const int w2 = wings[j];
                for (std::size_t k = j + 1; k < wings.size(); ++k) {
                    const int w3 = wings[k];
                    const unsigned int um = pm | cand_[w1] | cand_[w2] | cand_[w3];
                    if (um != pm) continue;

                    unsigned int zMask = pm & cand_[w1] & cand_[w2] & cand_[w3];
                    while (zMask) {
                        const unsigned int one = zMask & (~zMask + 1U);
                        const int z = firstDigit(one);
                        zMask &= (zMask - 1U);

                        std::vector<int> zCells;
                        if (pm & bit(z)) zCells.push_back(pivot);
                        if (cand_[w1] & bit(z)) zCells.push_back(w1);
                        if (cand_[w2] & bit(z)) zCells.push_back(w2);
                        if (cand_[w3] & bit(z)) zCells.push_back(w3);
                        if (zCells.size() < 3U) continue;

                        bool covered = true;
                        unsigned int rest = pm & ~bit(z);
                        while (rest) {
                            const unsigned int od = rest & (~rest + 1U);
                            const int d = firstDigit(od);
                            rest &= (rest - 1U);

                            int cnt = 0;
                            if (pm & bit(d)) ++cnt;
                            if (cand_[w1] & bit(d)) ++cnt;
                            if (cand_[w2] & bit(d)) ++cnt;
                            if (cand_[w3] & bit(d)) ++cnt;
                            if (cnt < 2) {
                                covered = false;
                                break;
                            }
                        }
                        if (!covered) continue;

                        int local = 0;
                        for (int c = 0; c < NN_; ++c) {
                            if (c == pivot || c == w1 || c == w2 || c == w3) continue;
                            if (grid_[c] != 0 || (cand_[c] & bit(z)) == 0U) continue;
                            bool seesAll = true;
                            for (int s : zCells) {
                                if (!isPeerCell(c, s)) { seesAll = false; break; }
                            }
                            if (!seesAll) continue;
                            bool changed = false;
                            if (!removeCandidate(c, z, changed)) return false;
                            if (changed) ++local;
                        }
                        if (local > 0) {
                            std::ostringstream ss;
                            ss << "WXYZWing: pivot " << cellName(pivot) << " remove " << z
                               << " from " << local << " cell(s)";
                            pushDebugLog(ss.str());
                            n = local;
                            return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applySkyscraper(int& n) {
    n = 0;

    for (int d = 1; d <= N_; ++d) {
        std::vector<std::pair<int, std::array<int, 2>>> rowPairs;
        for (int r = 0; r < N_; ++r) {
            std::vector<int> cols;
            for (int c = 0; c < N_; ++c) {
                const int idx = r * N_ + c;
                if (grid_[idx] == 0 && (cand_[idx] & bit(d))) cols.push_back(c);
            }
            if (cols.size() == 2U) rowPairs.push_back({r, {cols[0], cols[1]}});
        }

        for (std::size_t i = 0; i < rowPairs.size(); ++i) {
            for (std::size_t j = i + 1; j < rowPairs.size(); ++j) {
                const int r1 = rowPairs[i].first, r2 = rowPairs[j].first;
                const std::array<int, 2>& a = rowPairs[i].second;
                const std::array<int, 2>& b = rowPairs[j].second;
                int shared = -1, oa = -1, ob = -1;
                for (int x : a) {
                    for (int y : b) {
                        if (x == y) shared = x;
                    }
                }
                if (shared < 0) continue;
                oa = (a[0] == shared) ? a[1] : a[0];
                ob = (b[0] == shared) ? b[1] : b[0];
                const int roof1 = r1 * N_ + oa;
                const int roof2 = r2 * N_ + ob;

                int local = 0;
                for (int c = 0; c < NN_; ++c) {
                    if (c == roof1 || c == roof2) continue;
                    if (grid_[c] != 0 || (cand_[c] & bit(d)) == 0U) continue;
                    if (!isPeerCell(c, roof1) || !isPeerCell(c, roof2)) continue;
                    bool changed = false;
                    if (!removeCandidate(c, d, changed)) return false;
                    if (changed) ++local;
                }
                if (local > 0) { n = local; return true; }
            }
        }

        std::vector<std::pair<int, std::array<int, 2>>> colPairs;
        for (int c = 0; c < N_; ++c) {
            std::vector<int> rows;
            for (int r = 0; r < N_; ++r) {
                const int idx = r * N_ + c;
                if (grid_[idx] == 0 && (cand_[idx] & bit(d))) rows.push_back(r);
            }
            if (rows.size() == 2U) colPairs.push_back({c, {rows[0], rows[1]}});
        }
        for (std::size_t i = 0; i < colPairs.size(); ++i) {
            for (std::size_t j = i + 1; j < colPairs.size(); ++j) {
                const int c1 = colPairs[i].first, c2 = colPairs[j].first;
                const std::array<int, 2>& a = colPairs[i].second;
                const std::array<int, 2>& b = colPairs[j].second;
                int shared = -1, oa = -1, ob = -1;
                for (int x : a) {
                    for (int y : b) {
                        if (x == y) shared = x;
                    }
                }
                if (shared < 0) continue;
                oa = (a[0] == shared) ? a[1] : a[0];
                ob = (b[0] == shared) ? b[1] : b[0];
                const int roof1 = oa * N_ + c1;
                const int roof2 = ob * N_ + c2;

                int local = 0;
                for (int c = 0; c < NN_; ++c) {
                    if (c == roof1 || c == roof2) continue;
                    if (grid_[c] != 0 || (cand_[c] & bit(d)) == 0U) continue;
                    if (!isPeerCell(c, roof1) || !isPeerCell(c, roof2)) continue;
                    bool changed = false;
                    if (!removeCandidate(c, d, changed)) return false;
                    if (changed) ++local;
                }
                if (local > 0) { n = local; return true; }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyTwoStringKite(int& n) {
    n = 0;
    for (int d = 1; d <= N_; ++d) {
        std::vector<std::tuple<int, int, int>> rowPairs;
        std::vector<std::tuple<int, int, int>> colPairs;
        for (int r = 0; r < N_; ++r) {
            std::vector<int> cols;
            for (int c = 0; c < N_; ++c) {
                const int idx = r * N_ + c;
                if (grid_[idx] == 0 && (cand_[idx] & bit(d))) cols.push_back(c);
            }
            if (cols.size() == 2U) rowPairs.push_back({r, cols[0], cols[1]});
        }
        for (int c = 0; c < N_; ++c) {
            std::vector<int> rows;
            for (int r = 0; r < N_; ++r) {
                const int idx = r * N_ + c;
                if (grid_[idx] == 0 && (cand_[idx] & bit(d))) rows.push_back(r);
            }
            if (rows.size() == 2U) colPairs.push_back({c, rows[0], rows[1]});
        }

        for (const auto& rp : rowPairs) {
            const int r = std::get<0>(rp);
            const int c1 = std::get<1>(rp), c2 = std::get<2>(rp);
            for (const auto& cp : colPairs) {
                const int c = std::get<0>(cp);
                const int r1 = std::get<1>(cp), r2 = std::get<2>(cp);
                if (c != c1 && c != c2) continue;
                if (r != r1 && r != r2) continue;

                const int otherCol = (c == c1) ? c2 : c1;
                const int otherRow = (r == r1) ? r2 : r1;
                const int a = r * N_ + otherCol;
                const int b = otherRow * N_ + c;

                int local = 0;
                for (int i = 0; i < NN_; ++i) {
                    if (i == a || i == b) continue;
                    if (grid_[i] != 0 || (cand_[i] & bit(d)) == 0U) continue;
                    if (!isPeerCell(i, a) || !isPeerCell(i, b)) continue;
                    bool changed = false;
                    if (!removeCandidate(i, d, changed)) return false;
                    if (changed) ++local;
                }
                if (local > 0) { n = local; return true; }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applySimpleColoring(int& n) {
    n = 0;
    for (int d = 1; d <= N_; ++d) {
        std::vector<std::vector<int>> g(NN_);
        std::vector<char> hasNode(NN_, 0);
        for (const auto& h : houses_) {
            std::vector<int> w;
            for (int i : h) if (grid_[i] == 0 && (cand_[i] & bit(d))) w.push_back(i);
            if (w.size() == 2U) {
                const int a = w[0], b = w[1];
                g[a].push_back(b);
                g[b].push_back(a);
                hasNode[a] = hasNode[b] = 1;
            }
        }

        std::vector<int> color(NN_, -1);
        std::vector<int> comp(NN_, -1);
        std::vector<std::array<std::vector<int>, 2>> compNodes;

        int compId = 0;
        for (int s = 0; s < NN_; ++s) {
            if (!hasNode[s] || color[s] != -1) continue;
            compNodes.push_back({});
            std::vector<int> q = {s};
            color[s] = 0;
            comp[s] = compId;
            for (std::size_t qi = 0; qi < q.size(); ++qi) {
                const int u = q[qi];
                compNodes[compId][color[u]].push_back(u);
                for (int v : g[u]) {
                    if (color[v] == -1) {
                        color[v] = 1 - color[u];
                        comp[v] = compId;
                        q.push_back(v);
                    }
                }
            }
            ++compId;
        }

        int local = 0;
        for (int i = 0; i < NN_; ++i) {
            if (grid_[i] != 0 || (cand_[i] & bit(d)) == 0U || color[i] != -1) continue;
            std::map<int, int> seen;
            for (int p : peers_[i]) {
                if (color[p] == -1 || comp[p] < 0) continue;
                seen[comp[p]] |= (1 << color[p]);
            }
            bool removed = false;
            for (const auto& it : seen) {
                if (it.second == 3) {
                    bool changed = false;
                    if (!removeCandidate(i, d, changed)) return false;
                    if (changed) { ++local; removed = true; }
                    break;
                }
            }
            if (removed) continue;
        }
        if (local > 0) { n = local; return true; }

        for (int cid = 0; cid < compId; ++cid) {
            for (int clr = 0; clr < 2; ++clr) {
                bool badColor = false;
                for (const auto& h : houses_) {
                    int cnt = 0;
                    for (int i : h) {
                        if (comp[i] == cid && color[i] == clr && grid_[i] == 0 && (cand_[i] & bit(d))) {
                            ++cnt;
                            if (cnt >= 2) { badColor = true; break; }
                        }
                    }
                    if (badColor) break;
                }
                if (!badColor) continue;

                int removed = 0;
                for (int i : compNodes[cid][clr]) {
                    bool changed = false;
                    if (!removeCandidate(i, d, changed)) return false;
                    if (changed) ++removed;
                }
                if (removed > 0) { n = removed; return true; }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyUniqueRectangleType1(int& n) {
    n = 0;
    for (int r1 = 0; r1 < N_; ++r1) {
        for (int r2 = r1 + 1; r2 < N_; ++r2) {
            for (int c1 = 0; c1 < N_; ++c1) {
                for (int c2 = c1 + 1; c2 < N_; ++c2) {
                    const std::array<int, 4> cells = {
                        r1 * N_ + c1, r1 * N_ + c2, r2 * N_ + c1, r2 * N_ + c2
                    };
                    bool ok = true;
                    for (int idx : cells) {
                        if (grid_[idx] != 0 || bits(cand_[idx]) < 2) { ok = false; break; }
                    }
                    if (!ok) continue;

                    std::set<unsigned int> pairMasks;
                    for (int idx : cells) if (bits(cand_[idx]) == 2) pairMasks.insert(cand_[idx]);
                    for (unsigned int pairMask : pairMasks) {
                        int exact = 0;
                        int extraCell = -1;
                        bool valid = true;
                        for (int idx : cells) {
                            const unsigned int m = cand_[idx];
                            if ((m & pairMask) != pairMask) { valid = false; break; }
                            if (m == pairMask) {
                                ++exact;
                            } else if ((m & ~pairMask) != 0U && extraCell == -1) {
                                extraCell = idx;
                            } else {
                                valid = false;
                                break;
                            }
                        }
                        if (!valid || exact != 3 || extraCell < 0) continue;

                        int local = 0;
                        unsigned int rm = pairMask;
                        while (rm) {
                            const unsigned int one = rm & (~rm + 1U);
                            bool changed = false;
                            if (!removeCandidate(extraCell, firstDigit(one), changed)) return false;
                            if (changed) ++local;
                            rm &= (rm - 1U);
                        }
                        if (local > 0) { n = local; return true; }
                    }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyUniqueRectangleType2to6(int& n) {
    n = 0;

    auto digitsInRow = [&](int r, int d) {
        std::vector<int> out;
        for (int c = 0; c < N_; ++c) {
            const int idx = r * N_ + c;
            if (grid_[idx] == 0 && (cand_[idx] & bit(d))) out.push_back(idx);
        }
        return out;
    };
    auto digitsInCol = [&](int c, int d) {
        std::vector<int> out;
        for (int r = 0; r < N_; ++r) {
            const int idx = r * N_ + c;
            if (grid_[idx] == 0 && (cand_[idx] & bit(d))) out.push_back(idx);
        }
        return out;
    };
    auto removeFromSharedLine = [&](int roofA, int roofB, int z, int& local) -> bool {
        if (row(roofA) == row(roofB)) {
            const int r = row(roofA);
            for (int idx : houses_[r]) {
                if (idx == roofA || idx == roofB || grid_[idx] != 0) continue;
                bool changed = false;
                if (!removeCandidate(idx, z, changed)) return false;
                if (changed) ++local;
            }
            return true;
        }
        if (col(roofA) == col(roofB)) {
            const int c = col(roofA);
            for (int idx : houses_[N_ + c]) {
                if (idx == roofA || idx == roofB || grid_[idx] != 0) continue;
                bool changed = false;
                if (!removeCandidate(idx, z, changed)) return false;
                if (changed) ++local;
            }
            return true;
        }
        return true;
    };

    for (int r1 = 0; r1 < N_; ++r1) {
        for (int r2 = r1 + 1; r2 < N_; ++r2) {
            for (int c1 = 0; c1 < N_; ++c1) {
                for (int c2 = c1 + 1; c2 < N_; ++c2) {
                    const int a = r1 * N_ + c1;
                    const int b = r1 * N_ + c2;
                    const int c = r2 * N_ + c1;
                    const int d = r2 * N_ + c2;
                    const std::array<int, 4> cells = {a, b, c, d};

                    bool allUnsolved = true;
                    for (int idx : cells) {
                        if (grid_[idx] != 0 || bits(cand_[idx]) < 2) {
                            allUnsolved = false;
                            break;
                        }
                    }
                    if (!allUnsolved) continue;

                    const unsigned int common = cand_[a] & cand_[b] & cand_[c] & cand_[d];
                    const std::vector<int> commonDigits = digitsFromMask(common);
                    if (commonDigits.size() < 2U) continue;

                    for (std::size_t i = 0; i < commonDigits.size(); ++i) {
                        for (std::size_t j = i + 1; j < commonDigits.size(); ++j) {
                            const int p1 = commonDigits[i];
                            const int p2 = commonDigits[j];
                            const unsigned int pairMask = bit(p1) | bit(p2);

                            bool pairInAll = true;
                            for (int idx : cells) {
                                if ((cand_[idx] & pairMask) != pairMask) {
                                    pairInAll = false;
                                    break;
                                }
                            }
                            if (!pairInAll) continue;

                            std::vector<int> floors;
                            std::vector<int> roofs;
                            for (int idx : cells) {
                                if (cand_[idx] == pairMask) floors.push_back(idx);
                                else roofs.push_back(idx);
                            }
                            if (floors.size() != 2U || roofs.size() != 2U) continue;

                            const int roofA = roofs[0], roofB = roofs[1];
                            const unsigned int extraA = cand_[roofA] & ~pairMask;
                            const unsigned int extraB = cand_[roofB] & ~pairMask;
                            if (extraA == 0U || extraB == 0U) continue;

                            // Hidden UR (konserwatywnie): para (p1,p2) jest ukryta w obu wierszach i kolumnach
                            // prostokata, wiec dodatki poza para mozna usunac z tych 4 komorek.
                            auto inSet2 = [&](const std::vector<int>& where, int u, int v) {
                                return where.size() == 2U &&
                                       ((where[0] == u && where[1] == v) || (where[0] == v && where[1] == u));
                            };
                            auto whereInRow = [&](int rr, int dig) {
                                std::vector<int> where;
                                for (int idx : houses_[rr]) {
                                    if (grid_[idx] == 0 && (cand_[idx] & bit(dig))) where.push_back(idx);
                                }
                                return where;
                            };
                            auto whereInCol = [&](int cc, int dig) {
                                std::vector<int> where;
                                for (int idx : houses_[N_ + cc]) {
                                    if (grid_[idx] == 0 && (cand_[idx] & bit(dig))) where.push_back(idx);
                                }
                                return where;
                            };
                            bool hiddenUr = true;
                            for (int dig : {p1, p2}) {
                                if (!inSet2(whereInRow(r1, dig), a, b)) hiddenUr = false;
                                if (!inSet2(whereInRow(r2, dig), c, d)) hiddenUr = false;
                                if (!inSet2(whereInCol(c1, dig), a, c)) hiddenUr = false;
                                if (!inSet2(whereInCol(c2, dig), b, d)) hiddenUr = false;
                            }
                            if (hiddenUr) {
                                int local = 0;
                                for (int idx : cells) {
                                    unsigned int extras = cand_[idx] & ~pairMask;
                                    while (extras) {
                                        const unsigned int one = extras & (~extras + 1U);
                                        bool changed = false;
                                        if (!removeCandidate(idx, firstDigit(one), changed)) return false;
                                        if (changed) ++local;
                                        extras &= (extras - 1U);
                                    }
                                }
                                if (local > 0) {
                                    pushDebugLog("HiddenUR: remove extras outside pair in rectangle");
                                    n = local;
                                    return true;
                                }
                            }

                            // Type 2: dwa roofy w jednej linii, ten sam dodatkowy kandydat.
                            if ((row(roofA) == row(roofB) || col(roofA) == col(roofB)) &&
                                bits(extraA) == 1 && extraA == extraB) {
                                const int z = firstDigit(extraA);
                                int local = 0;
                                if (!removeFromSharedLine(roofA, roofB, z, local)) return false;
                                if (local > 0) { n = local; return true; }
                            }

                            // Type 4: silne lacze na jednej cyfrze pary w domu wspolnym dla roofow.
                            const std::vector<int> pairDigits = {p1, p2};
                            for (int p : pairDigits) {
                                const int other = (p == p1) ? p2 : p1;
                                std::vector<int> housesToCheck;
                                if (row(roofA) == row(roofB)) housesToCheck.push_back(row(roofA));
                                if (col(roofA) == col(roofB)) housesToCheck.push_back(N_ + col(roofA));
                                if (box(row(roofA), col(roofA)) == box(row(roofB), col(roofB))) {
                                    housesToCheck.push_back(2 * N_ + box(row(roofA), col(roofA)));
                                }
                                for (int h : housesToCheck) {
                                    std::vector<int> where;
                                    for (int idx : houses_[h]) {
                                        if (grid_[idx] == 0 && (cand_[idx] & bit(p))) where.push_back(idx);
                                    }
                                    if (where.size() != 2U) continue;
                                    const bool sameTwo =
                                        ((where[0] == roofA && where[1] == roofB) ||
                                         (where[0] == roofB && where[1] == roofA));
                                    if (!sameTwo) continue;

                                    int local = 0;
                                    for (int r : roofs) {
                                        bool changed = false;
                                        if (!removeCandidate(r, other, changed)) return false;
                                        if (changed) ++local;
                                    }
                                    if (local > 0) { n = local; return true; }
                                }
                            }

                            // Type 3 (konserwatywny): eliminacje dodatkow z roofow tylko jesli zalozenie
                            // roof=dodatkowa_cyfra prowadzi do sprzecznosci na poziomie propagacji logicznej.
                            if (row(roofA) == row(roofB) || col(roofA) == col(roofB)) {
                                const unsigned int extrasUnion = extraA | extraB;
                                if (bits(extrasUnion) >= 2) {
                                    int local = 0;
                                    for (int r : roofs) {
                                        unsigned int em = cand_[r] & ~pairMask;
                                        while (em) {
                                            const unsigned int one = em & (~em + 1U);
                                            const int z = firstDigit(one);
                                            if (!hasLogicalSupportWithAssignments({{r, z}})) {
                                                bool changed = false;
                                                if (!removeCandidate(r, z, changed)) return false;
                                                if (changed) ++local;
                                            }
                                            em &= (em - 1U);
                                        }
                                    }
                                    if (local > 0) {
                                        std::ostringstream ss;
                                        ss << "UR Type3: remove roof extras in "
                                           << cellName(roofA) << "/" << cellName(roofB);
                                        pushDebugLog(ss.str());
                                        n = local;
                                        return true;
                                    }
                                }
                            }

                            // Type 5 (konserwatywny): jesli cyfra pary jest silnie zwiazana z oboma floorami
                            // (koniugacje w obu osiach dla kazdego floora), usun z peerow obu floorow.
                            auto strongAroundFloor = [&](int floorCell, int digit) -> bool {
                                int rowCnt = 0, colCnt = 0;
                                const int rr = row(floorCell), cc = col(floorCell);
                                for (int idx : houses_[rr]) {
                                    if (grid_[idx] == 0 && (cand_[idx] & bit(digit))) ++rowCnt;
                                }
                                for (int idx : houses_[N_ + cc]) {
                                    if (grid_[idx] == 0 && (cand_[idx] & bit(digit))) ++colCnt;
                                }
                                return (rowCnt == 2 && colCnt == 2);
                            };
                            for (int p : pairDigits) {
                                if (!strongAroundFloor(floors[0], p) || !strongAroundFloor(floors[1], p)) continue;
                                int local = 0;
                                for (int idx = 0; idx < NN_; ++idx) {
                                    if (idx == floors[0] || idx == floors[1]) continue;
                                    if (grid_[idx] != 0 || (cand_[idx] & bit(p)) == 0U) continue;
                                    if (!isPeerCell(idx, floors[0]) || !isPeerCell(idx, floors[1])) continue;
                                    bool changed = false;
                                    if (!removeCandidate(idx, p, changed)) return false;
                                    if (changed) ++local;
                                }
                                if (local > 0) {
                                    std::ostringstream ss;
                                    ss << "UR Type5: remove " << p << " from peers of floors "
                                       << cellName(floors[0]) << "/" << cellName(floors[1]);
                                    pushDebugLog(ss.str());
                                    n = local;
                                    return true;
                                }
                            }

                            // Type 6 (wariant X-link): diagonalne roofy + silne lacza tej samej cyfry pary na obu wierszach i kolumnach.
                            const bool diagonalRoofs =
                                ((roofA == a && roofB == d) || (roofA == d && roofB == a) ||
                                 (roofA == b && roofB == c) || (roofA == c && roofB == b));
                            if (diagonalRoofs) {
                                for (int p : pairDigits) {
                                    const std::vector<int> r1pos = digitsInRow(r1, p);
                                    const std::vector<int> r2pos = digitsInRow(r2, p);
                                    const std::vector<int> c1pos = digitsInCol(c1, p);
                                    const std::vector<int> c2pos = digitsInCol(c2, p);
                                    if (r1pos.size() != 2U || r2pos.size() != 2U || c1pos.size() != 2U || c2pos.size() != 2U) continue;
                                    bool rowOk = true, colOk = true;
                                    for (int idx : r1pos) if (!(idx == a || idx == b)) rowOk = false;
                                    for (int idx : r2pos) if (!(idx == c || idx == d)) rowOk = false;
                                    for (int idx : c1pos) if (!(idx == a || idx == c)) colOk = false;
                                    for (int idx : c2pos) if (!(idx == b || idx == d)) colOk = false;
                                    if (!rowOk || !colOk) continue;

                                    const int other = (p == p1) ? p2 : p1;
                                    int local = 0;
                                    for (int r : roofs) {
                                        bool changed = false;
                                        if (!removeCandidate(r, other, changed)) return false;
                                        if (changed) ++local;
                                    }
                                    if (local > 0) { n = local; return true; }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyUniqueLoop(int& n) {
    n = 0;

    for (int p1 = 1; p1 <= N_; ++p1) {
        for (int p2 = p1 + 1; p2 <= N_; ++p2) {
            const unsigned int pairMask = bit(p1) | bit(p2);
            std::vector<int> nodes;
            std::vector<int> idByCell(NN_, -1);
            for (int idx = 0; idx < NN_; ++idx) {
                if (grid_[idx] != 0) continue;
                if ((cand_[idx] & pairMask) != pairMask) continue;
                const int bc = bits(cand_[idx]);
                if (bc < 2 || bc > 3) continue;
                idByCell[idx] = static_cast<int>(nodes.size());
                nodes.push_back(idx);
            }
            if (nodes.size() < 6U) continue;

            const int M = static_cast<int>(nodes.size());
            std::vector<std::vector<int>> adj(M);
            auto addEdge = [&](int a, int b) {
                if (a < 0 || b < 0 || a == b) return;
                if (std::find(adj[a].begin(), adj[a].end(), b) == adj[a].end()) adj[a].push_back(b);
                if (std::find(adj[b].begin(), adj[b].end(), a) == adj[b].end()) adj[b].push_back(a);
            };
            for (int h = 0; h < 3 * N_; ++h) {
                std::vector<int> where;
                for (int idx : houses_[h]) {
                    if (idByCell[idx] >= 0 && (cand_[idx] & pairMask) == pairMask) where.push_back(idByCell[idx]);
                }
                if (where.size() == 2U) addEdge(where[0], where[1]);
            }

            std::vector<char> used(M, 0);
            std::vector<int> path;
            std::function<bool(int, int)> dfs = [&](int start, int u) -> bool {
                if (path.size() > 14U) return false;
                for (int v : adj[u]) {
                    if (v == start) {
                        const std::size_t L = path.size();
                        if (L < 6U || (L % 2U) != 0U) continue;

                        int extraNode = -1;
                        bool ok = true;
                        for (int nid : path) {
                            const unsigned int extras = cand_[nodes[nid]] & ~pairMask;
                            if (extras == 0U) continue;
                            if (bits(extras) > 1 || extraNode != -1) { ok = false; break; }
                            extraNode = nid;
                        }
                        if (!ok || extraNode < 0) continue;

                        int local = 0;
                        bool ch1 = false, ch2 = false;
                        if (!removeCandidate(nodes[extraNode], p1, ch1)) return false;
                        if (ch1) ++local;
                        if (!removeCandidate(nodes[extraNode], p2, ch2)) return false;
                        if (ch2) ++local;
                        if (local > 0) {
                            std::ostringstream ss;
                            ss << "UniqueLoop: cycle length " << L
                               << " remove {" << p1 << "," << p2 << "} from " << cellName(nodes[extraNode]);
                            pushDebugLog(ss.str());
                            n = local;
                            return true;
                        }
                        continue;
                    }
                    if (used[v]) continue;
                    used[v] = 1;
                    path.push_back(v);
                    if (dfs(start, v)) return true;
                    path.pop_back();
                    used[v] = 0;
                }
                return false;
            };

            for (int s = 0; s < M; ++s) {
                std::fill(used.begin(), used.end(), 0);
                path.clear();
                used[s] = 1;
                path.push_back(s);
                if (dfs(s, s)) return true;
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyBivalueOddagon(int& n) {
    n = 0;

    for (int p1 = 1; p1 <= N_; ++p1) {
        for (int p2 = p1 + 1; p2 <= N_; ++p2) {
            const unsigned int pairMask = bit(p1) | bit(p2);
            std::vector<int> nodes;
            std::vector<int> idByCell(NN_, -1);
            for (int idx = 0; idx < NN_; ++idx) {
                if (grid_[idx] != 0) continue;
                if ((cand_[idx] & pairMask) != pairMask) continue;
                const int bc = bits(cand_[idx]);
                if (bc < 2 || bc > 3) continue;
                idByCell[idx] = static_cast<int>(nodes.size());
                nodes.push_back(idx);
            }
            if (nodes.size() < 5U) continue;

            const int M = static_cast<int>(nodes.size());
            std::vector<std::vector<int>> adj(M);
            auto addEdge = [&](int a, int b) {
                if (a < 0 || b < 0 || a == b) return;
                if (std::find(adj[a].begin(), adj[a].end(), b) == adj[a].end()) adj[a].push_back(b);
                if (std::find(adj[b].begin(), adj[b].end(), a) == adj[b].end()) adj[b].push_back(a);
            };
            for (int h = 0; h < 3 * N_; ++h) {
                std::vector<int> where;
                for (int idx : houses_[h]) {
                    if (idByCell[idx] >= 0 && (cand_[idx] & pairMask) == pairMask) where.push_back(idByCell[idx]);
                }
                if (where.size() == 2U) addEdge(where[0], where[1]);
            }

            std::vector<char> used(M, 0);
            std::vector<int> path;
            std::function<bool(int, int)> dfs = [&](int start, int u) -> bool {
                if (path.size() > 13U) return false;
                for (int v : adj[u]) {
                    if (v == start) {
                        const std::size_t L = path.size();
                        if (L < 5U || (L % 2U) == 0U) continue;

                        int extraNode = -1;
                        bool ok = true;
                        for (int nid : path) {
                            const unsigned int extras = cand_[nodes[nid]] & ~pairMask;
                            if (extras == 0U) continue;
                            if (bits(extras) > 1 || extraNode != -1) { ok = false; break; }
                            extraNode = nid;
                        }
                        if (!ok || extraNode < 0) continue;

                        int local = 0;
                        bool ch1 = false, ch2 = false;
                        if (!removeCandidate(nodes[extraNode], p1, ch1)) return false;
                        if (ch1) ++local;
                        if (!removeCandidate(nodes[extraNode], p2, ch2)) return false;
                        if (ch2) ++local;
                        if (local > 0) {
                            std::ostringstream ss;
                            ss << "BivalueOddagon: cycle length " << L
                               << " remove {" << p1 << "," << p2 << "} from " << cellName(nodes[extraNode]);
                            pushDebugLog(ss.str());
                            n = local;
                            return true;
                        }
                        continue;
                    }
                    if (used[v]) continue;
                    used[v] = 1;
                    path.push_back(v);
                    if (dfs(start, v)) return true;
                    path.pop_back();
                    used[v] = 0;
                }
                return false;
            };

            for (int s = 0; s < M; ++s) {
                std::fill(used.begin(), used.end(), 0);
                path.clear();
                used[s] = 1;
                path.push_back(s);
                if (dfs(s, s)) return true;
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyBUGPlus1(int& n) {
    n = 0;
    int bugCell = -1;
    for (int i = 0; i < NN_; ++i) {
        if (grid_[i] != 0) continue;
        const int bc = bits(cand_[i]);
        if (bc == 2) continue;
        if (bc == 3 && bugCell == -1) {
            bugCell = i;
            continue;
        }
        return false;
    }
    if (bugCell < 0) return false;

    const int r = row(bugCell), c = col(bugCell), b = box(r, c);
    for (int d : digitsFromMask(cand_[bugCell])) {
        int rowCnt = 0, colCnt = 0, boxCnt = 0;
        for (int i : houses_[r]) if (grid_[i] == 0 && (cand_[i] & bit(d))) ++rowCnt;
        for (int i : houses_[N_ + c]) if (grid_[i] == 0 && (cand_[i] & bit(d))) ++colCnt;
        for (int i : houses_[2 * N_ + b]) if (grid_[i] == 0 && (cand_[i] & bit(d))) ++boxCnt;
        if ((rowCnt % 2) == 1 && (colCnt % 2) == 1 && (boxCnt % 2) == 1) {
            if (!assignValue(bugCell, d)) return false;
            n = 1;
            return true;
        }
    }
    return false;
}

bool SudokuAnalyzer::applyRemotePairs(int& n) {
    n = 0;
    std::map<unsigned int, std::vector<int>> byPair;
    for (int i = 0; i < NN_; ++i) {
        if (grid_[i] == 0 && bits(cand_[i]) == 2) {
            byPair[cand_[i]].push_back(i);
        }
    }

    for (const auto& entry : byPair) {
        const unsigned int pairMask = entry.first;
        const std::vector<int>& nodes = entry.second;
        if (nodes.size() < 4U) continue;
        std::map<int, int> nodePos;
        for (std::size_t i = 0; i < nodes.size(); ++i) nodePos[nodes[i]] = static_cast<int>(i);

        std::vector<std::vector<int>> g(nodes.size());
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            for (std::size_t j = i + 1; j < nodes.size(); ++j) {
                if (!isPeerCell(nodes[i], nodes[j])) continue;
                g[i].push_back(static_cast<int>(j));
                g[j].push_back(static_cast<int>(i));
            }
        }

        std::vector<int> color(nodes.size(), -1), comp(nodes.size(), -1);
        int compId = 0;
        std::vector<std::array<std::vector<int>, 2>> compNodes;
        for (std::size_t s = 0; s < nodes.size(); ++s) {
            if (color[s] != -1) continue;
            compNodes.push_back({});
            std::vector<int> q = {static_cast<int>(s)};
            color[s] = 0;
            comp[s] = compId;
            for (std::size_t qi = 0; qi < q.size(); ++qi) {
                const int u = q[qi];
                compNodes[compId][color[u]].push_back(nodes[u]);
                for (int v : g[u]) {
                    if (color[v] == -1) {
                        color[v] = 1 - color[u];
                        comp[v] = compId;
                        q.push_back(v);
                    }
                }
            }
            ++compId;
        }

        int local = 0;
        for (int i = 0; i < NN_; ++i) {
            if (grid_[i] != 0) continue;
            if ((cand_[i] & pairMask) == 0U) continue;
            std::map<int, int> seen;
            for (std::size_t p = 0; p < nodes.size(); ++p) {
                if (!isPeerCell(i, nodes[p])) continue;
                seen[comp[p]] |= (1 << color[p]);
            }
            bool validComp = false;
            for (const auto& it : seen) {
                if (it.second == 3) { validComp = true; break; }
            }
            if (!validComp) continue;

            unsigned int rm = pairMask & cand_[i];
            while (rm) {
                const unsigned int one = rm & (~rm + 1U);
                bool changed = false;
                if (!removeCandidate(i, firstDigit(one), changed)) return false;
                if (changed) ++local;
                rm &= (rm - 1U);
            }
        }
        if (local > 0) { n = local; return true; }
    }
    return false;
}

bool SudokuAnalyzer::applyWWing(int& n) {
    n = 0;
    std::vector<std::vector<std::pair<int, int>>> conjugates(N_ + 1);
    for (int d = 1; d <= N_; ++d) {
        for (const auto& h : houses_) {
            std::vector<int> w;
            for (int i : h) if (grid_[i] == 0 && (cand_[i] & bit(d))) w.push_back(i);
            if (w.size() == 2U) conjugates[d].push_back({w[0], w[1]});
        }
    }

    std::vector<int> bivalueCells;
    for (int i = 0; i < NN_; ++i) if (grid_[i] == 0 && bits(cand_[i]) == 2) bivalueCells.push_back(i);

    for (std::size_t i = 0; i < bivalueCells.size(); ++i) {
        const int p = bivalueCells[i];
        for (std::size_t j = i + 1; j < bivalueCells.size(); ++j) {
            const int q = bivalueCells[j];
            if (cand_[p] != cand_[q] || isPeerCell(p, q)) continue;
            const std::vector<int> ds = digitsFromMask(cand_[p]);
            if (ds.size() != 2U) continue;

            for (int linkDigit : ds) {
                const int elimDigit = (ds[0] == linkDigit) ? ds[1] : ds[0];
                for (const auto& lk : conjugates[linkDigit]) {
                    const int a = lk.first, b = lk.second;
                    const bool ok1 = isPeerCell(p, a) && isPeerCell(q, b);
                    const bool ok2 = isPeerCell(p, b) && isPeerCell(q, a);
                    if (!ok1 && !ok2) continue;

                    int local = 0;
                    for (int c = 0; c < NN_; ++c) {
                        if (c == p || c == q) continue;
                        if (grid_[c] != 0 || (cand_[c] & bit(elimDigit)) == 0U) continue;
                        if (!isPeerCell(c, p) || !isPeerCell(c, q)) continue;
                        bool changed = false;
                        if (!removeCandidate(c, elimDigit, changed)) return false;
                        if (changed) ++local;
                    }
                    if (local > 0) { n = local; return true; }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyFinnedXWingSashimi(int& n) {
    n = 0;

    auto process = [&](bool rowBased) -> bool {
        for (int d = 1; d <= N_; ++d) {
            for (int l1 = 0; l1 < N_; ++l1) {
                std::vector<int> p1;
                for (int p = 0; p < N_; ++p) {
                    const int idx = rowBased ? (l1 * N_ + p) : (p * N_ + l1);
                    if (grid_[idx] == 0 && (cand_[idx] & bit(d))) p1.push_back(p);
                }
                if (p1.size() != 2U) continue;

                for (int l2 = 0; l2 < N_; ++l2) {
                    if (l2 == l1) continue;
                    std::vector<int> p2;
                    for (int p = 0; p < N_; ++p) {
                        const int idx = rowBased ? (l2 * N_ + p) : (p * N_ + l2);
                        if (grid_[idx] == 0 && (cand_[idx] & bit(d))) p2.push_back(p);
                    }
                    if (p2.size() < 3U || p2.size() > 4U) continue;

                    std::set<int> baseSet(p1.begin(), p1.end());
                    std::vector<int> extras;
                    for (int p : p2) {
                        if (!baseSet.count(p)) extras.push_back(p);
                    }
                    if (extras.empty()) continue;

                    int finBox = -1;
                    bool sameBox = true;
                    for (int ep : extras) {
                        const int b = rowBased ? box(l2, ep) : box(ep, l2);
                        if (finBox < 0) finBox = b;
                        else if (finBox != b) { sameBox = false; break; }
                    }
                    if (!sameBox || finBox < 0) continue;

                    for (int basePos : p1) {
                        const int baseInL2Box = rowBased ? box(l2, basePos) : box(basePos, l2);
                        if (baseInL2Box != finBox) continue;

                        int local = 0;
                        const std::vector<int>& targetHouse = rowBased ? houses_[N_ + basePos] : houses_[basePos];
                        for (int idx : targetHouse) {
                            const int line = rowBased ? row(idx) : col(idx);
                            if (line == l1 || line == l2) continue;
                            if ((cand_[idx] & bit(d)) == 0U || grid_[idx] != 0) continue;
                            if (box(row(idx), col(idx)) != finBox) continue;
                            bool changed = false;
                            if (!removeCandidate(idx, d, changed)) return false;
                            if (changed) ++local;
                        }
                        if (local > 0) { n = local; return true; }
                    }
                }
            }
        }
        return false;
    };

    if (process(true)) return true;
    return process(false);
}

bool SudokuAnalyzer::applyFinnedFish(int size, int& n) {
    n = 0;
    if (size < 2 || size > N_) return false;

    auto process = [&](bool rowBased) -> bool {
        for (int d = 1; d <= N_; ++d) {
            std::vector<std::vector<int>> positions(N_);
            std::vector<int> eligible;
            for (int line = 0; line < N_; ++line) {
                std::vector<int> pos;
                for (int p = 0; p < N_; ++p) {
                    const int idx = rowBased ? (line * N_ + p) : (p * N_ + line);
                    if (grid_[idx] == 0 && (cand_[idx] & bit(d))) pos.push_back(p);
                }
                positions[line] = std::move(pos);
                const int cnt = static_cast<int>(positions[line].size());
                if (cnt >= 2 && cnt <= size + 2) eligible.push_back(line);
            }
            if (static_cast<int>(eligible.size()) < size) continue;

            bool found = false;
            forEachCombo(eligible, size, [&](const std::vector<int>& lines) {
                if (found || contradiction_) return;

                std::vector<char> inLine(N_, 0), inUnionPos(N_, 0);
                std::vector<int> unionPos;
                for (int line : lines) {
                    inLine[line] = 1;
                    for (int p : positions[line]) {
                        if (!inUnionPos[p]) {
                            inUnionPos[p] = 1;
                            unionPos.push_back(p);
                        }
                    }
                }
                const int uc = static_cast<int>(unionPos.size());
                if (uc <= size || uc > size + 2) return;
                if (uc < size) return;

                forEachCombo(unionPos, size, [&](const std::vector<int>& coverPos) {
                    if (found || contradiction_) return;
                    std::vector<char> inCover(N_, 0);
                    for (int p : coverPos) inCover[p] = 1;

                    std::vector<int> finIndices;
                    bool valid = true;
                    for (int line : lines) {
                        int baseCnt = 0;
                        int finCnt = 0;
                        for (int p : positions[line]) {
                            if (inCover[p]) {
                                ++baseCnt;
                            } else {
                                ++finCnt;
                                const int finIdx = rowBased ? (line * N_ + p) : (p * N_ + line);
                                finIndices.push_back(finIdx);
                            }
                        }
                        if (baseCnt < 2 || finCnt > 2) {
                            valid = false;
                            break;
                        }
                    }
                    if (!valid || finIndices.empty()) return;

                    int finBox = -1;
                    for (int idx : finIndices) {
                        const int b = box(row(idx), col(idx));
                        if (finBox < 0) finBox = b;
                        else if (finBox != b) { valid = false; break; }
                    }
                    if (!valid || finBox < 0) return;

                    int local = 0;
                    for (int p : coverPos) {
                        const std::vector<int>& coverHouse = rowBased ? houses_[N_ + p] : houses_[p];
                        for (int idx : coverHouse) {
                            const int line = rowBased ? row(idx) : col(idx);
                            if (inLine[line]) continue;
                            if (grid_[idx] != 0 || (cand_[idx] & bit(d)) == 0U) continue;
                            if (box(row(idx), col(idx)) != finBox) continue;
                            bool seesAllFins = true;
                            for (int f : finIndices) {
                                if (!isPeerCell(idx, f)) { seesAllFins = false; break; }
                            }
                            if (!seesAllFins) continue;
                            bool changed = false;
                            if (!removeCandidate(idx, d, changed)) return;
                            if (changed) ++local;
                        }
                    }
                    if (local > 0) {
                        std::ostringstream ss;
                        ss << "FinnedFish(" << size << "): remove " << d
                           << " in fin box " << (finBox + 1)
                           << " from " << local << " cell(s)";
                        pushDebugLog(ss.str());
                        n = local;
                        found = true;
                    }
                });
            });
            if (found) return true;
        }
        return false;
    };

    if (process(true)) return true;
    return process(false);
}

bool SudokuAnalyzer::applyEmptyRectangle(int& n) {
    n = 0;
    for (int d = 1; d <= N_; ++d) {
        for (int b = 0; b < N_; ++b) {
            const auto& boxCells = houses_[2 * N_ + b];
            std::vector<int> candInBox;
            std::set<int> rowsInBox, colsInBox;
            for (int idx : boxCells) {
                if (grid_[idx] == 0 && (cand_[idx] & bit(d))) {
                    candInBox.push_back(idx);
                    rowsInBox.insert(row(idx));
                    colsInBox.insert(col(idx));
                }
            }
            if (candInBox.size() < 3U) continue;

            for (int r : rowsInBox) {
                for (int c : colsInBox) {
                    const int crossIdx = r * N_ + c;
                    if (box(r, c) != b) continue;
                    if (grid_[crossIdx] == 0 && (cand_[crossIdx] & bit(d))) continue;

                    bool crossShape = true;
                    std::vector<int> rowArm, colArm;
                    for (int idx : candInBox) {
                        if (row(idx) != r && col(idx) != c) { crossShape = false; break; }
                        if (row(idx) == r && col(idx) != c) rowArm.push_back(idx);
                        if (col(idx) == c && row(idx) != r) colArm.push_back(idx);
                    }
                    if (!crossShape || rowArm.empty() || colArm.empty()) continue;

                    std::vector<int> rowCandidates;
                    for (int idx : houses_[r]) {
                        if (grid_[idx] == 0 && (cand_[idx] & bit(d))) rowCandidates.push_back(idx);
                    }
                    std::vector<int> colCandidates;
                    for (int idx : houses_[N_ + c]) {
                        if (grid_[idx] == 0 && (cand_[idx] & bit(d))) colCandidates.push_back(idx);
                    }
                    if (rowCandidates.size() != 2U || colCandidates.size() != 2U) continue;

                    int rowInside = -1, rowOutside = -1;
                    for (int idx : rowCandidates) {
                        if (box(row(idx), col(idx)) == b) rowInside = idx;
                        else rowOutside = idx;
                    }
                    int colInside = -1, colOutside = -1;
                    for (int idx : colCandidates) {
                        if (box(row(idx), col(idx)) == b) colInside = idx;
                        else colOutside = idx;
                    }
                    if (rowInside < 0 || rowOutside < 0 || colInside < 0 || colOutside < 0) continue;

                    const int elimIdx = row(colOutside) * N_ + col(rowOutside);
                    if (row(elimIdx) != row(colOutside) || col(elimIdx) != col(rowOutside)) continue;
                    if (elimIdx == rowOutside || elimIdx == colOutside) continue;
                    if (grid_[elimIdx] != 0 || (cand_[elimIdx] & bit(d)) == 0U) continue;

                    bool changed = false;
                    if (!removeCandidate(elimIdx, d, changed)) return false;
                    if (changed) { n = 1; return true; }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyAvoidableRectangle(int& n) {
    n = 0;
    for (int r1 = 0; r1 < N_; ++r1) {
        for (int r2 = r1 + 1; r2 < N_; ++r2) {
            for (int c1 = 0; c1 < N_; ++c1) {
                for (int c2 = c1 + 1; c2 < N_; ++c2) {
                    const std::array<int, 4> cells = {
                        r1 * N_ + c1, r1 * N_ + c2, r2 * N_ + c1, r2 * N_ + c2
                    };

                    bool hasGiven = false;
                    for (int idx : cells) {
                        if (b_.cells[idx].revealed) { hasGiven = true; break; }
                    }
                    if (!hasGiven) continue;

                    std::set<unsigned int> pairMasks;
                    for (int idx : cells) {
                        if (grid_[idx] == 0 && bits(cand_[idx]) == 2) {
                            pairMasks.insert(cand_[idx]);
                        }
                    }

                    for (unsigned int pairMask : pairMasks) {
                        int exact = 0;
                        int extraCell = -1;
                        bool valid = true;
                        for (int idx : cells) {
                            if (grid_[idx] != 0) {
                                const int v = grid_[idx];
                                if ((pairMask & bit(v)) == 0U) { valid = false; break; }
                                ++exact;
                                continue;
                            }
                            const unsigned int m = cand_[idx];
                            if ((m & pairMask) != pairMask) { valid = false; break; }
                            if (m == pairMask) {
                                ++exact;
                            } else if ((m & ~pairMask) != 0U && extraCell == -1) {
                                extraCell = idx;
                            } else {
                                valid = false;
                                break;
                            }
                        }
                        if (!valid || exact != 3 || extraCell < 0) continue;

                        int local = 0;
                        unsigned int rm = pairMask;
                        while (rm) {
                            const unsigned int one = rm & (~rm + 1U);
                            bool changed = false;
                            if (!removeCandidate(extraCell, firstDigit(one), changed)) return false;
                            if (changed) ++local;
                            rm &= (rm - 1U);
                        }
                        if (local > 0) { n = local; return true; }
                    }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyXChain(int& n) {
    n = 0;
    for (int d = 1; d <= N_; ++d) {
        std::vector<int> nodes;
        for (int i = 0; i < NN_; ++i) {
            if (grid_[i] == 0 && (cand_[i] & bit(d))) nodes.push_back(i);
        }
        const int M = static_cast<int>(nodes.size());
        if (M < 4) continue;

        std::map<int, int> pos;
        for (int i = 0; i < M; ++i) pos[nodes[i]] = i;

        std::vector<std::vector<unsigned char>> edge(M, std::vector<unsigned char>(M, 0));
        for (const auto& h : houses_) {
            std::vector<int> hnodes;
            for (int idx : h) {
                const auto it = pos.find(idx);
                if (it != pos.end()) hnodes.push_back(it->second);
            }
            if (hnodes.size() < 2U) continue;
            const bool strongHouse = (hnodes.size() == 2U);
            for (std::size_t i = 0; i < hnodes.size(); ++i) {
                for (std::size_t j = i + 1; j < hnodes.size(); ++j) {
                    const int a = hnodes[i], b = hnodes[j];
                    if (strongHouse) {
                        edge[a][b] = edge[b][a] = 2U;
                    } else if ((edge[a][b] & 2U) == 0U) {
                        edge[a][b] = edge[b][a] = static_cast<unsigned char>(edge[a][b] | 1U);
                    }
                }
            }
        }

        const int maxDepth = 9;
        for (int s = 0; s < M; ++s) {
            std::vector<char> used(M, 0);
            used[s] = 1;
            std::function<bool(int, bool, int)> dfs = [&](int u, bool expectStrong, int depth) -> bool {
                if (depth >= maxDepth) return false;
                for (int v = 0; v < M; ++v) {
                    if (used[v] || v == u) continue;
                    const unsigned char t = edge[u][v];
                    const bool ok = expectStrong ? ((t & 2U) != 0U) : ((t & 1U) != 0U);
                    if (!ok) continue;

                    const int newDepth = depth + 1;
                    const bool edgeWasStrong = expectStrong;
                    if (edgeWasStrong && newDepth >= 3) {
                        int local = 0;
                        for (int i = 0; i < NN_; ++i) {
                            if (i == nodes[s] || i == nodes[v]) continue;
                            if (grid_[i] != 0 || (cand_[i] & bit(d)) == 0U) continue;
                            if (!isPeerCell(i, nodes[s]) || !isPeerCell(i, nodes[v])) continue;
                            bool changed = false;
                            if (!removeCandidate(i, d, changed)) return true;
                            if (changed) ++local;
                        }
                        if (local > 0) { n = local; return true; }
                    }

                    used[v] = 1;
                    if (dfs(v, !expectStrong, newDepth)) return true;
                    used[v] = 0;
                }
                return false;
            };

            if (dfs(s, true, 0)) return true;
        }
    }
    return false;
}

bool SudokuAnalyzer::applyXYChain(int& n) {
    n = 0;
    std::vector<int> cells;
    for (int i = 0; i < NN_; ++i) {
        if (grid_[i] == 0 && bits(cand_[i]) == 2) cells.push_back(i);
    }
    if (cells.size() < 3U) return false;

    auto otherDigit = [&](unsigned int m, int d) -> int {
        const unsigned int rest = m & ~bit(d);
        return firstDigit(rest);
    };

    const int maxLen = 10;
    for (std::size_t si = 0; si < cells.size(); ++si) {
        const int sCell = cells[si];
        const std::vector<int> sd = digitsFromMask(cand_[sCell]);
        if (sd.size() != 2U) continue;

        for (int z : sd) {
            const int firstShared = (sd[0] == z) ? sd[1] : sd[0];
            std::vector<char> used(cells.size(), 0);
            used[si] = 1;

            std::function<bool(int, int, int)> dfs = [&](int curIdx, int sharedDigit, int len) -> bool {
                if (len >= maxLen) return false;
                const int curCell = cells[curIdx];
                for (std::size_t ni = 0; ni < cells.size(); ++ni) {
                    if (used[ni] || ni == static_cast<std::size_t>(curIdx)) continue;
                    const int nxtCell = cells[ni];
                    if (!isPeerCell(curCell, nxtCell)) continue;
                    if ((cand_[nxtCell] & bit(sharedDigit)) == 0U) continue;

                    const int nextShared = otherDigit(cand_[nxtCell], sharedDigit);
                    if (nextShared == 0) continue;
                    const int newLen = len + 1;

                    if (newLen >= 3 && (cand_[nxtCell] & bit(z))) {
                        int local = 0;
                        for (int i = 0; i < NN_; ++i) {
                            if (i == sCell || i == nxtCell) continue;
                            if (grid_[i] != 0 || (cand_[i] & bit(z)) == 0U) continue;
                            if (!isPeerCell(i, sCell) || !isPeerCell(i, nxtCell)) continue;
                            bool changed = false;
                            if (!removeCandidate(i, z, changed)) return true;
                            if (changed) ++local;
                        }
                        if (local > 0) { n = local; return true; }
                    }

                    used[ni] = 1;
                    if (dfs(static_cast<int>(ni), nextShared, newLen)) return true;
                    used[ni] = 0;
                }
                return false;
            };

            if (dfs(static_cast<int>(si), firstShared, 1)) return true;
        }
    }
    return false;
}

bool SudokuAnalyzer::applyAIC(int& n) {
    n = 0;

    struct Node {
        int cell = -1;
        int digit = 0;
    };

    std::vector<Node> nodes;
    std::vector<int> nodeByCellDigit(NN_ * (N_ + 1), -1);
    auto nodeKey = [&](int cell, int digit) { return cell * (N_ + 1) + digit; };

    for (int cell = 0; cell < NN_; ++cell) {
        if (grid_[cell] != 0) continue;
        unsigned int m = cand_[cell];
        while (m) {
            const unsigned int one = m & (~m + 1U);
            const int d = firstDigit(one);
            const int id = static_cast<int>(nodes.size());
            nodes.push_back({cell, d});
            nodeByCellDigit[nodeKey(cell, d)] = id;
            m &= (m - 1U);
        }
    }
    if (nodes.size() < 4U) return false;

    std::vector<std::vector<std::pair<int, unsigned char>>> adj(nodes.size());
    auto addEdge = [&](int a, int b, unsigned char edgeType) {
        if (a < 0 || b < 0 || a == b) return;
        bool found = false;
        for (auto& pr : adj[a]) {
            if (pr.first == b) {
                pr.second = static_cast<unsigned char>(pr.second | edgeType);
                found = true;
                break;
            }
        }
        if (!found) adj[a].push_back({b, edgeType});
        found = false;
        for (auto& pr : adj[b]) {
            if (pr.first == a) {
                pr.second = static_cast<unsigned char>(pr.second | edgeType);
                found = true;
                break;
            }
        }
        if (!found) adj[b].push_back({a, edgeType});
    };

    for (int cell = 0; cell < NN_; ++cell) {
        if (grid_[cell] != 0) continue;
        std::vector<int> ds = digitsFromMask(cand_[cell]);
        if (ds.size() < 2U) continue;
        const unsigned char t = (ds.size() == 2U) ? static_cast<unsigned char>(3U) : static_cast<unsigned char>(1U);
        for (std::size_t i = 0; i < ds.size(); ++i) {
            for (std::size_t j = i + 1; j < ds.size(); ++j) {
                const int a = nodeByCellDigit[nodeKey(cell, ds[i])];
                const int b = nodeByCellDigit[nodeKey(cell, ds[j])];
                addEdge(a, b, t);
            }
        }
    }

    for (const auto& h : houses_) {
        for (int d = 1; d <= N_; ++d) {
            std::vector<int> houseNodes;
            for (int idx : h) {
                const int nid = nodeByCellDigit[nodeKey(idx, d)];
                if (nid >= 0) houseNodes.push_back(nid);
            }
            if (houseNodes.size() < 2U) continue;
            const unsigned char t = (houseNodes.size() == 2U) ? static_cast<unsigned char>(3U) : static_cast<unsigned char>(1U);
            for (std::size_t i = 0; i < houseNodes.size(); ++i) {
                for (std::size_t j = i + 1; j < houseNodes.size(); ++j) {
                    addEdge(houseNodes[i], houseNodes[j], t);
                }
            }
        }
    }

    const int maxDepth = 12;
    for (int s = 0; s < static_cast<int>(nodes.size()); ++s) {
        std::vector<char> used(nodes.size(), 0);
        used[s] = 1;

        std::function<bool(int, bool, int)> dfs = [&](int u, bool expectStrong, int depth) -> bool {
            if (depth >= maxDepth) return false;
            for (const auto& e : adj[u]) {
                const int v = e.first;
                const unsigned char et = e.second;
                if (used[v]) continue;
                if (expectStrong && (et & 2U) == 0U) continue;
                if (!expectStrong && (et & 1U) == 0U) continue;

                const int newDepth = depth + 1;
                const bool edgeWasStrong = expectStrong;
                if (edgeWasStrong && newDepth >= 3 &&
                    nodes[v].digit == nodes[s].digit &&
                    nodes[v].cell != nodes[s].cell) {
                    const int d = nodes[s].digit;
                    int local = 0;
                    for (int i = 0; i < NN_; ++i) {
                        if (i == nodes[s].cell || i == nodes[v].cell) continue;
                        if (grid_[i] != 0 || (cand_[i] & bit(d)) == 0U) continue;
                        if (!isPeerCell(i, nodes[s].cell) || !isPeerCell(i, nodes[v].cell)) continue;
                        bool changed = false;
                        if (!removeCandidate(i, d, changed)) return false;
                        if (changed) ++local;
                    }
                    if (local > 0) {
                        n = local;
                        return true;
                    }
                }

                used[v] = 1;
                if (dfs(v, !expectStrong, newDepth)) return true;
                used[v] = 0;
            }
            return false;
        };

        if (dfs(s, true, 0)) return true;
    }

    return false;
}

bool SudokuAnalyzer::applyContinuousNiceLoop(int& n) {
    n = 0;

    struct Node {
        int cell = -1;
        int digit = 0;
    };

    std::vector<Node> nodes;
    std::vector<int> nodeByCellDigit(NN_ * (N_ + 1), -1);
    auto nodeKey = [&](int cell, int digit) { return cell * (N_ + 1) + digit; };

    for (int cell = 0; cell < NN_; ++cell) {
        if (grid_[cell] != 0) continue;
        unsigned int m = cand_[cell];
        while (m) {
            const unsigned int one = m & (~m + 1U);
            const int d = firstDigit(one);
            const int id = static_cast<int>(nodes.size());
            nodes.push_back({cell, d});
            nodeByCellDigit[nodeKey(cell, d)] = id;
            m &= (m - 1U);
        }
    }
    if (nodes.size() < 4U) return false;

    const int M = static_cast<int>(nodes.size());
    std::vector<std::vector<unsigned char>> edge(M, std::vector<unsigned char>(M, 0U));
    auto addEdge = [&](int a, int b, unsigned char edgeType) {
        if (a < 0 || b < 0 || a == b) return;
        edge[a][b] = static_cast<unsigned char>(edge[a][b] | edgeType);
        edge[b][a] = static_cast<unsigned char>(edge[b][a] | edgeType);
    };

    for (int cell = 0; cell < NN_; ++cell) {
        if (grid_[cell] != 0) continue;
        const std::vector<int> ds = digitsFromMask(cand_[cell]);
        if (ds.size() < 2U) continue;
        const unsigned char t = (ds.size() == 2U) ? static_cast<unsigned char>(3U) : static_cast<unsigned char>(1U);
        for (std::size_t i = 0; i < ds.size(); ++i) {
            for (std::size_t j = i + 1; j < ds.size(); ++j) {
                addEdge(nodeByCellDigit[nodeKey(cell, ds[i])], nodeByCellDigit[nodeKey(cell, ds[j])], t);
            }
        }
    }
    for (const auto& h : houses_) {
        for (int d = 1; d <= N_; ++d) {
            std::vector<int> houseNodes;
            for (int idx : h) {
                const int nid = nodeByCellDigit[nodeKey(idx, d)];
                if (nid >= 0) houseNodes.push_back(nid);
            }
            if (houseNodes.size() < 2U) continue;
            const unsigned char t = (houseNodes.size() == 2U) ? static_cast<unsigned char>(3U) : static_cast<unsigned char>(1U);
            for (std::size_t i = 0; i < houseNodes.size(); ++i) {
                for (std::size_t j = i + 1; j < houseNodes.size(); ++j) {
                    addEdge(houseNodes[i], houseNodes[j], t);
                }
            }
        }
    }

    auto tryWeakEdgeElim = [&](int a, int b) -> bool {
        if (nodes[a].digit != nodes[b].digit) return false;
        if (nodes[a].cell == nodes[b].cell) return false;
        const int d = nodes[a].digit;
        int local = 0;
        for (int i = 0; i < NN_; ++i) {
            if (i == nodes[a].cell || i == nodes[b].cell) continue;
            if (grid_[i] != 0 || (cand_[i] & bit(d)) == 0U) continue;
            if (!isPeerCell(i, nodes[a].cell) || !isPeerCell(i, nodes[b].cell)) continue;
            bool changed = false;
            if (!removeCandidate(i, d, changed)) return false;
            if (changed) ++local;
        }
        if (local > 0) {
            n = local;
            std::ostringstream ss;
            ss << "ContinuousNiceLoop: weak-link closure removes " << d
               << " from peers of " << cellName(nodes[a].cell)
               << " and " << cellName(nodes[b].cell);
            pushDebugLog(ss.str());
            return true;
        }
        return false;
    };

    for (int u = 0; u < M; ++u) {
        for (int v = 0; v < M; ++v) {
            if (v == u || (edge[u][v] & 2U) == 0U) continue; // strong
            for (int w = 0; w < M; ++w) {
                if (w == u || w == v || (edge[v][w] & 1U) == 0U) continue; // weak
                for (int x = 0; x < M; ++x) {
                    if (x == u || x == v || x == w) continue;
                    if ((edge[w][x] & 2U) == 0U) continue; // strong
                    if ((edge[x][u] & 1U) == 0U) continue; // weak closes loop

                    if (tryWeakEdgeElim(v, w)) return true;
                    if (tryWeakEdgeElim(x, u)) return true;
                }
            }
        }
    }

    return false;
}

bool SudokuAnalyzer::applySKLoop(int& n) {
    n = 0;
    if (N_ != 9) return false;

    auto strongInHouse = [&](const std::vector<int>& house, int d, int a, int b) {
        std::vector<int> where;
        for (int idx : house) {
            if (grid_[idx] == 0 && (cand_[idx] & bit(d))) where.push_back(idx);
        }
        return where.size() == 2U &&
               ((where[0] == a && where[1] == b) || (where[0] == b && where[1] == a));
    };

    for (int r1 = 0; r1 < N_; ++r1) {
        for (int r2 = r1 + 1; r2 < N_; ++r2) {
            for (int c1 = 0; c1 < N_; ++c1) {
                for (int c2 = c1 + 1; c2 < N_; ++c2) {
                    const int a = r1 * N_ + c1;
                    const int b = r1 * N_ + c2;
                    const int c = r2 * N_ + c1;
                    const int d = r2 * N_ + c2;
                    if (grid_[a] != 0 || grid_[b] != 0 || grid_[c] != 0 || grid_[d] != 0) continue;
                    if (bits(cand_[a]) != 2 || bits(cand_[b]) != 2 || bits(cand_[c]) != 2 || bits(cand_[d]) != 2) continue;

                    const unsigned int ab = cand_[a] & cand_[b];
                    const unsigned int bd = cand_[b] & cand_[d];
                    const unsigned int dc = cand_[d] & cand_[c];
                    const unsigned int ca = cand_[c] & cand_[a];
                    if (bits(ab) != 1 || bits(bd) != 1 || bits(dc) != 1 || bits(ca) != 1) continue;

                    const int dab = firstDigit(ab);
                    const int dbd = firstDigit(bd);
                    const int ddc = firstDigit(dc);
                    const int dca = firstDigit(ca);

                    if ((cand_[a] & ~(bit(dab) | bit(dca))) != 0U) continue;
                    if ((cand_[b] & ~(bit(dab) | bit(dbd))) != 0U) continue;
                    if ((cand_[d] & ~(bit(dbd) | bit(ddc))) != 0U) continue;
                    if ((cand_[c] & ~(bit(ddc) | bit(dca))) != 0U) continue;

                    if (!strongInHouse(houses_[r1], dab, a, b)) continue;
                    if (!strongInHouse(houses_[N_ + c2], dbd, b, d)) continue;
                    if (!strongInHouse(houses_[r2], ddc, d, c)) continue;
                    if (!strongInHouse(houses_[N_ + c1], dca, c, a)) continue;

                    int local = 0;
                    auto elimFromPeers = [&](int u, int v, int dig) -> bool {
                        for (int idx = 0; idx < NN_; ++idx) {
                            if (idx == u || idx == v || grid_[idx] != 0) continue;
                            if ((cand_[idx] & bit(dig)) == 0U) continue;
                            if (!isPeerCell(idx, u) || !isPeerCell(idx, v)) continue;
                            bool changed = false;
                            if (!removeCandidate(idx, dig, changed)) return false;
                            if (changed) ++local;
                        }
                        return true;
                    };
                    if (!elimFromPeers(a, b, dab)) return false;
                    if (!elimFromPeers(b, d, dbd)) return false;
                    if (!elimFromPeers(d, c, ddc)) return false;
                    if (!elimFromPeers(c, a, dca)) return false;

                    if (local > 0) {
                        std::ostringstream ss;
                        ss << "SKLoop: rectangle r" << (r1 + 1) << "/r" << (r2 + 1)
                           << " c" << (c1 + 1) << "/c" << (c2 + 1)
                           << " removed " << local << " candidate(s)";
                        pushDebugLog(ss.str());
                        n = local;
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyALSXZ(int& n) {
    n = 0;

    struct ALS {
        std::vector<int> cells;
        unsigned int mask = 0U;
    };

    std::vector<ALS> alss;
    std::set<std::string> seen;

    auto addAls = [&](std::vector<int> cells) {
        if (cells.empty() || cells.size() > 3U) return;
        std::sort(cells.begin(), cells.end());
        std::ostringstream key;
        for (int c : cells) key << c << ",";
        if (!seen.insert(key.str()).second) return;

        unsigned int m = 0U;
        for (int c : cells) {
            if (grid_[c] != 0) return;
            m |= cand_[c];
        }
        if (bits(m) != static_cast<int>(cells.size()) + 1) return;
        alss.push_back({std::move(cells), m});
    };

    for (int i = 0; i < NN_; ++i) {
        if (grid_[i] == 0 && bits(cand_[i]) == 2) {
            addAls({i});
        }
    }

    for (const auto& h : houses_) {
        std::vector<int> unsolved;
        for (int i : h) if (grid_[i] == 0) unsolved.push_back(i);
        if (unsolved.size() >= 2U) {
            forEachCombo(unsolved, 2, [&](const std::vector<int>& cs) { addAls(cs); });
        }
        if (unsolved.size() >= 3U) {
            forEachCombo(unsolved, 3, [&](const std::vector<int>& cs) { addAls(cs); });
        }
    }

    auto hasCell = [&](const ALS& a, int c) {
        return std::find(a.cells.begin(), a.cells.end(), c) != a.cells.end();
    };
    auto seesAll = [&](int c, const std::vector<int>& xs) {
        for (int x : xs) if (!isPeerCell(c, x)) return false;
        return true;
    };

    for (std::size_t i = 0; i < alss.size(); ++i) {
        for (std::size_t j = i + 1; j < alss.size(); ++j) {
            const ALS& A = alss[i];
            const ALS& B = alss[j];

            bool disjoint = true;
            for (int c : A.cells) {
                if (hasCell(B, c)) { disjoint = false; break; }
            }
            if (!disjoint) continue;

            const unsigned int common = A.mask & B.mask;
            if (bits(common) < 2) continue;

            for (int x : digitsFromMask(common)) {
                std::vector<int> ax, bx;
                for (int c : A.cells) if (cand_[c] & bit(x)) ax.push_back(c);
                for (int c : B.cells) if (cand_[c] & bit(x)) bx.push_back(c);
                if (ax.empty() || bx.empty()) continue;

                bool restricted = true;
                for (int ca : ax) {
                    for (int cb : bx) {
                        if (!isPeerCell(ca, cb)) {
                            restricted = false;
                            break;
                        }
                    }
                    if (!restricted) break;
                }
                if (!restricted) continue;

                unsigned int zMask = common & ~bit(x);
                while (zMask) {
                    const unsigned int one = zMask & (~zMask + 1U);
                    const int z = firstDigit(one);
                    zMask &= (zMask - 1U);

                    std::vector<int> az, bz;
                    for (int c : A.cells) if (cand_[c] & bit(z)) az.push_back(c);
                    for (int c : B.cells) if (cand_[c] & bit(z)) bz.push_back(c);
                    if (az.empty() || bz.empty()) continue;

                    int local = 0;
                    for (int c = 0; c < NN_; ++c) {
                        if (grid_[c] != 0 || (cand_[c] & bit(z)) == 0U) continue;
                        if (hasCell(A, c) || hasCell(B, c)) continue;
                        if (!seesAll(c, az) || !seesAll(c, bz)) continue;
                        bool changed = false;
                        if (!removeCandidate(c, z, changed)) return false;
                        if (changed) ++local;
                    }
                    if (local > 0) {
                        n = local;
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

bool SudokuAnalyzer::applyALSXYWing(int& n) {
    n = 0;

    struct ALS {
        std::vector<int> cells;
        unsigned int mask = 0U;
    };
    std::vector<ALS> alss;
    std::set<std::string> seen;
    auto addAls = [&](std::vector<int> cells) {
        if (cells.empty() || cells.size() > 3U) return;
        std::sort(cells.begin(), cells.end());
        std::ostringstream key;
        for (int c : cells) key << c << ",";
        if (!seen.insert(key.str()).second) return;
        unsigned int m = 0U;
        for (int c : cells) {
            if (grid_[c] != 0) return;
            m |= cand_[c];
        }
        if (bits(m) != static_cast<int>(cells.size()) + 1) return;
        alss.push_back({std::move(cells), m});
    };
    for (int i = 0; i < NN_; ++i) if (grid_[i] == 0 && bits(cand_[i]) == 2) addAls({i});
    for (const auto& h : houses_) {
        std::vector<int> uns;
        for (int i : h) if (grid_[i] == 0) uns.push_back(i);
        if (uns.size() >= 2U) forEachCombo(uns, 2, [&](const std::vector<int>& cs) { addAls(cs); });
        if (uns.size() >= 3U) forEachCombo(uns, 3, [&](const std::vector<int>& cs) { addAls(cs); });
    }
    if (alss.size() > 260U) alss.resize(260U);

    auto hasCell = [&](const ALS& a, int c) {
        return std::find(a.cells.begin(), a.cells.end(), c) != a.cells.end();
    };
    auto disjoint3 = [&](const ALS& A, const ALS& B, const ALS& C) {
        for (int c : A.cells) if (hasCell(B, c) || hasCell(C, c)) return false;
        for (int c : B.cells) if (hasCell(C, c)) return false;
        return true;
    };
    auto digitCells = [&](const ALS& a, int d) {
        std::vector<int> out;
        for (int c : a.cells) if (cand_[c] & bit(d)) out.push_back(c);
        return out;
    };
    auto isRCC = [&](const ALS& A, const ALS& B, int d) {
        const std::vector<int> da = digitCells(A, d);
        const std::vector<int> db = digitCells(B, d);
        if (da.empty() || db.empty()) return false;
        for (int ca : da) for (int cb : db) if (!isPeerCell(ca, cb)) return false;
        return true;
    };
    auto seesAll = [&](int cell, const std::vector<int>& xs) {
        for (int x : xs) if (!isPeerCell(cell, x)) return false;
        return true;
    };

    for (std::size_t i = 0; i < alss.size(); ++i) {
        for (std::size_t j = 0; j < alss.size(); ++j) {
            if (j == i) continue;
            for (std::size_t k = 0; k < alss.size(); ++k) {
                if (k == i || k == j) continue;
                const ALS& A = alss[i];
                const ALS& B = alss[j];
                const ALS& C = alss[k];
                if (!disjoint3(A, B, C)) continue;

                const unsigned int ab = A.mask & B.mask;
                const unsigned int ac = A.mask & C.mask;
                const unsigned int bc = B.mask & C.mask;
                if (ab == 0U || ac == 0U || bc == 0U) continue;

                for (int x : digitsFromMask(ab)) {
                    if (!isRCC(A, B, x)) continue;
                    for (int y : digitsFromMask(ac)) {
                        if (y == x || !isRCC(A, C, y)) continue;
                        unsigned int zMask = bc & ~bit(x) & ~bit(y);
                        while (zMask) {
                            const unsigned int one = zMask & (~zMask + 1U);
                            const int z = firstDigit(one);
                            zMask &= (zMask - 1U);
                            const std::vector<int> bz = digitCells(B, z);
                            const std::vector<int> cz = digitCells(C, z);
                            if (bz.empty() || cz.empty()) continue;

                            int local = 0;
                            for (int cell = 0; cell < NN_; ++cell) {
                                if (grid_[cell] != 0 || (cand_[cell] & bit(z)) == 0U) continue;
                                if (hasCell(A, cell) || hasCell(B, cell) || hasCell(C, cell)) continue;
                                if (!seesAll(cell, bz) || !seesAll(cell, cz)) continue;
                                bool changed = false;
                                if (!removeCandidate(cell, z, changed)) return false;
                                if (changed) ++local;
                            }
                            if (local > 0) {
                                std::ostringstream ss;
                                ss << "ALS-XY-Wing: remove " << z << " from " << local << " cell(s)";
                                pushDebugLog(ss.str());
                                n = local;
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyALSChain(int& n) {
    n = 0;

    struct ALS {
        std::vector<int> cells;
        unsigned int mask = 0U;
    };
    std::vector<ALS> alss;
    std::set<std::string> seen;
    auto addAls = [&](std::vector<int> cells) {
        if (cells.empty() || cells.size() > 3U) return;
        std::sort(cells.begin(), cells.end());
        std::ostringstream key;
        for (int c : cells) key << c << ",";
        if (!seen.insert(key.str()).second) return;
        unsigned int m = 0U;
        for (int c : cells) {
            if (grid_[c] != 0) return;
            m |= cand_[c];
        }
        if (bits(m) != static_cast<int>(cells.size()) + 1) return;
        alss.push_back({std::move(cells), m});
    };
    for (int i = 0; i < NN_; ++i) if (grid_[i] == 0 && bits(cand_[i]) == 2) addAls({i});
    for (const auto& h : houses_) {
        std::vector<int> uns;
        for (int i : h) if (grid_[i] == 0) uns.push_back(i);
        if (uns.size() >= 2U) forEachCombo(uns, 2, [&](const std::vector<int>& cs) { addAls(cs); });
        if (uns.size() >= 3U) forEachCombo(uns, 3, [&](const std::vector<int>& cs) { addAls(cs); });
    }
    if (alss.size() > 260U) alss.resize(260U);

    auto hasCell = [&](const ALS& a, int c) {
        return std::find(a.cells.begin(), a.cells.end(), c) != a.cells.end();
    };
    auto disjoint3 = [&](const ALS& A, const ALS& B, const ALS& C) {
        for (int c : A.cells) if (hasCell(B, c) || hasCell(C, c)) return false;
        for (int c : B.cells) if (hasCell(C, c)) return false;
        return true;
    };
    auto digitCells = [&](const ALS& a, int d) {
        std::vector<int> out;
        for (int c : a.cells) if (cand_[c] & bit(d)) out.push_back(c);
        return out;
    };
    auto isRCC = [&](const ALS& A, const ALS& B, int d) {
        const std::vector<int> da = digitCells(A, d);
        const std::vector<int> db = digitCells(B, d);
        if (da.empty() || db.empty()) return false;
        for (int ca : da) for (int cb : db) if (!isPeerCell(ca, cb)) return false;
        return true;
    };
    auto seesAll = [&](int cell, const std::vector<int>& xs) {
        for (int x : xs) if (!isPeerCell(cell, x)) return false;
        return true;
    };

    for (std::size_t i = 0; i < alss.size(); ++i) {
        for (std::size_t j = 0; j < alss.size(); ++j) {
            if (j == i) continue;
            for (std::size_t k = 0; k < alss.size(); ++k) {
                if (k == i || k == j) continue;
                const ALS& A = alss[i];
                const ALS& B = alss[j];
                const ALS& C = alss[k];
                if (!disjoint3(A, B, C)) continue;

                const unsigned int ab = A.mask & B.mask;
                const unsigned int bc = B.mask & C.mask;
                const unsigned int ac = A.mask & C.mask;
                if (ab == 0U || bc == 0U || ac == 0U) continue;

                for (int x : digitsFromMask(ab)) {
                    if (!isRCC(A, B, x)) continue;
                    for (int y : digitsFromMask(bc)) {
                        if (y == x || !isRCC(B, C, y)) continue;
                        unsigned int zMask = ac & ~bit(x) & ~bit(y);
                        while (zMask) {
                            const unsigned int one = zMask & (~zMask + 1U);
                            const int z = firstDigit(one);
                            zMask &= (zMask - 1U);
                            const std::vector<int> az = digitCells(A, z);
                            const std::vector<int> cz = digitCells(C, z);
                            if (az.empty() || cz.empty()) continue;

                            int local = 0;
                            for (int cell = 0; cell < NN_; ++cell) {
                                if (grid_[cell] != 0 || (cand_[cell] & bit(z)) == 0U) continue;
                                if (hasCell(A, cell) || hasCell(B, cell) || hasCell(C, cell)) continue;
                                if (!seesAll(cell, az) || !seesAll(cell, cz)) continue;
                                bool changed = false;
                                if (!removeCandidate(cell, z, changed)) return false;
                                if (changed) ++local;
                            }
                            if (local > 0) {
                                std::ostringstream ss;
                                ss << "ALS-Chain: remove " << z << " from " << local << " cell(s)";
                                pushDebugLog(ss.str());
                                n = local;
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyDeathBlossom(int& n) {
    n = 0;
    if (N_ != 9) return false;

    struct ALS {
        std::vector<int> cells;
        unsigned int mask = 0U;
    };
    std::vector<ALS> alss;
    std::set<std::string> seen;
    auto addAls = [&](std::vector<int> cells) {
        if (cells.empty() || cells.size() > 3U) return;
        std::sort(cells.begin(), cells.end());
        std::ostringstream key;
        for (int c : cells) key << c << ",";
        if (!seen.insert(key.str()).second) return;
        unsigned int m = 0U;
        for (int c : cells) {
            if (grid_[c] != 0) return;
            m |= cand_[c];
        }
        if (bits(m) != static_cast<int>(cells.size()) + 1) return;
        alss.push_back({std::move(cells), m});
    };
    for (int i = 0; i < NN_; ++i) if (grid_[i] == 0 && bits(cand_[i]) == 2) addAls({i});
    for (const auto& h : houses_) {
        std::vector<int> uns;
        for (int i : h) if (grid_[i] == 0) uns.push_back(i);
        if (uns.size() >= 2U) forEachCombo(uns, 2, [&](const std::vector<int>& cs) { addAls(cs); });
        if (uns.size() >= 3U) forEachCombo(uns, 3, [&](const std::vector<int>& cs) { addAls(cs); });
    }
    if (alss.size() > 220U) alss.resize(220U);

    auto hasCell = [&](const ALS& a, int c) {
        return std::find(a.cells.begin(), a.cells.end(), c) != a.cells.end();
    };
    auto digitCells = [&](const ALS& a, int d) {
        std::vector<int> out;
        for (int c : a.cells) if (cand_[c] & bit(d)) out.push_back(c);
        return out;
    };
    auto rccPivot = [&](int pivot, const ALS& a, int d) {
        const std::vector<int> ds = digitCells(a, d);
        if (ds.empty()) return false;
        for (int c : ds) if (!isPeerCell(pivot, c)) return false;
        return true;
    };
    auto seesAll = [&](int cell, const std::vector<int>& xs) {
        for (int x : xs) if (!isPeerCell(cell, x)) return false;
        return true;
    };

    for (int pivot = 0; pivot < NN_; ++pivot) {
        if (grid_[pivot] != 0 || bits(cand_[pivot]) != 3) continue;
        const std::vector<int> pd = digitsFromMask(cand_[pivot]);
        if (pd.size() != 3U) continue;

        for (int z = 1; z <= N_; ++z) {
            if ((cand_[pivot] & bit(z)) != 0U) continue;

            std::vector<int> p0, p1, p2;
            for (int ai = 0; ai < static_cast<int>(alss.size()); ++ai) {
                const ALS& A = alss[ai];
                if (hasCell(A, pivot)) continue;
                if ((A.mask & bit(z)) == 0U) continue;
                if ((A.mask & bit(pd[0])) && rccPivot(pivot, A, pd[0])) p0.push_back(ai);
                if ((A.mask & bit(pd[1])) && rccPivot(pivot, A, pd[1])) p1.push_back(ai);
                if ((A.mask & bit(pd[2])) && rccPivot(pivot, A, pd[2])) p2.push_back(ai);
            }
            if (p0.empty() || p1.empty() || p2.empty()) continue;

            for (int a : p0) for (int b : p1) for (int c : p2) {
                if (a == b || b == c || a == c) continue;
                const ALS& A0 = alss[a];
                const ALS& A1 = alss[b];
                const ALS& A2 = alss[c];

                bool disjoint = true;
                for (int v : A0.cells) if (hasCell(A1, v) || hasCell(A2, v)) { disjoint = false; break; }
                if (!disjoint) continue;
                for (int v : A1.cells) if (hasCell(A2, v)) { disjoint = false; break; }
                if (!disjoint) continue;

                const std::vector<int> z0 = digitCells(A0, z);
                const std::vector<int> z1 = digitCells(A1, z);
                const std::vector<int> z2 = digitCells(A2, z);
                if (z0.empty() || z1.empty() || z2.empty()) continue;

                int local = 0;
                for (int cell = 0; cell < NN_; ++cell) {
                    if (grid_[cell] != 0 || (cand_[cell] & bit(z)) == 0U) continue;
                    if (cell == pivot || hasCell(A0, cell) || hasCell(A1, cell) || hasCell(A2, cell)) continue;
                    if (!seesAll(cell, z0) || !seesAll(cell, z1) || !seesAll(cell, z2)) continue;
                    bool changed = false;
                    if (!removeCandidate(cell, z, changed)) return false;
                    if (changed) ++local;
                }
                if (local > 0) {
                    std::ostringstream ss;
                    ss << "DeathBlossom: pivot " << cellName(pivot) << " remove " << z
                       << " from " << local << " cell(s)";
                    pushDebugLog(ss.str());
                    n = local;
                    return true;
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applySueDeCoq(int& n) {
    n = 0;

    auto tryMode = [&](bool rowMode) -> bool {
        for (int b = 0; b < N_; ++b) {
            const int bpr = N_ / BC_;
            const int sr = (b / bpr) * BR_;
            const int sc = (b % bpr) * BC_;

            for (int off = 0; off < (rowMode ? BR_ : BC_); ++off) {
                const int line = (rowMode ? (sr + off) : (sc + off));
                std::vector<int> I;
                for (int d1 = 0; d1 < (rowMode ? BC_ : BR_); ++d1) {
                    const int r = rowMode ? line : (sr + d1);
                    const int c = rowMode ? (sc + d1) : line;
                    const int idx = r * N_ + c;
                    if (grid_[idx] == 0) I.push_back(idx);
                }
                if (I.size() != 2U) continue;

                unsigned int mI = 0U;
                for (int c : I) mI |= cand_[c];

                std::vector<int> linePool;
                const std::vector<int>& lineHouse = rowMode ? houses_[line] : houses_[N_ + line];
                for (int idx : lineHouse) {
                    if (grid_[idx] != 0) continue;
                    if (box(row(idx), col(idx)) == b) continue;
                    linePool.push_back(idx);
                }
                std::vector<int> boxPool;
                const std::vector<int>& boxHouse = houses_[2 * N_ + b];
                for (int idx : boxHouse) {
                    if (grid_[idx] != 0) continue;
                    if ((rowMode && row(idx) == line) || (!rowMode && col(idx) == line)) continue;
                    boxPool.push_back(idx);
                }

                if (linePool.empty() || boxPool.empty()) continue;

                for (int kL = 1; kL <= std::min(2, static_cast<int>(linePool.size())); ++kL) {
                    for (int kB = 1; kB <= std::min(2, static_cast<int>(boxPool.size())); ++kB) {
                        bool found = false;
                        forEachCombo(linePool, kL, [&](const std::vector<int>& Lset) {
                            if (found || contradiction_) return;
                            unsigned int mL = 0U;
                            for (int c : Lset) mL |= cand_[c];

                            forEachCombo(boxPool, kB, [&](const std::vector<int>& Bset) {
                                if (found || contradiction_) return;
                                unsigned int mB = 0U;
                                for (int c : Bset) mB |= cand_[c];

                                if ((mL & ~mI) == 0U || (mB & ~mI) == 0U) return;
                                if (((mL & mB) & ~mI) != 0U) return;

                                const unsigned int mIL = mI | mL;
                                const unsigned int mIB = mI | mB;
                                const unsigned int mAll = mI | mL | mB;
                                if (bits(mIL) != static_cast<int>(I.size()) + static_cast<int>(Lset.size())) return;
                                if (bits(mIB) != static_cast<int>(I.size()) + static_cast<int>(Bset.size())) return;
                                if (bits(mAll) != static_cast<int>(I.size()) + static_cast<int>(Lset.size()) + static_cast<int>(Bset.size())) return;

                                std::set<int> Iset(I.begin(), I.end());
                                std::set<int> Lmark(Lset.begin(), Lset.end());
                                std::set<int> Bmark(Bset.begin(), Bset.end());

                                int local = 0;
                                for (int idx : lineHouse) {
                                    if (grid_[idx] != 0) continue;
                                    if (Iset.count(idx) || Lmark.count(idx)) continue;
                                    unsigned int rm = cand_[idx] & mIL;
                                    while (rm) {
                                        const unsigned int one = rm & (~rm + 1U);
                                        bool changed = false;
                                        if (!removeCandidate(idx, firstDigit(one), changed)) return;
                                        if (changed) ++local;
                                        rm &= (rm - 1U);
                                    }
                                }
                                for (int idx : boxHouse) {
                                    if (grid_[idx] != 0) continue;
                                    if (Iset.count(idx) || Bmark.count(idx)) continue;
                                    unsigned int rm = cand_[idx] & mIB;
                                    while (rm) {
                                        const unsigned int one = rm & (~rm + 1U);
                                        bool changed = false;
                                        if (!removeCandidate(idx, firstDigit(one), changed)) return;
                                        if (changed) ++local;
                                        rm &= (rm - 1U);
                                    }
                                }
                                if (local > 0) {
                                    n = local;
                                    found = true;
                                }
                            });
                        });
                        if (found) return true;
                    }
                }
            }
        }
        return false;
    };

    if (tryMode(true)) return true;
    return tryMode(false);
}

bool SudokuAnalyzer::applyMSLS(int& n) {
    n = 0;
    if (N_ != 9) return false;

    auto tryMode = [&](bool rowMode) -> bool {
        for (int b = 0; b < N_; ++b) {
            const int bpr = N_ / BC_;
            const int sr = (b / bpr) * BR_;
            const int sc = (b % bpr) * BC_;

            for (int off = 0; off < (rowMode ? BR_ : BC_); ++off) {
                const int line = rowMode ? (sr + off) : (sc + off);

                std::vector<int> inter;
                for (int k = 0; k < (rowMode ? BC_ : BR_); ++k) {
                    const int r = rowMode ? line : (sr + k);
                    const int c = rowMode ? (sc + k) : line;
                    const int idx = r * N_ + c;
                    if (grid_[idx] == 0) inter.push_back(idx);
                }
                if (inter.size() < 2U) continue;

                const std::vector<int>& lineHouse = rowMode ? houses_[line] : houses_[N_ + line];
                const std::vector<int>& boxHouse = houses_[2 * N_ + b];

                std::vector<int> linePool;
                for (int idx : lineHouse) {
                    if (grid_[idx] != 0) continue;
                    if (box(row(idx), col(idx)) == b) continue;
                    linePool.push_back(idx);
                }
                std::vector<int> boxPool;
                for (int idx : boxHouse) {
                    if (grid_[idx] != 0) continue;
                    if ((rowMode && row(idx) == line) || (!rowMode && col(idx) == line)) continue;
                    boxPool.push_back(idx);
                }
                if (linePool.empty() || boxPool.empty()) continue;

                const int iMax = std::min(3, static_cast<int>(inter.size()));
                for (int isz = 2; isz <= iMax; ++isz) {
                    bool found = false;
                    forEachCombo(inter, isz, [&](const std::vector<int>& Iset) {
                        if (found || contradiction_) return;
                        unsigned int mI = 0U;
                        for (int c : Iset) mI |= cand_[c];

                        for (int kL = 1; kL <= std::min(3, static_cast<int>(linePool.size())); ++kL) {
                            forEachCombo(linePool, kL, [&](const std::vector<int>& Lset) {
                                if (found || contradiction_) return;
                                unsigned int mL = 0U;
                                for (int c : Lset) mL |= cand_[c];

                                for (int kB = 1; kB <= std::min(3, static_cast<int>(boxPool.size())); ++kB) {
                                    forEachCombo(boxPool, kB, [&](const std::vector<int>& Bset) {
                                        if (found || contradiction_) return;
                                        unsigned int mB = 0U;
                                        for (int c : Bset) mB |= cand_[c];

                                        if ((mL & ~mI) == 0U || (mB & ~mI) == 0U) return;
                                        if (((mL & mB) & ~mI) != 0U) return;

                                        const unsigned int mIL = mI | mL;
                                        const unsigned int mIB = mI | mB;
                                        const unsigned int mAll = mI | mL | mB;
                                        if (bits(mIL) != static_cast<int>(Iset.size()) + static_cast<int>(Lset.size())) return;
                                        if (bits(mIB) != static_cast<int>(Iset.size()) + static_cast<int>(Bset.size())) return;
                                        if (bits(mAll) != static_cast<int>(Iset.size()) + static_cast<int>(Lset.size()) + static_cast<int>(Bset.size())) return;

                                        std::set<int> Imark(Iset.begin(), Iset.end());
                                        std::set<int> Lmark(Lset.begin(), Lset.end());
                                        std::set<int> Bmark(Bset.begin(), Bset.end());

                                        int local = 0;
                                        for (int idx : lineHouse) {
                                            if (grid_[idx] != 0) continue;
                                            if (Imark.count(idx) || Lmark.count(idx)) continue;
                                            unsigned int rm = cand_[idx] & mIL;
                                            while (rm) {
                                                const unsigned int one = rm & (~rm + 1U);
                                                bool changed = false;
                                                if (!removeCandidate(idx, firstDigit(one), changed)) return;
                                                if (changed) ++local;
                                                rm &= (rm - 1U);
                                            }
                                        }
                                        for (int idx : boxHouse) {
                                            if (grid_[idx] != 0) continue;
                                            if (Imark.count(idx) || Bmark.count(idx)) continue;
                                            unsigned int rm = cand_[idx] & mIB;
                                            while (rm) {
                                                const unsigned int one = rm & (~rm + 1U);
                                                bool changed = false;
                                                if (!removeCandidate(idx, firstDigit(one), changed)) return;
                                                if (changed) ++local;
                                                rm &= (rm - 1U);
                                            }
                                        }
                                        if (local > 0) {
                                            std::ostringstream ss;
                                            ss << "MSLS: line/box sector elimination removed " << local << " candidate(s)";
                                            pushDebugLog(ss.str());
                                            n = local;
                                            found = true;
                                        }
                                    });
                                }
                            });
                        }
                    });
                    if (found) return true;
                }
            }
        }
        return false;
    };

    if (tryMode(true)) return true;
    return tryMode(false);
}

bool SudokuAnalyzer::applyExocet(int& n) {
    n = 0;
    if (N_ < 6) return false;

    std::vector<int> unsolved;
    for (int i = 0; i < NN_; ++i) {
        if (grid_[i] == 0 && bits(cand_[i]) >= 2) unsolved.push_back(i);
    }
    if (unsolved.size() < 6U) return false;

    struct Perm {
        int b1d = 0;
        int b2d = 0;
        int expect1 = 0;
        int expect2 = 0;
    };

    for (std::size_t i = 0; i < unsolved.size(); ++i) {
        for (std::size_t j = i + 1; j < unsolved.size(); ++j) {
            const int b1 = unsolved[i];
            const int b2 = unsolved[j];
            if (box(row(b1), col(b1)) != box(row(b2), col(b2))) continue;
            const int baseBox = box(row(b1), col(b1));

            const unsigned int common = cand_[b1] & cand_[b2];
            const std::vector<int> cd = digitsFromMask(common);
            if (cd.size() < 2U) continue;

            for (std::size_t d1 = 0; d1 < cd.size(); ++d1) {
                for (std::size_t d2 = d1 + 1; d2 < cd.size(); ++d2) {
                    const int x = cd[d1];
                    const int y = cd[d2];
                    const unsigned int pairMask = bit(x) | bit(y);
                    auto buildBranch = [&](int baseCell, bool useRow) -> std::vector<int> {
                        std::vector<int> out;
                        if (useRow) {
                            const int r = row(baseCell);
                            for (int idx : houses_[r]) {
                                if (grid_[idx] != 0) continue;
                                if (box(row(idx), col(idx)) == baseBox) continue;
                                if ((cand_[idx] & pairMask) == 0U) continue;
                                out.push_back(idx);
                            }
                        } else {
                            const int c = col(baseCell);
                            for (int idx : houses_[N_ + c]) {
                                if (grid_[idx] != 0) continue;
                                if (box(row(idx), col(idx)) == baseBox) continue;
                                if ((cand_[idx] & pairMask) == 0U) continue;
                                out.push_back(idx);
                            }
                        }
                        std::sort(out.begin(), out.end());
                        out.erase(std::unique(out.begin(), out.end()), out.end());
                        return out;
                    };

                    const std::vector<std::pair<bool, bool>> mapModes = {
                        {false, false}, // col-col
                        {true, true},   // row-row
                        {true, false},  // row-col
                        {false, true}   // col-row
                    };

                    for (const auto& mode : mapModes) {
                        const std::vector<int> branch1 = buildBranch(b1, mode.first);
                        const std::vector<int> branch2 = buildBranch(b2, mode.second);
                        if (branch1.size() < 2U || branch2.size() < 2U) continue;

                        bool disjoint = true;
                        for (int t : branch1) {
                            if (std::find(branch2.begin(), branch2.end(), t) != branch2.end()) {
                                disjoint = false;
                                break;
                            }
                        }
                        if (!disjoint) continue;

                        std::set<int> targetSet(branch1.begin(), branch1.end());
                        targetSet.insert(branch2.begin(), branch2.end());
                        if (targetSet.size() > 14U) continue;

                        const std::string modeName =
                            std::string(mode.first ? "row" : "col") + "-" + (mode.second ? "row" : "col");
                        auto cellListText = [&](const std::vector<int>& cells) -> std::string {
                            std::ostringstream ss;
                            for (std::size_t k = 0; k < cells.size(); ++k) {
                                if (k > 0) ss << ",";
                                ss << cellName(cells[k]);
                            }
                            return ss.str();
                        };

                        std::map<std::string, bool> supportCache;
                        auto supported = [&](std::vector<std::pair<int, int>> asg) -> bool {
                            std::sort(asg.begin(), asg.end(), [](const auto& A, const auto& B) {
                                if (A.first != B.first) return A.first < B.first;
                                return A.second < B.second;
                            });
                            std::ostringstream key;
                            for (const auto& p : asg) key << p.first << "=" << p.second << ";";
                            const std::string k = key.str();
                            const auto it = supportCache.find(k);
                            if (it != supportCache.end()) return it->second;
                            const bool ok = hasLogicalSupportWithAssignments(asg);
                            supportCache[k] = ok;
                            return ok;
                        };

                        const std::vector<Perm> perms = {
                            {x, y, y, x},
                            {y, x, x, y}
                        };
                        std::vector<Perm> feasiblePerms;
                        for (const Perm& p : perms) {
                            const std::vector<std::pair<int, int>> bases = {
                                {b1, p.b1d}, {b2, p.b2d}
                            };
                            if (!supported(bases)) continue;

                            bool support1 = false;
                            for (int t : branch1) {
                                if ((cand_[t] & bit(p.expect1)) == 0U) continue;
                                std::vector<std::pair<int, int>> asg = bases;
                                asg.push_back({t, p.expect1});
                                if (supported(asg)) { support1 = true; break; }
                            }
                            if (!support1) continue;

                            bool support2 = false;
                            for (int t : branch2) {
                                if ((cand_[t] & bit(p.expect2)) == 0U) continue;
                                std::vector<std::pair<int, int>> asg = bases;
                                asg.push_back({t, p.expect2});
                                if (supported(asg)) { support2 = true; break; }
                            }
                            if (!support2) continue;

                            feasiblePerms.push_back(p);
                        }
                        if (feasiblePerms.empty()) continue;

                        bool patternLogged = false;
                        auto ensurePatternLogged = [&]() {
                            if (patternLogged) return;
                            std::ostringstream permSs;
                            for (std::size_t pi = 0; pi < feasiblePerms.size(); ++pi) {
                                if (pi > 0) permSs << " | ";
                                permSs << "("
                                       << cellName(b1) << "=" << feasiblePerms[pi].b1d << ", "
                                       << cellName(b2) << "=" << feasiblePerms[pi].b2d << ")";
                            }
                            std::ostringstream head;
                            head << "Exocet pattern: base{" << cellName(b1) << "," << cellName(b2)
                                 << "}, pair={" << x << "," << y << "}, mode=" << modeName
                                 << ", branch1={" << cellListText(branch1) << "}, branch2={"
                                 << cellListText(branch2) << "}, feasible=" << permSs.str();
                            pushDebugLog(head.str());
                            patternLogged = true;
                        };

                        int local = 0;
                        auto eliminateInBranch = [&](const std::vector<int>& branch, int which) -> bool {
                            for (int t : branch) {
                                if (grid_[t] != 0) continue;
                                const std::vector<int> dlist = digitsFromMask(cand_[t]);
                                for (int d : dlist) {
                                    bool supp = false;
                                    for (const Perm& p : feasiblePerms) {
                                        const int expected = (which == 1) ? p.expect1 : p.expect2;
                                        if ((d == x || d == y) && d != expected) continue;
                                        std::vector<std::pair<int, int>> asg = {
                                            {b1, p.b1d}, {b2, p.b2d}, {t, d}
                                        };
                                        if (supported(asg)) { supp = true; break; }
                                    }
                                    if (!supp) {
                                        // Konserwatywnie: usuwaj tylko kandydaty bez wsparcia implikacyjnego.
                                        if (supported({{t, d}})) continue;
                                        bool changed = false;
                                        if (!removeCandidate(t, d, changed)) return false;
                                        if (changed) {
                                            ensurePatternLogged();
                                            std::ostringstream ss;
                                            ss << "Exocet: remove " << d << " from " << cellName(t)
                                               << " (branch " << which << ", mode=" << modeName << ")";
                                            pushDebugLog(ss.str());
                                            ++local;
                                        }
                                    }
                                }
                            }
                            return true;
                        };

                        if (!eliminateInBranch(branch1, 1)) return false;
                        if (!eliminateInBranch(branch2, 2)) return false;

                        for (int t : targetSet) {
                            if (grid_[t] != 0) continue;
                            unsigned int extras = cand_[t] & ~pairMask;
                            while (extras) {
                                const unsigned int one = extras & (~extras + 1U);
                                const int z = firstDigit(one);
                                bool supp = false;
                                for (const Perm& p : feasiblePerms) {
                                    std::vector<std::pair<int, int>> asg = {
                                        {b1, p.b1d}, {b2, p.b2d}, {t, z}
                                    };
                                    if (supported(asg)) { supp = true; break; }
                                }
                                if (!supp) {
                                    if (supported({{t, z}})) {
                                        extras &= (extras - 1U);
                                        continue;
                                    }
                                    bool changed = false;
                                    if (!removeCandidate(t, z, changed)) return false;
                                    if (changed) {
                                        ensurePatternLogged();
                                        std::ostringstream ss;
                                        ss << "Exocet: remove " << z << " from " << cellName(t)
                                           << " (target extra, mode=" << modeName << ")";
                                        pushDebugLog(ss.str());
                                        ++local;
                                    }
                                }
                                extras &= (extras - 1U);
                            }
                        }

                        if (local > 0) {
                            n = local;
                            return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applySeniorExocet(int& n) {
    n = 0;
    if (N_ < 6) return false;

    std::vector<int> unsolved;
    for (int i = 0; i < NN_; ++i) {
        if (grid_[i] == 0 && bits(cand_[i]) >= 2) unsolved.push_back(i);
    }
    if (unsolved.size() < 6U) return false;

    struct Perm {
        int b1d = 0;
        int b2d = 0;
        int expect1 = 0;
        int expect2 = 0;
    };

    for (std::size_t i = 0; i < unsolved.size(); ++i) {
        for (std::size_t j = i + 1; j < unsolved.size(); ++j) {
            const int b1 = unsolved[i];
            const int b2 = unsolved[j];
            if (box(row(b1), col(b1)) != box(row(b2), col(b2))) continue;
            const int baseBox = box(row(b1), col(b1));

            const unsigned int common = cand_[b1] & cand_[b2];
            const std::vector<int> cd = digitsFromMask(common);
            if (cd.size() < 2U) continue;

            for (std::size_t d1 = 0; d1 < cd.size(); ++d1) {
                for (std::size_t d2 = d1 + 1; d2 < cd.size(); ++d2) {
                    const int x = cd[d1];
                    const int y = cd[d2];
                    const unsigned int pairMask = bit(x) | bit(y);

                    auto buildBranch = [&](int baseCell, bool useRow) -> std::vector<int> {
                        std::vector<int> out;
                        if (useRow) {
                            const int r = row(baseCell);
                            for (int idx : houses_[r]) {
                                if (grid_[idx] != 0) continue;
                                if (box(row(idx), col(idx)) == baseBox) continue;
                                if ((cand_[idx] & pairMask) == 0U) continue;
                                out.push_back(idx);
                            }
                        } else {
                            const int c = col(baseCell);
                            for (int idx : houses_[N_ + c]) {
                                if (grid_[idx] != 0) continue;
                                if (box(row(idx), col(idx)) == baseBox) continue;
                                if ((cand_[idx] & pairMask) == 0U) continue;
                                out.push_back(idx);
                            }
                        }
                        std::sort(out.begin(), out.end());
                        out.erase(std::unique(out.begin(), out.end()), out.end());
                        return out;
                    };

                    const std::vector<std::pair<bool, bool>> mapModes = {
                        {false, false}, {true, true}, {true, false}, {false, true}
                    };
                    for (const auto& mode : mapModes) {
                        const std::vector<int> branch1 = buildBranch(b1, mode.first);
                        const std::vector<int> branch2 = buildBranch(b2, mode.second);
                        if (branch1.size() < 2U || branch2.size() < 2U) continue;

                        bool disjoint = true;
                        for (int t : branch1) {
                            if (std::find(branch2.begin(), branch2.end(), t) != branch2.end()) {
                                disjoint = false;
                                break;
                            }
                        }
                        if (!disjoint) continue;

                        std::map<std::string, bool> supportCache;
                        auto supported = [&](std::vector<std::pair<int, int>> asg) -> bool {
                            std::sort(asg.begin(), asg.end(), [](const auto& A, const auto& B) {
                                if (A.first != B.first) return A.first < B.first;
                                return A.second < B.second;
                            });
                            std::ostringstream key;
                            for (const auto& p : asg) key << p.first << "=" << p.second << ";";
                            const std::string k = key.str();
                            const auto it = supportCache.find(k);
                            if (it != supportCache.end()) return it->second;
                            const bool ok = hasLogicalSupportWithAssignments(asg);
                            supportCache[k] = ok;
                            return ok;
                        };

                        const std::vector<Perm> perms = {
                            {x, y, y, x},
                            {y, x, x, y}
                        };
                        std::vector<Perm> feasiblePerms;
                        for (const Perm& p : perms) {
                            const std::vector<std::pair<int, int>> bases = {
                                {b1, p.b1d}, {b2, p.b2d}
                            };
                            if (!supported(bases)) continue;

                            bool support1 = false;
                            for (int t : branch1) {
                                if ((cand_[t] & bit(p.expect1)) == 0U) continue;
                                std::vector<std::pair<int, int>> asg = bases;
                                asg.push_back({t, p.expect1});
                                if (supported(asg)) { support1 = true; break; }
                            }
                            if (!support1) continue;

                            bool support2 = false;
                            for (int t : branch2) {
                                if ((cand_[t] & bit(p.expect2)) == 0U) continue;
                                std::vector<std::pair<int, int>> asg = bases;
                                asg.push_back({t, p.expect2});
                                if (supported(asg)) { support2 = true; break; }
                            }
                            if (!support2) continue;

                            feasiblePerms.push_back(p);
                        }
                        if (feasiblePerms.size() != 1U) continue;
                        const Perm p = feasiblePerms[0];

                        int local = 0;
                        const int bad1 = (p.b1d == x) ? y : x;
                        const int bad2 = (p.b2d == x) ? y : x;
                        bool ch1 = false, ch2 = false;
                        if (!removeCandidate(b1, bad1, ch1)) return false;
                        if (ch1) ++local;
                        if (!removeCandidate(b2, bad2, ch2)) return false;
                        if (ch2) ++local;

                        if (local > 0) {
                            std::ostringstream ss;
                            ss << "SeniorExocet: fixed base pair "
                               << cellName(b1) << "=" << p.b1d << ", "
                               << cellName(b2) << "=" << p.b2d;
                            pushDebugLog(ss.str());
                            n = local;
                            return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

bool SudokuAnalyzer::applyPatternOverlayMethod(int& n) {
    n = 0;
    if (N_ != 9) return false;

    for (int d = 1; d <= N_; ++d) {
        std::vector<std::vector<int>> rowOpts(N_);
        bool impossible = false;
        for (int r = 0; r < N_; ++r) {
            int fixedCol = -1;
            for (int c = 0; c < N_; ++c) {
                if (grid_[r * N_ + c] == d) {
                    fixedCol = c;
                    break;
                }
            }
            if (fixedCol >= 0) {
                rowOpts[r].push_back(fixedCol);
            } else {
                for (int c = 0; c < N_; ++c) {
                    const int idx = r * N_ + c;
                    if (grid_[idx] == 0 && (cand_[idx] & bit(d))) rowOpts[r].push_back(c);
                }
            }
            if (rowOpts[r].empty()) {
                impossible = true;
                break;
            }
        }
        if (impossible) continue;

        std::vector<int> rowOrder(N_);
        for (int r = 0; r < N_; ++r) rowOrder[r] = r;
        std::sort(rowOrder.begin(), rowOrder.end(), [&](int a, int b) {
            return rowOpts[a].size() < rowOpts[b].size();
        });

        std::vector<int> chosenCol(N_, -1);
        std::vector<char> canBeTrue(NN_, 0);
        int patterns = 0;
        const int maxPatterns = 20000;
        bool overflow = false;

        std::function<void(int, int, int)> dfs = [&](int depth, int usedCols, int usedBoxes) {
            if (overflow) return;
            if (depth == N_) {
                ++patterns;
                for (int r = 0; r < N_; ++r) {
                    if (chosenCol[r] >= 0) canBeTrue[r * N_ + chosenCol[r]] = 1;
                }
                if (patterns >= maxPatterns) overflow = true;
                return;
            }

            const int r = rowOrder[depth];
            for (int c : rowOpts[r]) {
                const int b = box(r, c);
                if (usedCols & (1 << c)) continue;
                if (usedBoxes & (1 << b)) continue;
                chosenCol[r] = c;

                bool futureOk = true;
                for (int nd = depth + 1; nd < N_; ++nd) {
                    const int rr = rowOrder[nd];
                    bool any = false;
                    for (int cc : rowOpts[rr]) {
                        const int bb = box(rr, cc);
                        if ((usedCols & (1 << cc)) != 0) continue;
                        if ((usedBoxes & (1 << bb)) != 0) continue;
                        if (cc == c || bb == b) continue;
                        any = true;
                        break;
                    }
                    if (!any) {
                        futureOk = false;
                        break;
                    }
                }
                if (futureOk) dfs(depth + 1, usedCols | (1 << c), usedBoxes | (1 << b));
                chosenCol[r] = -1;
                if (overflow) return;
            }
        };
        dfs(0, 0, 0);
        if (overflow || patterns <= 0) continue;

        int local = 0;
        for (int idx = 0; idx < NN_; ++idx) {
            if (grid_[idx] != 0 || (cand_[idx] & bit(d)) == 0U) continue;
            if (canBeTrue[idx]) continue;
            bool changed = false;
            if (!removeCandidate(idx, d, changed)) return false;
            if (changed) ++local;
        }
        if (local > 0) {
            std::ostringstream ss;
            ss << "POM: digit " << d << " removed from " << local << " cell(s), patterns=" << patterns;
            pushDebugLog(ss.str());
            n = local;
            return true;
        }
    }
    return false;
}

bool SudokuAnalyzer::applyForcingChains(int& n) {
    n = 0;
    if (N_ != 9) return false;

    std::map<std::string, bool> supportCache;
    auto supported = [&](std::vector<std::pair<int, int>> asg) -> bool {
        std::sort(asg.begin(), asg.end(), [](const auto& A, const auto& B) {
            if (A.first != B.first) return A.first < B.first;
            return A.second < B.second;
        });
        std::ostringstream key;
        for (const auto& p : asg) key << p.first << "=" << p.second << ";";
        const std::string k = key.str();
        const auto it = supportCache.find(k);
        if (it != supportCache.end()) return it->second;
        const bool ok = hasLogicalSupportWithAssignments(asg);
        supportCache[k] = ok;
        return ok;
    };

    auto evalBranches = [&](const std::vector<std::pair<int, int>>& branches) -> bool {
        std::vector<std::pair<int, int>> feasible;
        feasible.reserve(branches.size());
        for (const auto& b : branches) {
            if (supported({b})) feasible.push_back(b);
        }
        if (feasible.size() < 2U) return false;

        std::vector<int> probeCells;
        probeCells.reserve(NN_);
        for (int i = 0; i < NN_; ++i) {
            if (grid_[i] == 0 && cand_[i] != 0U) {
                probeCells.push_back(i);
            }
        }

        for (int c : probeCells) {
            const std::vector<int> dlist = digitsFromMask(cand_[c]);
            for (int d : dlist) {
                bool supportedAnywhere = false;
                for (const auto& br : feasible) {
                    bool supp = false;
                    if (br.first == c) {
                        supp = (br.second == d);
                    } else {
                        supp = supported({br, {c, d}});
                    }
                    if (supp) {
                        supportedAnywhere = true;
                        break;
                    }
                }
                if (!supportedAnywhere) {
                    bool changed = false;
                    if (!removeCandidate(c, d, changed)) return false;
                    if (changed) {
                        std::ostringstream ss;
                        ss << "ForcingChains: remove " << d << " from " << cellName(c)
                           << " (all branches contradicted)";
                        pushDebugLog(ss.str());
                        n = 1;
                        return true;
                    }
                }
            }
        }
        return false;
    };

    for (int cell = 0; cell < NN_; ++cell) {
        if (grid_[cell] != 0) continue;
        const int bc = bits(cand_[cell]);
        if (bc < 2 || bc > 3) continue;
        std::vector<std::pair<int, int>> branches;
        for (int d : digitsFromMask(cand_[cell])) {
            branches.push_back({cell, d});
        }
        if (evalBranches(branches)) return true;
    }

    for (int h = 0; h < 3 * N_; ++h) {
        for (int d = 1; d <= N_; ++d) {
            std::vector<int> where;
            for (int idx : houses_[h]) {
                if (grid_[idx] == 0 && (cand_[idx] & bit(d))) where.push_back(idx);
            }
            if (where.size() < 2U || where.size() > 3U) continue;
            std::vector<std::pair<int, int>> branches;
            for (int idx : where) {
                branches.push_back({idx, d});
            }
            if (evalBranches(branches)) return true;
        }
    }

    return false;
}

bool SudokuAnalyzer::applyGroupedXCycle(int& n) {
    n = 0;

    auto eliminateSeeingBoth = [&](int a, int b, int d) -> bool {
        int local = 0;
        for (int i = 0; i < NN_; ++i) {
            if (i == a || i == b || grid_[i] != 0) continue;
            if ((cand_[i] & bit(d)) == 0U) continue;
            if (!isPeerCell(i, a) || !isPeerCell(i, b)) continue;
            bool changed = false;
            if (!removeCandidate(i, d, changed)) return false;
            if (changed) ++local;
        }
        if (local > 0) {
            n = local;
            return true;
        }
        return false;
    };

    for (int d = 1; d <= N_; ++d) {
        const int bpr = N_ / BC_;
        for (int b = 0; b < N_; ++b) {
            const int sr = (b / bpr) * BR_;
            const int sc = (b % bpr) * BC_;

            for (int rr = 0; rr < BR_; ++rr) {
                const int r = sr + rr;
                std::vector<int> group;
                for (int cc = 0; cc < BC_; ++cc) {
                    const int c = sc + cc;
                    const int idx = r * N_ + c;
                    if (grid_[idx] == 0 && (cand_[idx] & bit(d))) group.push_back(idx);
                }
                if (group.size() < 2U) continue;

                std::vector<int> rowAll;
                for (int idx : houses_[r]) if (grid_[idx] == 0 && (cand_[idx] & bit(d))) rowAll.push_back(idx);
                if (rowAll.size() != group.size() + 1U) continue;
                int outsideRow = -1;
                for (int idx : rowAll) {
                    if (std::find(group.begin(), group.end(), idx) == group.end()) {
                        outsideRow = idx;
                        break;
                    }
                }
                if (outsideRow < 0) continue;

                std::vector<int> boxAll;
                for (int idx : houses_[2 * N_ + b]) if (grid_[idx] == 0 && (cand_[idx] & bit(d))) boxAll.push_back(idx);
                if (boxAll.size() != group.size() + 1U) continue;
                int outsideBox = -1;
                for (int idx : boxAll) {
                    if (std::find(group.begin(), group.end(), idx) == group.end()) {
                        outsideBox = idx;
                        break;
                    }
                }
                if (outsideBox < 0 || outsideBox == outsideRow) continue;

                if (eliminateSeeingBoth(outsideRow, outsideBox, d)) return true;
            }

            for (int cc = 0; cc < BC_; ++cc) {
                const int c = sc + cc;
                std::vector<int> group;
                for (int rr = 0; rr < BR_; ++rr) {
                    const int r = sr + rr;
                    const int idx = r * N_ + c;
                    if (grid_[idx] == 0 && (cand_[idx] & bit(d))) group.push_back(idx);
                }
                if (group.size() < 2U) continue;

                std::vector<int> colAll;
                for (int idx : houses_[N_ + c]) if (grid_[idx] == 0 && (cand_[idx] & bit(d))) colAll.push_back(idx);
                if (colAll.size() != group.size() + 1U) continue;
                int outsideCol = -1;
                for (int idx : colAll) {
                    if (std::find(group.begin(), group.end(), idx) == group.end()) {
                        outsideCol = idx;
                        break;
                    }
                }
                if (outsideCol < 0) continue;

                std::vector<int> boxAll;
                for (int idx : houses_[2 * N_ + b]) if (grid_[idx] == 0 && (cand_[idx] & bit(d))) boxAll.push_back(idx);
                if (boxAll.size() != group.size() + 1U) continue;
                int outsideBox = -1;
                for (int idx : boxAll) {
                    if (std::find(group.begin(), group.end(), idx) == group.end()) {
                        outsideBox = idx;
                        break;
                    }
                }
                if (outsideBox < 0 || outsideBox == outsideCol) continue;

                if (eliminateSeeingBoth(outsideCol, outsideBox, d)) return true;
            }
        }
    }

    return false;
}

bool SudokuAnalyzer::applyGroupedAIC(int& n) {
    n = 0;
    if (N_ != 9) return false;

    auto branchText = [&](const std::vector<std::pair<int, int>>& branches) -> std::string {
        std::ostringstream ss;
        for (std::size_t i = 0; i < branches.size(); ++i) {
            if (i > 0) ss << " | ";
            ss << cellName(branches[i].first) << "=" << branches[i].second;
        }
        return ss.str();
    };

    std::map<std::string, bool> supportCache;
    auto supported = [&](std::vector<std::pair<int, int>> asg) -> bool {
        std::sort(asg.begin(), asg.end(), [](const auto& A, const auto& B) {
            if (A.first != B.first) return A.first < B.first;
            return A.second < B.second;
        });
        std::ostringstream key;
        for (const auto& p : asg) key << p.first << "=" << p.second << ";";
        const std::string k = key.str();
        const auto it = supportCache.find(k);
        if (it != supportCache.end()) return it->second;
        const bool ok = hasLogicalSupportWithAssignments(asg);
        supportCache[k] = ok;
        return ok;
    };

    auto evalBranches = [&](const std::vector<std::pair<int, int>>& branches) -> bool {
        std::vector<std::pair<int, int>> feasible;
        for (const auto& b : branches) {
            if (supported({b})) feasible.push_back(b);
        }
        if (feasible.size() < 2U) return false;

        for (int c = 0; c < NN_; ++c) {
            if (grid_[c] != 0) continue;
            for (int d : digitsFromMask(cand_[c])) {
                bool supp = false;
                for (const auto& br : feasible) {
                    if (br.first == c) {
                        if (br.second == d) { supp = true; break; }
                    } else if (supported({br, {c, d}})) {
                        supp = true;
                        break;
                    }
                }
                if (!supp) {
                    bool changed = false;
                    if (!removeCandidate(c, d, changed)) return false;
                    if (changed) {
                        std::ostringstream ss;
                        ss << "GroupedAIC: branches={" << branchText(feasible)
                           << "} -> remove " << d << " from " << cellName(c);
                        pushDebugLog(ss.str());
                        n = 1;
                        return true;
                    }
                }
            }
        }
        return false;
    };

    const int bpr = N_ / BC_;
    for (int d = 1; d <= N_; ++d) {
        for (int b = 0; b < N_; ++b) {
            const int sr = (b / bpr) * BR_;
            const int sc = (b % bpr) * BC_;

            for (int rr = 0; rr < BR_; ++rr) {
                const int r = sr + rr;
                std::vector<int> group;
                for (int cc = 0; cc < BC_; ++cc) {
                    const int idx = r * N_ + (sc + cc);
                    if (grid_[idx] == 0 && (cand_[idx] & bit(d))) group.push_back(idx);
                }
                if (group.size() >= 2U && group.size() <= 4U) {
                    std::vector<std::pair<int, int>> branches;
                    for (int idx : group) branches.push_back({idx, d});
                    if (evalBranches(branches)) return true;
                }
            }

            for (int cc = 0; cc < BC_; ++cc) {
                const int c = sc + cc;
                std::vector<int> group;
                for (int rr = 0; rr < BR_; ++rr) {
                    const int idx = (sr + rr) * N_ + c;
                    if (grid_[idx] == 0 && (cand_[idx] & bit(d))) group.push_back(idx);
                }
                if (group.size() >= 2U && group.size() <= 4U) {
                    std::vector<std::pair<int, int>> branches;
                    for (int idx : group) branches.push_back({idx, d});
                    if (evalBranches(branches)) return true;
                }
            }
        }
    }

    return false;
}

bool SudokuAnalyzer::applyThreeDMedusa(int& n) {
    n = 0;

    struct Node { int cell = -1; int digit = 0; };
    std::vector<Node> nodes;
    std::vector<int> nodeId(NN_ * (N_ + 1), -1);
    auto key = [&](int cell, int digit) { return cell * (N_ + 1) + digit; };

    for (int c = 0; c < NN_; ++c) {
        if (grid_[c] != 0) continue;
        for (int d : digitsFromMask(cand_[c])) {
            const int id = static_cast<int>(nodes.size());
            nodes.push_back({c, d});
            nodeId[key(c, d)] = id;
        }
    }
    if (nodes.empty()) return false;

    std::vector<std::vector<int>> g(nodes.size());
    auto addEdge = [&](int a, int b) {
        if (a < 0 || b < 0 || a == b) return;
        g[a].push_back(b);
        g[b].push_back(a);
    };

    for (int c = 0; c < NN_; ++c) {
        if (grid_[c] != 0 || bits(cand_[c]) != 2) continue;
        const std::vector<int> ds = digitsFromMask(cand_[c]);
        addEdge(nodeId[key(c, ds[0])], nodeId[key(c, ds[1])]);
    }
    for (int h = 0; h < 3 * N_; ++h) {
        for (int d = 1; d <= N_; ++d) {
            std::vector<int> where;
            for (int idx : houses_[h]) {
                if (grid_[idx] == 0 && (cand_[idx] & bit(d))) where.push_back(idx);
            }
            if (where.size() == 2U) {
                addEdge(nodeId[key(where[0], d)], nodeId[key(where[1], d)]);
            }
        }
    }

    std::vector<int> comp(nodes.size(), -1), color(nodes.size(), -1);
    int compCnt = 0;
    for (int s = 0; s < static_cast<int>(nodes.size()); ++s) {
        if (comp[s] != -1) continue;
        std::vector<int> q = {s};
        comp[s] = compCnt;
        color[s] = 0;
        for (std::size_t qi = 0; qi < q.size(); ++qi) {
            const int u = q[qi];
            for (int v : g[u]) {
                if (comp[v] == -1) {
                    comp[v] = compCnt;
                    color[v] = 1 - color[u];
                    q.push_back(v);
                }
            }
        }
        ++compCnt;
    }

    std::vector<std::array<bool, 2>> bad(compCnt, {false, false});

    for (int c = 0; c < NN_; ++c) {
        std::map<int, std::array<int, 2>> cnt;
        for (int d : digitsFromMask(cand_[c])) {
            const int id = nodeId[key(c, d)];
            if (id < 0) continue;
            cnt[comp[id]][color[id]]++;
        }
        for (const auto& kv : cnt) {
            for (int col = 0; col < 2; ++col) {
                if (kv.second[col] >= 2) bad[kv.first][col] = true;
            }
        }
    }

    for (int h = 0; h < 3 * N_; ++h) {
        for (int d = 1; d <= N_; ++d) {
            std::map<int, std::array<int, 2>> cnt;
            for (int c : houses_[h]) {
                const int id = nodeId[key(c, d)];
                if (id < 0) continue;
                cnt[comp[id]][color[id]]++;
            }
            for (const auto& kv : cnt) {
                for (int col = 0; col < 2; ++col) {
                    if (kv.second[col] >= 2) bad[kv.first][col] = true;
                }
            }
        }
    }

    int local = 0;
    for (int id = 0; id < static_cast<int>(nodes.size()); ++id) {
        if (!bad[comp[id]][color[id]]) continue;
        bool changed = false;
        if (!removeCandidate(nodes[id].cell, nodes[id].digit, changed)) return false;
        if (changed) ++local;
    }
    if (local > 0) { n = local; return true; }

    for (int c = 0; c < NN_; ++c) {
        if (grid_[c] != 0) continue;
        for (int d : digitsFromMask(cand_[c])) {
            std::map<int, std::array<bool, 2>> seen;
            for (int id = 0; id < static_cast<int>(nodes.size()); ++id) {
                if (nodes[id].digit != d) continue;
                if (!isPeerCell(c, nodes[id].cell)) continue;
                seen[comp[id]][color[id]] = true;
            }
            bool remove = false;
            for (const auto& kv : seen) {
                if (kv.second[0] && kv.second[1]) { remove = true; break; }
            }
            if (remove) {
                bool changed = false;
                if (!removeCandidate(c, d, changed)) return false;
                if (changed) { n = 1; return true; }
            }
        }
    }

    return false;
}

void SudokuAnalyzer::logicalSolve() {
    bool progress = true;
    while (!contradiction_ && !solved() && progress) {
        progress = false;
        int n = 0;
        if (applyNakedSingles(n))      { use(Strategy::NakedSingle, n); progress = true; continue; }
        if (applyHiddenSingles(n))     { use(Strategy::HiddenSingle, n); progress = true; continue; }
        if (applyNakedSubset(2, n))    { use(Strategy::NakedPair, n); progress = true; continue; }
        if (applyHiddenSubset(2, n))   { use(Strategy::HiddenPair, n); progress = true; continue; }
        if (applyPointingPairsTriples(n)) { use(Strategy::PointingPairsTriples, n); progress = true; continue; }
        if (applyBoxLineReduction(n))  { use(Strategy::BoxLineReduction, n); progress = true; continue; }
        if (applyNakedSubset(3, n))    { use(Strategy::NakedTriple, n); progress = true; continue; }
        if (applyHiddenSubset(3, n))   { use(Strategy::HiddenTriple, n); progress = true; continue; }
        if (applyNakedSubset(4, n))    { use(Strategy::NakedQuad, n); progress = true; continue; }
        if (applyHiddenSubset(4, n))   { use(Strategy::HiddenQuad, n); progress = true; continue; }
        if (applyFish(2, n))           { use(Strategy::XWing, n); progress = true; continue; }
        if (applyYWing(n))             { use(Strategy::YWing, n); progress = true; continue; }
        if (applyXYZWing(n))           { use(Strategy::XYZWing, n); progress = true; continue; }
        if (applyWXYZWing(n))          { use(Strategy::WXYZWing, n); progress = true; continue; }
        if (applyFish(3, n))           { use(Strategy::Swordfish, n); progress = true; continue; }
        if (applyFinnedFish(3, n))     { use(Strategy::FinnedSwordfish, n); progress = true; continue; }
        if (applyFrankenMutantFish(2, n)) { use(Strategy::FrankenMutantFish, n); progress = true; continue; }
        if (applyKrakenFish(n))        { use(Strategy::KrakenFish, n); progress = true; continue; }
        if (applySkyscraper(n))        { use(Strategy::Skyscraper, n); progress = true; continue; }
        if (applyTwoStringKite(n))     { use(Strategy::TwoStringKite, n); progress = true; continue; }
        if (applySimpleColoring(n))    { use(Strategy::SimpleColoring, n); progress = true; continue; }
        if (applyThreeDMedusa(n))      { use(Strategy::ThreeDMedusa, n); progress = true; continue; }
        if (applyFish(4, n))           { use(Strategy::Jellyfish, n); progress = true; continue; }
        if (applyFinnedXWingSashimi(n)) { use(Strategy::FinnedXWingSashimi, n); progress = true; continue; }
        if (applyFinnedFish(4, n))     { use(Strategy::FinnedJellyfish, n); progress = true; continue; }
        if (applyFrankenMutantFish(3, n)) { use(Strategy::FrankenMutantFish, n); progress = true; continue; }
        if (applyEmptyRectangle(n))    { use(Strategy::EmptyRectangle, n); progress = true; continue; }
        if (applyUniqueRectangleType1(n)) { use(Strategy::UniqueRectangle, n); progress = true; continue; }
        if (applyUniqueRectangleType2to6(n)) { use(Strategy::UniqueRectangle, n); progress = true; continue; }
        if (applyUniqueLoop(n))        { use(Strategy::UniqueLoop, n); progress = true; continue; }
        if (applyBivalueOddagon(n))    { use(Strategy::BivalueOddagon, n); progress = true; continue; }
        if (applyAvoidableRectangle(n)) { use(Strategy::AvoidableRectangle, n); progress = true; continue; }
        if (applyBUGPlus1(n))          { use(Strategy::BUGPlus1, n); progress = true; continue; }
        if (applyRemotePairs(n))       { use(Strategy::RemotePairs, n); progress = true; continue; }
        if (applyWWing(n))             { use(Strategy::WWing, n); progress = true; continue; }
        if (applyGroupedXCycle(n))     { use(Strategy::GroupedXCycle, n); progress = true; continue; }
        if (applyXChain(n))            { use(Strategy::XChain, n); progress = true; continue; }
        if (applyXYChain(n))           { use(Strategy::XYChain, n); progress = true; continue; }
        if (applyGroupedAIC(n))        { use(Strategy::GroupedAIC, n); progress = true; continue; }
        if (applyAIC(n))               { use(Strategy::AIC, n); progress = true; continue; }
        if (applyContinuousNiceLoop(n)) { use(Strategy::ContinuousNiceLoop, n); progress = true; continue; }
        if (applySKLoop(n))            { use(Strategy::SKLoop, n); progress = true; continue; }
        if (applyALSXZ(n))             { use(Strategy::ALSXZ, n); progress = true; continue; }
        if (applyALSXYWing(n))         { use(Strategy::ALSXYWing, n); progress = true; continue; }
        if (applyALSChain(n))          { use(Strategy::ALSChain, n); progress = true; continue; }
        if (applyDeathBlossom(n))      { use(Strategy::DeathBlossom, n); progress = true; continue; }
        if (applySueDeCoq(n))          { use(Strategy::SueDeCoq, n); progress = true; continue; }
        if (applyMSLS(n))              { use(Strategy::MSLS, n); progress = true; continue; }
        if (applyExocet(n))            { use(Strategy::Exocet, n); progress = true; continue; }
        if (applySeniorExocet(n))      { use(Strategy::SeniorExocet, n); progress = true; continue; }
        if (applyPatternOverlayMethod(n)) { use(Strategy::PatternOverlayMethod, n); progress = true; continue; }
        if (applyForcingChains(n))     { use(Strategy::ForcingChains, n); progress = true; continue; }
    }
}

AnalysisReport SudokuAnalyzer::run() {
    AnalysisReport r;
    debug_logic_logs_.clear();
    debug_logic_truncated_ = false;
    r.initial_clues = cluesCount();
    if (!contradiction_) logicalSolve();
    r.contradiction = contradiction_;
    r.solved_logically = (!contradiction_ && solved());

    if (!r.contradiction && !r.solved_logically) {
        const BacktrackingSolveStats bt = solveWithBacktracking(b_, grid_);
        r.solved_with_backtracking = bt.solved;
        r.backtracking_nodes = bt.nodes;
        r.backtracking_decisions = bt.decisions;
        r.backtracking_backtracks = bt.backtracks;
        if (bt.solved) {
            const int btUsage = static_cast<int>(std::max(1LL, bt.decisions));
            use(Strategy::Backtracking, btUsage);
        }
    }
    r.requires_guessing = (!contradiction_ && !r.solved_logically);
    r.strategy_usage = usage_;
    r.debug_logic_logs = debug_logic_logs_;

    if (r.contradiction) {
        r.hardest_strategy = "Sprzeczna plansza (bledne dane wejsciowe)";
        r.hardest_rank = 0;
    } else if (r.solved_with_backtracking) {
        r.hardest_strategy = "Backtracking";
        r.hardest_rank = strategyRank(Strategy::Backtracking);
    } else if (r.requires_guessing) {
        r.hardest_strategy = "Wymaga zgadywania/Backtrackingu";
        r.hardest_rank = 100;
    } else if (hardest_rank_ == 0) {
        r.hardest_strategy = "Brak (same wskazowki)";
        r.hardest_rank = 0;
    } else {
        r.hardest_strategy = hardest_name_;
        r.hardest_rank = hardest_rank_;
    }
    return r;
}

BacktrackingCounter::BacktrackingCounter(int br, int bc, int n, std::vector<int> grid)
    : BR_(br), BC_(bc), N_(n), NN_(n * n), all_((1U << n) - 1U), grid_(std::move(grid)) {}

unsigned int BacktrackingCounter::allowed(int idx) const {
    if (grid_[idx] != 0) return bit(grid_[idx]);
    const int r = row(idx), c = col(idx), b = box(r, c);
    unsigned int m = all_;
    for (int i = 0; i < N_; ++i) {
        const int ri = r * N_ + i, ci = i * N_ + c;
        if (grid_[ri] != 0) m &= ~bit(grid_[ri]);
        if (grid_[ci] != 0) m &= ~bit(grid_[ci]);
    }
    const int bpr = N_ / BC_;
    const int sr = (b / bpr) * BR_, sc = (b % bpr) * BC_;
    for (int dr = 0; dr < BR_; ++dr)
        for (int dc = 0; dc < BC_; ++dc) {
            const int i = (sr + dr) * N_ + (sc + dc);
            if (grid_[i] != 0) m &= ~bit(grid_[i]);
        }
    return m;
}

bool BacktrackingCounter::validState() const {
    for (int r = 0; r < N_; ++r) {
        unsigned int seen = 0U;
        for (int c = 0; c < N_; ++c) {
            const int v = grid_[r * N_ + c];
            if (v == 0) continue;
            const unsigned int b = bit(v);
            if (seen & b) return false;
            seen |= b;
        }
    }
    for (int c = 0; c < N_; ++c) {
        unsigned int seen = 0U;
        for (int r = 0; r < N_; ++r) {
            const int v = grid_[r * N_ + c];
            if (v == 0) continue;
            const unsigned int b = bit(v);
            if (seen & b) return false;
            seen |= b;
        }
    }
    const int bpr = N_ / BC_;
    for (int b = 0; b < N_; ++b) {
        unsigned int seen = 0U;
        const int sr = (b / bpr) * BR_, sc = (b % bpr) * BC_;
        for (int dr = 0; dr < BR_; ++dr) {
            for (int dc = 0; dc < BC_; ++dc) {
                const int v = grid_[(sr + dr) * N_ + (sc + dc)];
                if (v == 0) continue;
                const unsigned int bt = bit(v);
                if (seen & bt) return false;
                seen |= bt;
            }
        }
    }
    return true;
}

void BacktrackingCounter::search() {
    if (solutions_ >= limit_) return;
    int best = -1, bestc = INT_MAX;
    unsigned int bm = 0U;
    for (int i = 0; i < NN_; ++i) {
        if (grid_[i] != 0) continue;
        const unsigned int m = allowed(i);
        const int c = bits(m);
        if (c == 0) return;
        if (c < bestc) {
            best = i; bestc = c; bm = m;
            if (c == 1) break;
        }
    }
    if (best == -1) { ++solutions_; return; }
    while (bm) {
        const unsigned int one = bm & (~bm + 1U);
        grid_[best] = firstDigit(one);
        search();
        grid_[best] = 0;
        if (solutions_ >= limit_) return;
        bm &= (bm - 1U);
    }
}

int BacktrackingCounter::countSolutions(int limit) {
    limit_ = limit;
    solutions_ = 0;
    if (!validState()) return 0;
    search();
    return solutions_;
}

BacktrackingSolver::BacktrackingSolver(int br, int bc, int n, std::vector<int> grid)
    : BR_(br), BC_(bc), N_(n), NN_(n * n), all_((1U << n) - 1U), grid_(std::move(grid)) {}

unsigned int BacktrackingSolver::allowed(int idx) const {
    if (grid_[idx] != 0) return bit(grid_[idx]);
    const int r = row(idx), c = col(idx), b = box(r, c);
    unsigned int m = all_;
    for (int i = 0; i < N_; ++i) {
        const int ri = r * N_ + i, ci = i * N_ + c;
        if (grid_[ri] != 0) m &= ~bit(grid_[ri]);
        if (grid_[ci] != 0) m &= ~bit(grid_[ci]);
    }
    const int bpr = N_ / BC_;
    const int sr = (b / bpr) * BR_, sc = (b % bpr) * BC_;
    for (int dr = 0; dr < BR_; ++dr) {
        for (int dc = 0; dc < BC_; ++dc) {
            const int i = (sr + dr) * N_ + (sc + dc);
            if (grid_[i] != 0) m &= ~bit(grid_[i]);
        }
    }
    return m;
}

bool BacktrackingSolver::validState() const {
    for (int r = 0; r < N_; ++r) {
        unsigned int seen = 0U;
        for (int c = 0; c < N_; ++c) {
            const int v = grid_[r * N_ + c];
            if (v == 0) continue;
            const unsigned int b = bit(v);
            if (seen & b) return false;
            seen |= b;
        }
    }
    for (int c = 0; c < N_; ++c) {
        unsigned int seen = 0U;
        for (int r = 0; r < N_; ++r) {
            const int v = grid_[r * N_ + c];
            if (v == 0) continue;
            const unsigned int b = bit(v);
            if (seen & b) return false;
            seen |= b;
        }
    }
    const int bpr = N_ / BC_;
    for (int b = 0; b < N_; ++b) {
        unsigned int seen = 0U;
        const int sr = (b / bpr) * BR_, sc = (b % bpr) * BC_;
        for (int dr = 0; dr < BR_; ++dr) {
            for (int dc = 0; dc < BC_; ++dc) {
                const int v = grid_[(sr + dr) * N_ + (sc + dc)];
                if (v == 0) continue;
                const unsigned int bt = bit(v);
                if (seen & bt) return false;
                seen |= bt;
            }
        }
    }
    return true;
}

bool BacktrackingSolver::search() {
    ++stats_.nodes;

    int best = -1;
    int bestc = INT_MAX;
    unsigned int bm = 0U;
    for (int i = 0; i < NN_; ++i) {
        if (grid_[i] != 0) continue;
        const unsigned int m = allowed(i);
        const int c = bits(m);
        if (c == 0) return false;
        if (c < bestc) {
            best = i;
            bestc = c;
            bm = m;
            if (c == 1) break;
        }
    }

    if (best == -1) {
        return true;
    }

    while (bm) {
        const unsigned int one = bm & (~bm + 1U);
        grid_[best] = firstDigit(one);
        ++stats_.decisions;
        if (search()) {
            return true;
        }
        grid_[best] = 0;
        ++stats_.backtracks;
        bm &= (bm - 1U);
    }

    return false;
}

BacktrackingSolveStats BacktrackingSolver::solve() {
    stats_ = BacktrackingSolveStats{};
    if (!validState()) {
        return stats_;
    }
    stats_.solved = search();
    return stats_;
}

static BacktrackingSolveStats solveWithBacktracking(const SudokuBoard& b, const std::vector<int>& initialGrid) {
    BacktrackingSolver solver(b.block_rows, b.block_cols, b.side_size, initialGrid);
    return solver.solve();
}

static SudokuBoard parseSudokuLine(const std::string& line) {
    SudokuBoard b;
    std::stringstream ss(line);
    std::string seg;
    std::vector<std::string> t;
    while (std::getline(ss, seg, ',')) {
        const std::string x = trim(seg);
        t.push_back(x);
    }
    if (t.size() < 4U) { b.error = "Za malo tokenow"; return b; }
    if (!parseLLStrict(t[0], b.seed)) { b.error = "Niepoprawny seed"; return b; }
    if (!parseIntStrict(t[1], b.block_rows) || !parseIntStrict(t[2], b.block_cols)) {
        b.error = "Niepoprawne Rows/Cols"; return b;
    }
    if (b.block_rows <= 0 || b.block_cols <= 0) { b.error = "Rows/Cols musza byc > 0"; return b; }
    b.side_size = b.block_rows * b.block_cols;
    b.total_cells = b.side_size * b.side_size;
    if (b.side_size <= 0 || b.side_size > 16) { b.error = "Nieobslugiwany rozmiar"; return b; }
    if (static_cast<int>(t.size()) < 3 + b.total_cells) { b.error = "Za malo danych"; return b; }

    b.cells.reserve(b.total_cells);
    for (int i = 0; i < b.total_cells; ++i) {
        const std::string tok = t[3 + i];
        Cell c;
        if (tok.empty() || tok == "0" || tok == "x" || tok == "X") {
            c.value = 0; c.revealed = false;
        } else if (tok[0] == 't' || tok[0] == 'T') {
            int v = 0;
            if (!parseIntStrict(tok.substr(1), v) || v < 1 || v > b.side_size) {
                b.error = "Niepoprawna wskazowka"; b.cells.clear(); return b;
            }
            c.value = v; c.revealed = true;
        } else {
            c.value = 0; c.revealed = false;
        }
        b.cells.push_back(c);
    }
    b.valid = true;
    return b;
}

static int countSolutionsWithBacktracking(const SudokuBoard& b, int limit) {
    std::vector<int> g(b.total_cells, 0);
    for (int i = 0; i < b.total_cells; ++i) if (b.cells[i].revealed) g[i] = b.cells[i].value;
    BacktrackingCounter bt(b.block_rows, b.block_cols, b.side_size, g);
    return bt.countSolutions(limit);
}

static std::string selectFolderModern() {
    std::string path;
    IFileOpenDialog* dialog = nullptr;
    HRESULT hr = CoCreateInstance(CLSID_FileOpenDialog, nullptr, CLSCTX_ALL,
                                  IID_IFileOpenDialog, reinterpret_cast<void**>(&dialog));
    if (FAILED(hr) || dialog == nullptr) return path;

    DWORD opts = 0;
    if (SUCCEEDED(dialog->GetOptions(&opts))) dialog->SetOptions(opts | FOS_PICKFOLDERS);
    dialog->SetTitle(L"Wybierz folder z plikami Sudoku (.txt)");

    if (SUCCEEDED(dialog->Show(nullptr))) {
        IShellItem* item = nullptr;
        if (SUCCEEDED(dialog->GetResult(&item)) && item != nullptr) {
            PWSTR selected = nullptr;
            if (SUCCEEDED(item->GetDisplayName(SIGDN_FILESYSPATH, &selected)) && selected != nullptr) {
                const std::wstring ws(selected);
                path.assign(ws.begin(), ws.end());
                CoTaskMemFree(selected);
            }
            item->Release();
        }
    }
    dialog->Release();
    return path;
}

static bool isTxtFile(const fs::path& p) {
    std::string e = p.extension().string();
    std::transform(e.begin(), e.end(), e.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return e == ".txt";
}

static bool isPathWithin(const fs::path& path, const fs::path& parent) {
    std::string fullPath = fs::absolute(path).lexically_normal().generic_string();
    std::string fullParent = fs::absolute(parent).lexically_normal().generic_string();

    std::transform(fullPath.begin(), fullPath.end(), fullPath.begin(),
                   [](unsigned char c) { return static_cast<char>(tolower(c)); });
    std::transform(fullParent.begin(), fullParent.end(), fullParent.begin(),
                   [](unsigned char c) { return static_cast<char>(tolower(c)); });

    if (!fullParent.empty() && fullParent.back() == '/') {
        fullParent.pop_back();
    }
    if (fullPath == fullParent) {
        return true;
    }
    fullParent.push_back('/');
    return fullPath.rfind(fullParent, 0) == 0;
}

static std::vector<fs::path> collectTxtFilesRecursive(const fs::path& root, const fs::path& excludedRoot) {
    std::vector<fs::path> files;
    std::error_code ec;

    fs::recursive_directory_iterator end;
    fs::recursive_directory_iterator it(root, fs::directory_options::skip_permission_denied, ec);
    while (it != end) {
        if (ec) {
            ec.clear();
            it.increment(ec);
            continue;
        }

        const fs::path currentPath = it->path();
        if (it->is_directory(ec)) {
            if (!ec && isPathWithin(currentPath, excludedRoot)) {
                it.disable_recursion_pending();
            }
            ec.clear();
            it.increment(ec);
            continue;
        }

        if (!it->is_regular_file(ec)) {
            ec.clear();
            it.increment(ec);
            continue;
        }
        ec.clear();

        if (isPathWithin(currentPath, excludedRoot)) {
            it.increment(ec);
            continue;
        }
        if (isTxtFile(currentPath)) {
            files.push_back(currentPath);
        }
        it.increment(ec);
    }

    std::sort(files.begin(), files.end());
    return files;
}

static long long countNonEmptyLinesInTxtFiles(const std::vector<fs::path>& files) {
    long long total = 0;
    for (const fs::path& p : files) {
        std::ifstream in(p);
        if (!in.is_open()) {
            continue;
        }

        std::string line;
        while (std::getline(in, line)) {
            if (!trim(line).empty()) {
                ++total;
            }
        }
    }
    return total;
}

static std::string folderKeyFromRelativePath(const fs::path& rel) {
    const std::string key = rel.generic_string();
    if (key.empty() || key == ".") {
        return "ROOT";
    }
    return key;
}

static std::string sanitizeFileName(const std::string& name) {
    std::string out = name;
    for (char& ch : out) {
        const bool bad = (ch == '<' || ch == '>' || ch == ':' || ch == '"' || ch == '/' ||
                          ch == '\\' || ch == '|' || ch == '?' || ch == '*');
        if (bad) ch = '_';
    }
    if (out.empty()) out = "ROOT";
    return out;
}

static std::string csvEscape(const std::string& field) {
    std::string out = "\"";
    for (char ch : field) {
        if (ch == '"') {
            out += "\"\"";
        } else {
            out.push_back(ch);
        }
    }
    out.push_back('"');
    return out;
}

static int difficultyLevelFromReport(const AnalysisReport& report) {
    if (report.contradiction) {
        return 0;
    }
    if (report.requires_guessing || report.hardest_rank >= 100) {
        return 10;
    }

    // Scale logical difficulty directly from hardest strategy rank to 1-10.
    if (report.hardest_rank <= 1) return 1;
    if (report.hardest_rank >= 10) return 10;
    return report.hardest_rank;
}

static std::string difficultyTypeFromReport(const AnalysisReport& report) {
    if (report.contradiction) {
        return "Sprzeczne";
    }
    if (report.requires_guessing || report.hardest_rank >= 100) {
        return "Eksperckie";
    }
    const int level = difficultyLevelFromReport(report);
    if (level <= 2) return "Latwe";
    if (level <= 4) return "Srednie";
    if (level <= 6) return "Trudne";
    if (level <= 8) return "Bardzo trudne";
    return "Eksperckie";
}

static std::string boardTypeFromBoard(const SudokuBoard& board) {
    if (board.side_size <= 0 || board.block_rows <= 0 || board.block_cols <= 0) {
        return "Nieznany";
    }
    std::ostringstream ss;
    ss << board.side_size << "x" << board.side_size
       << " (" << board.block_rows << "x" << board.block_cols << ")";
    return ss.str();
}

static void appendInvalidPuzzleReport(FolderStats& stats, const std::string& sourceFile, int lineNo,
                                      const std::string& parseError) {
    PuzzleReportEntry entry;
    entry.source_file = sourceFile;
    entry.line_no = lineNo;
    entry.valid = false;
    entry.sudoku_type = "Niepoprawne dane";
    entry.parse_error = parseError;
    entry.hardest_strategy = "Brak - blad danych";
    entry.strategy_usage.clear();
    entry.debug_logic_logs.clear();
    stats.puzzle_reports.push_back(std::move(entry));
}

static void appendValidPuzzleReport(FolderStats& stats, const std::string& sourceFile, int lineNo,
                                    const SudokuBoard& board, const AnalysisReport& report) {
    PuzzleReportEntry entry;
    entry.source_file = sourceFile;
    entry.line_no = lineNo;
    entry.valid = true;
    entry.sudoku_type = difficultyTypeFromReport(report);
    entry.board_type = boardTypeFromBoard(board);
    entry.initial_clues = report.initial_clues;
    entry.difficulty_level = difficultyLevelFromReport(report);
    entry.solved_logically = report.solved_logically;
    entry.requires_guessing = report.requires_guessing;
    entry.solved_with_backtracking = report.solved_with_backtracking;
    entry.contradiction = report.contradiction;
    entry.solution_count = report.solution_count;
    entry.backtracking_nodes = report.backtracking_nodes;
    entry.backtracking_decisions = report.backtracking_decisions;
    entry.backtracking_backtracks = report.backtracking_backtracks;
    entry.strategy_usage = report.strategy_usage;
    entry.hardest_strategy = report.hardest_strategy;
    entry.debug_logic_logs = report.debug_logic_logs;
    stats.puzzle_reports.push_back(std::move(entry));
}

static void updateFolderStats(FolderStats& stats, const AnalysisReport& report) {
    ++stats.analyzed_puzzles;
    stats.clues_sum += report.initial_clues;

    if (report.contradiction) ++stats.contradictions;
    if (report.solved_logically) ++stats.solved_logically;
    if (report.requires_guessing) ++stats.requires_guessing;
    if (report.solved_with_backtracking) ++stats.solved_with_backtracking;

    if (report.solution_count == 1) ++stats.unique_solutions;
    else if (report.solution_count == 0) ++stats.no_solution;
    else ++stats.multiple_solutions;

    stats.backtracking_nodes_sum += report.backtracking_nodes;
    stats.backtracking_decisions_sum += report.backtracking_decisions;
    stats.backtracking_backtracks_sum += report.backtracking_backtracks;

    for (const auto& it : report.strategy_usage) {
        stats.strategy_usage[it.first] += it.second;
    }
    stats.hardest_histogram[report.hardest_strategy] += 1;

    const int difficulty = difficultyLevelFromReport(report);
    if (difficulty > 0) {
        stats.difficulty_sum += difficulty;
        ++stats.difficulty_count;
        if (difficulty > stats.max_difficulty) {
            stats.max_difficulty = difficulty;
        }
    }

    if (report.hardest_rank > stats.hardest_rank_seen) {
        stats.hardest_rank_seen = report.hardest_rank;
        stats.hardest_name_seen = report.hardest_strategy;
    }
}

static void writeFolderReport(const fs::path& outDir, const std::string& folderKey, const FolderStats& stats) {
    fs::path relFolder = stats.relative_folder;
    if (relFolder.empty() || relFolder == fs::path(".")) {
        relFolder = fs::path("ROOT");
    }
    const fs::path folderOutDir = outDir / relFolder;
    std::error_code ec;
    fs::create_directories(folderOutDir, ec);
    const fs::path outPath = folderOutDir / "statystyki_folder.txt";
    std::ofstream out(outPath);
    if (!out.is_open()) {
        std::cerr << "Nie mozna zapisac raportu folderu: " << outPath.string() << "\n";
        return;
    }

    out << "=== Statystyki Folderu Sudoku ===\n";
    out << "Folder: " << folderKey << "\n";
    out << "Niepuste wpisy (linie): " << stats.non_empty_lines << "\n";
    out << "Poprawnie przeanalizowane plansze: " << stats.analyzed_puzzles << "\n";
    out << "Bledne wpisy: " << stats.invalid_lines << "\n\n";

    out << "Wynik solvera:\n";
    out << "- Rozwiazane logicznie: " << stats.solved_logically << "\n";
    out << "- Wymaga zgadywania/backtrackingu: " << stats.requires_guessing << "\n";
    out << "- Rozwiazane z uzyciem backtrackingu: " << stats.solved_with_backtracking << "\n";
    out << "- Sprzeczne plansze: " << stats.contradictions << "\n\n";

    out << "Metryki backtrackingu (suma):\n";
    out << "- Decyzje: " << stats.backtracking_decisions_sum << "\n";
    out << "- Backtracki: " << stats.backtracking_backtracks_sum << "\n";
    out << "- Nody rekursji: " << stats.backtracking_nodes_sum << "\n\n";

    out << "Unikalnosc rozwiazania:\n";
    out << "- Unikalne: " << stats.unique_solutions << "\n";
    out << "- Wiele rozwiazan: " << stats.multiple_solutions << "\n";
    out << "- Brak rozwiazania: " << stats.no_solution << "\n\n";

    out << "Srednia liczba clues: ";
    if (stats.analyzed_puzzles > 0) {
        const double avg = static_cast<double>(stats.clues_sum) / static_cast<double>(stats.analyzed_puzzles);
        out << std::fixed << std::setprecision(2) << avg << "\n";
    } else {
        out << "0.00\n";
    }

    out << "Poziom trudnosci sudoku (1-10): ";
    if (stats.difficulty_count > 0) {
        const double avgDifficulty = static_cast<double>(stats.difficulty_sum) / static_cast<double>(stats.difficulty_count);
        out << std::fixed << std::setprecision(2) << avgDifficulty
            << " (max: " << stats.max_difficulty << ")\n";
    } else {
        out << "brak danych\n";
    }

    out << "\nUzycie strategii (suma):\n";
    const std::vector<Strategy> order = {
        Strategy::NakedSingle, Strategy::HiddenSingle,
        Strategy::NakedPair, Strategy::HiddenPair,
        Strategy::PointingPairsTriples, Strategy::BoxLineReduction,
        Strategy::NakedTriple, Strategy::HiddenTriple,
        Strategy::NakedQuad, Strategy::HiddenQuad,
        Strategy::XWing, Strategy::YWing, Strategy::XYZWing, Strategy::WXYZWing, Strategy::Swordfish, Strategy::Jellyfish, Strategy::FrankenMutantFish, Strategy::KrakenFish,
        Strategy::Skyscraper, Strategy::TwoStringKite, Strategy::SimpleColoring, Strategy::ThreeDMedusa,
        Strategy::FinnedXWingSashimi, Strategy::FinnedSwordfish, Strategy::FinnedJellyfish, Strategy::EmptyRectangle,
        Strategy::UniqueRectangle, Strategy::UniqueLoop, Strategy::BivalueOddagon, Strategy::AvoidableRectangle, Strategy::BUGPlus1,
        Strategy::RemotePairs, Strategy::WWing, Strategy::GroupedXCycle, Strategy::XChain, Strategy::XYChain,
        Strategy::GroupedAIC, Strategy::AIC, Strategy::ContinuousNiceLoop,
        Strategy::ALSXZ, Strategy::ALSXYWing, Strategy::ALSChain, Strategy::DeathBlossom,
        Strategy::SueDeCoq, Strategy::MSLS, Strategy::Exocet, Strategy::SeniorExocet, Strategy::SKLoop, Strategy::PatternOverlayMethod,
        Strategy::ForcingChains, Strategy::Backtracking
    };
    bool any = false;
    for (Strategy s : order) {
        const auto it = stats.strategy_usage.find(s);
        if (it == stats.strategy_usage.end() || it->second <= 0) continue;
        any = true;
        out << "- " << strategyName(s) << ": " << it->second << "\n";
    }
    if (!any) {
        out << "- Brak\n";
    }

    out << "\nStatus implementacji technik:\n";
    for (Strategy s : order) {
        out << "- " << strategyName(s) << ": " << strategyImplementationStatus(s) << "\n";
    }

    out << "\nNajtrudniejsza technika zaobserwowana w folderze: " << stats.hardest_name_seen << "\n";
    out << "Rozklad najtrudniejszej techniki (na plansze):\n";
    for (const auto& it : stats.hardest_histogram) {
        out << "- " << it.first << ": " << it.second << "\n";
    }

    out << "\n=== Raporty Sudoku w folderze (zbiorczo) ===\n";
    if (stats.puzzle_reports.empty()) {
        out << "Brak wpisow Sudoku.\n";
    }

    for (std::size_t i = 0; i < stats.puzzle_reports.size(); ++i) {
        const PuzzleReportEntry& entry = stats.puzzle_reports[i];
        out << "\n[" << (i + 1) << "] Plik: " << entry.source_file
            << " | Linia: " << entry.line_no << "\n";

        if (!entry.valid) {
            out << "Typ: " << entry.sudoku_type << "\n";
            out << "Blad: " << entry.parse_error << "\n";
            continue;
        }

        out << "Typ: " << entry.sudoku_type << "\n";
        out << "Rozmiar: " << entry.board_type << "\n";
        out << "Najtrudniejsza technika: " << entry.hardest_strategy << "\n";
        out << "Poziom trudnosci (1-10): " << entry.difficulty_level << "\n";
        out << "Liczba clues: " << entry.initial_clues << "\n";
        out << "Rozwiazane logicznie: " << (entry.solved_logically ? "TAK" : "NIE") << "\n";
        out << "Wymaga zgadywania: " << (entry.requires_guessing ? "TAK" : "NIE") << "\n";
        out << "Rozwiazane backtrackingiem: " << (entry.solved_with_backtracking ? "TAK" : "NIE") << "\n";
        out << "Sprzeczne: " << (entry.contradiction ? "TAK" : "NIE") << "\n";
        out << "Liczba rozwiazan: " << entry.solution_count << "\n";
        out << "Backtracking - decyzje: " << entry.backtracking_decisions << "\n";
        out << "Backtracking - backtracki: " << entry.backtracking_backtracks << "\n";
        out << "Backtracking - nody rekursji: " << entry.backtracking_nodes << "\n";
        out << "Uzycie metod:\n";
        for (Strategy s : order) {
            const auto it = entry.strategy_usage.find(s);
            const long long used = (it == entry.strategy_usage.end()) ? 0LL : static_cast<long long>(it->second);
            out << "  - " << strategyName(s) << ": " << used << "\n";
        }
        out << "Debug log (strategie zaawansowane):\n";
        if (entry.debug_logic_logs.empty()) {
            out << "  - Brak wpisow\n";
        } else {
            for (const std::string& line : entry.debug_logic_logs) {
                out << "  - " << line << "\n";
            }
        }
    }

    auto writeFilteredList = [&](const std::string& fileName, const std::string& title, auto predicate) {
        const fs::path listPath = folderOutDir / fileName;
        std::ofstream listOut(listPath);
        if (!listOut.is_open()) {
            std::cerr << "Nie mozna zapisac raportu: " << listPath.string() << "\n";
            return;
        }

        listOut << title << "\n";
        listOut << "Folder: " << folderKey << "\n\n";

        bool anyMatch = false;
        for (const PuzzleReportEntry& entry : stats.puzzle_reports) {
            if (!entry.valid) continue;
            if (!predicate(entry)) continue;
            anyMatch = true;
            listOut << entry.source_file << " | Linia: " << entry.line_no << "\n";
        }

        if (!anyMatch) {
            listOut << "Brak wpisow.\n";
        }
    };

    writeFilteredList(
        "sprzeczne_tak.txt",
        "=== Sudoku ze statusem: Sprzeczne = TAK ===",
        [](const PuzzleReportEntry& entry) { return entry.contradiction; }
    );
    writeFilteredList(
        "liczba_rozwiazan_wiecej_niz_1.txt",
        "=== Sudoku ze statusem: Liczba rozwiazan > 1 ===",
        [](const PuzzleReportEntry& entry) { return entry.solution_count > 1; }
    );
    writeFilteredList(
        "wymaga_zgadywania_tak.txt",
        "=== Sudoku ze statusem: Wymaga zgadywania = TAK ===",
        [](const PuzzleReportEntry& entry) { return entry.requires_guessing; }
    );
}

static void writeGlobalSummary(const fs::path& outDir, const std::map<std::string, FolderStats>& allStats,
                               long long txtFilesScanned) {
    const fs::path outPath = outDir / "podsumowanie_folderow.txt";
    std::ofstream out(outPath);
    if (!out.is_open()) {
        std::cerr << "Nie mozna zapisac podsumowania globalnego: " << outPath.string() << "\n";
        return;
    }

    long long allNonEmpty = 0;
    long long allInvalid = 0;
    long long allAnalyzed = 0;
    long long allSolved = 0;
    long long allGuess = 0;
    long long allSolvedWithBacktracking = 0;
    long long allContradictions = 0;
    long long allUnique = 0;
    long long allMultiple = 0;
    long long allNoSolution = 0;
    long long allBacktrackingNodes = 0;
    long long allBacktrackingDecisions = 0;
    long long allBacktrackingBacktracks = 0;
    long long allDifficultySum = 0;
    long long allDifficultyCount = 0;
    int allMaxDifficulty = 0;

    for (const auto& entry : allStats) {
        const FolderStats& st = entry.second;
        allNonEmpty += st.non_empty_lines;
        allInvalid += st.invalid_lines;
        allAnalyzed += st.analyzed_puzzles;
        allSolved += st.solved_logically;
        allGuess += st.requires_guessing;
        allSolvedWithBacktracking += st.solved_with_backtracking;
        allContradictions += st.contradictions;
        allUnique += st.unique_solutions;
        allMultiple += st.multiple_solutions;
        allNoSolution += st.no_solution;
        allBacktrackingNodes += st.backtracking_nodes_sum;
        allBacktrackingDecisions += st.backtracking_decisions_sum;
        allBacktrackingBacktracks += st.backtracking_backtracks_sum;
        allDifficultySum += st.difficulty_sum;
        allDifficultyCount += st.difficulty_count;
        if (st.max_difficulty > allMaxDifficulty) {
            allMaxDifficulty = st.max_difficulty;
        }
    }

    out << "=== Podsumowanie Globalne ===\n";
    out << "Liczba folderow: " << allStats.size() << "\n";
    out << "Przeskanowane pliki .txt: " << txtFilesScanned << "\n";
    out << "Niepuste wpisy (linie): " << allNonEmpty << "\n";
    out << "Poprawnie przeanalizowane plansze: " << allAnalyzed << "\n";
    out << "Bledne wpisy: " << allInvalid << "\n\n";
    out << "Rozwiazane logicznie: " << allSolved << "\n";
    out << "Wymaga zgadywania/backtrackingu: " << allGuess << "\n";
    out << "Rozwiazane z uzyciem backtrackingu: " << allSolvedWithBacktracking << "\n";
    out << "Sprzeczne plansze: " << allContradictions << "\n";
    out << "Unikalne rozwiazania: " << allUnique << "\n";
    out << "Wiele rozwiazan: " << allMultiple << "\n";
    out << "Brak rozwiazania: " << allNoSolution << "\n";
    out << "Backtracking - decyzje: " << allBacktrackingDecisions << "\n";
    out << "Backtracking - backtracki: " << allBacktrackingBacktracks << "\n";
    out << "Backtracking - nody rekursji: " << allBacktrackingNodes << "\n\n";
    out << "Poziom trudnosci sudoku (1-10): ";
    if (allDifficultyCount > 0) {
        const double avgDifficulty = static_cast<double>(allDifficultySum) / static_cast<double>(allDifficultyCount);
        out << std::fixed << std::setprecision(2) << avgDifficulty
            << " (max: " << allMaxDifficulty << ")\n\n";
    } else {
        out << "brak danych\n\n";
    }

    out << "=== Szczegoly per folder ===\n";
    for (const auto& entry : allStats) {
        const FolderStats& st = entry.second;
        const double avgDifficulty = (st.difficulty_count > 0)
            ? static_cast<double>(st.difficulty_sum) / static_cast<double>(st.difficulty_count)
            : 0.0;
        out << "- " << entry.first
            << " | analyzed=" << st.analyzed_puzzles
            << " | invalid=" << st.invalid_lines
            << " | logic=" << st.solved_logically
            << " | guess=" << st.requires_guessing
            << " | bt_solved=" << st.solved_with_backtracking
            << " | unique=" << st.unique_solutions
            << " | diff=" << std::fixed << std::setprecision(2) << avgDifficulty
            << " | hardest=" << st.hardest_name_seen
            << "\n";
    }
}

static void writeFolderCsv(const fs::path& outDir, const std::map<std::string, FolderStats>& allStats) {
    const fs::path csvPath = outDir / "statystyki_folderow.csv";
    std::ofstream out(csvPath);
    if (!out.is_open()) {
        std::cerr << "Nie mozna zapisac CSV: " << csvPath.string() << "\n";
        return;
    }

    out << "folder,non_empty_lines,analyzed_puzzles,invalid_lines,solved_logically,requires_guessing,"
        << "solved_with_backtracking,contradictions,unique_solutions,multiple_solutions,no_solution,"
        << "backtracking_decisions,backtracking_backtracks,backtracking_nodes,avg_clues,hardest_seen,"
        << "avg_difficulty,max_difficulty,"
        << "naked_single,hidden_single,naked_pair,hidden_pair,pointing_pairs_triples,box_line_reduction,"
        << "naked_triple,hidden_triple,naked_quad,hidden_quad,"
        << "x_wing,y_wing,xyz_wing,wxyz_wing,swordfish,jellyfish,franken_mutant_fish,kraken_fish,skyscraper,two_string_kite,simple_coloring,three_d_medusa,"
        << "finned_x_wing_sashimi,finned_swordfish,finned_jellyfish,empty_rectangle,unique_rectangle,unique_loop,bivalue_oddagon,avoidable_rectangle,bug_plus_1,"
        << "remote_pairs,w_wing,grouped_x_cycle,x_chain,xy_chain,grouped_aic,aic,continuous_nice_loop,"
        << "als_xz,als_xy_wing,als_chain,death_blossom,sue_de_coq,msls,exocet,senior_exocet,sk_loop,pattern_overlay_method,forcing_chains,backtracking\n";

    for (const auto& entry : allStats) {
        const std::string& folder = entry.first;
        const FolderStats& st = entry.second;
        const double avgClues = (st.analyzed_puzzles > 0)
            ? static_cast<double>(st.clues_sum) / static_cast<double>(st.analyzed_puzzles)
            : 0.0;
        const double avgDifficulty = (st.difficulty_count > 0)
            ? static_cast<double>(st.difficulty_sum) / static_cast<double>(st.difficulty_count)
            : 0.0;

        auto countUsage = [&](Strategy s) -> long long {
            const auto it = st.strategy_usage.find(s);
            return (it == st.strategy_usage.end()) ? 0LL : it->second;
        };

        out << csvEscape(folder) << ","
            << st.non_empty_lines << ","
            << st.analyzed_puzzles << ","
            << st.invalid_lines << ","
            << st.solved_logically << ","
            << st.requires_guessing << ","
            << st.solved_with_backtracking << ","
            << st.contradictions << ","
            << st.unique_solutions << ","
            << st.multiple_solutions << ","
            << st.no_solution << ","
            << st.backtracking_decisions_sum << ","
            << st.backtracking_backtracks_sum << ","
            << st.backtracking_nodes_sum << ","
            << std::fixed << std::setprecision(2) << avgClues << ","
            << csvEscape(st.hardest_name_seen) << ","
            << std::fixed << std::setprecision(2) << avgDifficulty << ","
            << st.max_difficulty << ","
            << countUsage(Strategy::NakedSingle) << ","
            << countUsage(Strategy::HiddenSingle) << ","
            << countUsage(Strategy::NakedPair) << ","
            << countUsage(Strategy::HiddenPair) << ","
            << countUsage(Strategy::PointingPairsTriples) << ","
            << countUsage(Strategy::BoxLineReduction) << ","
            << countUsage(Strategy::NakedTriple) << ","
            << countUsage(Strategy::HiddenTriple) << ","
            << countUsage(Strategy::NakedQuad) << ","
            << countUsage(Strategy::HiddenQuad) << ","
            << countUsage(Strategy::XWing) << ","
            << countUsage(Strategy::YWing) << ","
            << countUsage(Strategy::XYZWing) << ","
            << countUsage(Strategy::WXYZWing) << ","
            << countUsage(Strategy::Swordfish) << ","
            << countUsage(Strategy::Jellyfish) << ","
            << countUsage(Strategy::FrankenMutantFish) << ","
            << countUsage(Strategy::KrakenFish) << ","
            << countUsage(Strategy::Skyscraper) << ","
            << countUsage(Strategy::TwoStringKite) << ","
            << countUsage(Strategy::SimpleColoring) << ","
            << countUsage(Strategy::ThreeDMedusa) << ","
            << countUsage(Strategy::FinnedXWingSashimi) << ","
            << countUsage(Strategy::FinnedSwordfish) << ","
            << countUsage(Strategy::FinnedJellyfish) << ","
            << countUsage(Strategy::EmptyRectangle) << ","
            << countUsage(Strategy::UniqueRectangle) << ","
            << countUsage(Strategy::UniqueLoop) << ","
            << countUsage(Strategy::BivalueOddagon) << ","
            << countUsage(Strategy::AvoidableRectangle) << ","
            << countUsage(Strategy::BUGPlus1) << ","
            << countUsage(Strategy::RemotePairs) << ","
            << countUsage(Strategy::WWing) << ","
            << countUsage(Strategy::GroupedXCycle) << ","
            << countUsage(Strategy::XChain) << ","
            << countUsage(Strategy::XYChain) << ","
            << countUsage(Strategy::GroupedAIC) << ","
            << countUsage(Strategy::AIC) << ","
            << countUsage(Strategy::ContinuousNiceLoop) << ","
            << countUsage(Strategy::ALSXZ) << ","
            << countUsage(Strategy::ALSXYWing) << ","
            << countUsage(Strategy::ALSChain) << ","
            << countUsage(Strategy::DeathBlossom) << ","
            << countUsage(Strategy::SueDeCoq) << ","
            << countUsage(Strategy::MSLS) << ","
            << countUsage(Strategy::Exocet) << ","
            << countUsage(Strategy::SeniorExocet) << ","
            << countUsage(Strategy::SKLoop) << ","
            << countUsage(Strategy::PatternOverlayMethod) << ","
            << countUsage(Strategy::ForcingChains) << ","
            << countUsage(Strategy::Backtracking)
            << "\n";
    }
}

int main() {
    if (FAILED(CoInitialize(nullptr))) {
        MessageBoxA(nullptr, "Nie udalo sie zainicjalizowac COM.", "Blad", MB_ICONERROR | MB_OK);
        return 1;
    }

    IProgressDialog* progress = nullptr;
    try {
        const std::string input = selectFolderModern();
        if (input.empty()) {
            MessageBoxA(nullptr, "Anulowano wybor folderu.", "Koniec", MB_OK);
            CoUninitialize();
            return 0;
        }

        const fs::path inputRoot = fs::absolute(fs::path(input)).lexically_normal();
        const fs::path outDir = fs::current_path() / "analizy_folderow";
        const fs::path folderReportsDir = outDir / "raporty_folderow";
        const fs::path outDirAbs = fs::absolute(outDir).lexically_normal();
        fs::create_directories(outDir);
        std::error_code cleanupEc;
        for (fs::directory_iterator it(outDir, cleanupEc), end; !cleanupEc && it != end; it.increment(cleanupEc)) {
            if (!it->is_regular_file()) {
                continue;
            }
            const std::string name = it->path().filename().string();
            if (name.rfind("statystyki_folder_", 0) == 0 && it->path().extension() == ".txt") {
                std::error_code removeEc;
                fs::remove(it->path(), removeEc);
            }
        }
        fs::remove_all(folderReportsDir);
        fs::create_directories(folderReportsDir);

        const std::vector<fs::path> txtFiles = collectTxtFilesRecursive(inputRoot, outDirAbs);
        const long long txtFilesScanned = static_cast<long long>(txtFiles.size());
        const long long totalEntries = countNonEmptyLinesInTxtFiles(txtFiles);
        const ULONGLONG progressTotal = static_cast<ULONGLONG>(totalEntries > 0 ? totalEntries : 1);

        if (SUCCEEDED(CoCreateInstance(
                CLSID_ProgressDialog,
                nullptr,
                CLSCTX_INPROC_SERVER,
                IID_IProgressDialog,
                reinterpret_cast<void**>(&progress))) && progress != nullptr) {
            progress->SetTitle(L"Analiza Sudoku");
            const std::wstring line1 = L"Pliki .txt: " + std::to_wstring(txtFilesScanned);
            const std::wstring line2 = L"Do analizy (linie Sudoku): " + std::to_wstring(totalEntries);
            progress->SetLine(1, line1.c_str(), FALSE, nullptr);
            progress->SetLine(2, line2.c_str(), FALSE, nullptr);
            progress->SetLine(3, L"Start analizy...", FALSE, nullptr);
            progress->SetCancelMsg(L"Przerywanie analizy...", nullptr);
            progress->StartProgressDialog(nullptr, nullptr, PROGDLG_NORMAL | PROGDLG_AUTOTIME, nullptr);
            progress->SetProgress64(0ULL, progressTotal);
        }

        std::map<std::string, FolderStats> statsByFolder;
        long long totalInvalid = 0;
        long long processedEntries = 0;
        bool cancelled = false;
        long long totalAnalyzed = 0;

        for (const fs::path& currentPath : txtFiles) {
            std::ifstream in(currentPath);
            if (!in.is_open()) {
                continue;
            }
            const std::string sourceFileName = currentPath.filename().string();

            std::error_code ec;
            fs::path relFolder = fs::relative(currentPath.parent_path(), inputRoot, ec);
            if (ec || relFolder.empty()) {
                relFolder = ".";
            }
            const std::string folderKey = folderKeyFromRelativePath(relFolder);
            FolderStats& folderStats = statsByFolder[folderKey];
            folderStats.relative_folder = relFolder;

            std::string line;
            int line_no = 0;
            while (std::getline(in, line)) {
                ++line_no;
                const std::string clean = trim(line);
                if (clean.empty()) {
                    continue;
                }

                ++processedEntries;
                if (progress != nullptr && (processedEntries == 1 || processedEntries % 25 == 0 || processedEntries == totalEntries)) {
                    const std::wstring p1 = L"Przeanalizowano Sudoku: " + std::to_wstring(processedEntries)
                        + L" / " + std::to_wstring(totalEntries);
                    const std::wstring p2 = L"Plik: " + currentPath.filename().wstring();
                    progress->SetLine(1, p1.c_str(), FALSE, nullptr);
                    progress->SetLine(2, p2.c_str(), FALSE, nullptr);
                    progress->SetProgress64(static_cast<ULONGLONG>(processedEntries), progressTotal);
                }
                if (progress != nullptr && progress->HasUserCancelled()) {
                    cancelled = true;
                    break;
                }

                ++folderStats.non_empty_lines;
                SudokuBoard board = parseSudokuLine(clean);
                if (!board.valid) {
                    ++folderStats.invalid_lines;
                    ++totalInvalid;
                    appendInvalidPuzzleReport(folderStats, sourceFileName, line_no, board.error);
                    std::cerr << "Pominieto wpis (" << currentPath.string() << ":" << line_no
                              << "): " << board.error << "\n";
                    continue;
                }

                SudokuAnalyzer analyzer(board);
                AnalysisReport report = analyzer.run();
                report.solution_count = countSolutionsWithBacktracking(board, 2);
                report.unique_solution = (report.solution_count == 1);
                updateFolderStats(folderStats, report);
                appendValidPuzzleReport(folderStats, sourceFileName, line_no, board, report);
                ++totalAnalyzed;
            }

            if (cancelled) {
                break;
            }
        }

        if (progress != nullptr) {
            progress->SetProgress64(static_cast<ULONGLONG>(processedEntries), progressTotal);
            progress->StopProgressDialog();
            progress->Release();
            progress = nullptr;
        }

        for (const auto& entry : statsByFolder) {
            writeFolderReport(folderReportsDir, entry.first, entry.second);
        }
        writeGlobalSummary(outDir, statsByFolder, txtFilesScanned);
        writeFolderCsv(outDir, statsByFolder);

        std::ostringstream msg;
        if (cancelled) {
            msg << "Analiza przerwana przez uzytkownika.\n";
        }
        msg << "Przeskanowane pliki .txt: " << txtFilesScanned
            << "\nPrzetworzone wpisy Sudoku: " << processedEntries
            << "\nPoprawnie przeanalizowane plansze: " << totalAnalyzed
            << "\nFoldery ze statystykami: " << statsByFolder.size()
            << "\nBledne wpisy: " << totalInvalid
            << "\nRaporty zapisano w: " << outDir.string()
            << "\nRaporty folderowe: " << folderReportsDir.string()
            << "\nCSV: " << (outDir / "statystyki_folderow.csv").string();
        MessageBoxA(nullptr, msg.str().c_str(), "Analiza Sudoku zakonczona", MB_OK);
    } catch (const std::exception& ex) {
        if (progress != nullptr) {
            progress->StopProgressDialog();
            progress->Release();
            progress = nullptr;
        }
        MessageBoxA(nullptr, ex.what(), "Blad krytyczny", MB_ICONERROR | MB_OK);
        CoUninitialize();
        return 1;
    }

    CoUninitialize();
    return 0;
}
