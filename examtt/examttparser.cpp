#include "examttparser.h"
#include <ctime>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <set>
#include <string>
#include <iomanip>
#include <map>

#include "../emilibase.h"

#include "../SA/sa_exploration.h"

/**
 * @brief Useless
 * in place implementation
 */
struct SANoExploration : SAExploration {
    SANoExploration() : SAExploration(
        new emili::EmptyNeighBorHood(),
        new SABasicAcceptance(),
        new SAWhileTrueTermination(),
        "no"
    ){

    }

    virtual emili::Solution* nextSolution(emili::Solution *startingSolution, SAStatus&) override {
        return startingSolution;
    }
};

namespace prs {
namespace ExamTT {

const std::string

// Algos
    IG = "ig",
    ILS = "ils",
    TABU = "tabu",
    FIRST = "first",
    BEST = "best",
    VND = "vnd",
    SA = "sa",
    TEST_INIT = "stin",

// initial solution heuristics
    INITIAL_RANDOM = "random",

// Termination criteria
    TERMINATION_MAXSTEPS = "maxstep",
    TERMINATION_TIME = "time",
    TERMINATION_LOCMIN = "locmin",
    TERMINATION_ITERA = "iteration",
    TERMINATION_WTRUE = "true",
    TERMINATION_SOA = "soater",

// acceptance criteria
    ACCEPTANCE_PROB = "prob",
    ACCEPTANCE_METRO = "metropolis",
    ACCEPTANCE_PMETRO = "pmetro",
    ACCEPTANCE_IMPROVE_PLATEAU = "implat",
    ACCEPTANCE_TEST = "testacc",
    ACCEPTANCE_SOA = "soaacc",
    ACCEPTANCE_ALWAYS = "always",
    ACCEPTANCE_INTENSIFY = "intensify",
    ACCEPTANCE_DIVERSIFY = "diversify",
    ACCEPTANCE_IMPROVE = "improve",
    ACCEPTANCE_SA_METRO = "sa_metropolis",
    ACCEPTANCE_SA = "saacc";

std::set<std::string> PROBLEMS_DEF = {
    "ExamTT"
};

using std::string;
using std::vector;
using std::get;
using std::setw;

std::string ExamTTParser::info()
{
    ostringstream oss;
    oss << "\nUsage: ";
    oss << "EMILI INSTANCE_FILE_PATH ExamTT <LOCAL_SEARCH | ITERATED_LOCAL_SEARCH | TABU_SEARCH | VND_SEARCH> [rnds seed]";
    return oss.str();
}

void ExamTTParser::genericError(string name) {
    cerr << "ERROR: " << name << endl;
    cout << info() << endl;
    exit(-1);
}

void ExamTTParser::genericError(std::ostream& stream) {
    stream << endl << "ERROR" << endl;
    cout << info() << endl;
    exit(-1);
}

void ExamTTParser::genericError(std::ostringstream& stream) {
    cerr << endl << "ERROR: " << stream.str() << endl;
    cout << info() << endl;
    exit(-1);
}

void ExamTTParser::errorExpected(prs::TokenManager& tm, string name, const std::vector<string> &tokens)
{
    cerr << "'" << *tm << "' -> ERROR a " << name << " is expected : ";
    if(tokens.size() == 0)
        cerr << "<>";
    else {
        auto it = tokens.begin();
        cerr << "<" << *it++;
        while(it != tokens.end())
            cerr << " | " << *it++;
        cerr << ">";
    }
    cerr << endl;

    cout << info() << endl;
    exit(-1);
}

namespace {
bool checkTokenParams(prs::TokenManager& tm, std::string val, std::vector<std::string> const& params) {
    if(tm.peek() == val) {
        std::ostringstream oss;
        oss << val << " (" << params.size() << ") ";
        for(auto const& p : params)
            oss << p << " ";
        printTab(oss.str());
        tm.next();
        return true;
    } else {
        return false;
    }
}

int getIntParam(prs::TokenManager& tm, std::string val) {
    int i = tm.getInteger();
    std::ostringstream oss;
    oss << val << ":" << i << endl;
    printTab(oss.str());
    return i;
}

float getDecimalParam(prs::TokenManager& tm, std::string val) {
    float i = tm.getDecimal();
    std::ostringstream oss;
    oss << val << ":" << i << endl;
    printTab(oss.str());
    return i;
}

int getIntPercentParam(prs::TokenManager& tm, std::string val, int max) {
    int i;
    if(tm.checkToken("percent")) {
        i = tm.getInteger() * max / 100;
    } else {
        i = tm.getInteger();
    }
    std::ostringstream oss;
    oss << val << ":" << i;
    printTab(oss.str());
    return i;
}

float getFloatPercentParam(prs::TokenManager& tm, std::string val, float max) {
    float i;
    if(tm.checkToken("percent")) {
        i = tm.getDecimal() * max / 100;
    } else {
        i = tm.getDecimal();
    }
    std::ostringstream oss;
    oss << val << ":" << i << endl;
    printTab(oss.str());
    return i;
}

}

struct ArgParser {
    std::map<std::string, int> dataInt;
    std::map<std::string, bool> dataBool;
    std::map<std::string, float> dataFloat;
    std::vector<std::string> order;
    std::set<std::string> given;

    // add

    void addInt(std::string s, int d) {
        if(dataInt.count(s) || dataBool.count(s) || dataFloat.count(s))
            throw std::invalid_argument(" argument " + s + " already exist !");
        dataInt[s] = d;
        order.push_back(s);
    }

    void addFloat(std::string s, float d) {
        if(dataInt.count(s) || dataBool.count(s) || dataFloat.count(s))
            throw std::invalid_argument(" argument " + s + " already exist !");
        dataFloat[s] = d;
        order.push_back(s);
    }

    void addBool(std::string s, bool d) {
        if(dataInt.count(s) || dataBool.count(s) || dataFloat.count(s))
            throw std::invalid_argument(" argument " + s + " already exist !");
        dataBool[s] = d;
        order.push_back(s);
    }

    // get

    int Int(std::string s) {
        if(! dataInt.count(s))
            throw std::invalid_argument("Missing int " + s);
        return dataInt[s];
    }

    float Float(std::string s) {
        if(! dataFloat.count(s))
            throw std::invalid_argument("Missing int " + s);
        return dataFloat[s];
    }

    bool Bool(std::string s) {
        if(! dataBool.count(s))
            throw std::invalid_argument("Missing bool " + s);
        return dataBool[s];
    }

    bool has(std::string s) {
        if(!(dataInt.count(s) || dataBool.count(s) || dataFloat.count(s)))
            throw std::invalid_argument("Argument " + s + " does not exist");
        return given.count(s);
    }

    // methods

    void operator()(prs::TokenManager& tm) {
        parse(tm);
    }

    void parse(prs::TokenManager& tm) {
        /*
         * Simple Algo : for(;;) if(tm.checkToken("A")) A = tm.getInteger() else if()... else break;
         */

        for(;;) {
            bool found = false;

            for(auto& p : dataInt) {
                if(tm.checkToken(p.first)) {
                    p.second = tm.getInteger();
                    found = true;
                    given.insert(p.first);
                }
            }

            for(auto& p : dataFloat) {
                if(tm.checkToken(p.first)) {
                    p.second = tm.getDecimal();
                    found = true;
                    given.insert(p.first);
                }
            }

            for(auto& p : dataBool) {
                if(tm.checkToken(p.first)) {
                    p.second = true;
                    found = true;
                    given.insert(p.first);
                }
                else if(tm.checkToken("no-" + p.first)) {
                    p.second = false;
                    found = true;
                    given.insert(p.first);
                }
            }

            if(! found)
                break;
        }
    }

    /**
     * @brief Sort by type, then name
     */
    void printTypeAlpha() {
        for(auto& p : dataInt)
            cout << setw(20) << p.first << ": " << p.second << endl;
        for(auto& p : dataFloat)
            cout << setw(20) << p.first << ": " << p.second << endl;
        for(auto& p : dataBool)
            cout << setw(20) << p.first << ": " << boolalpha << p.second << endl;
    }

    /**
     * @brief Insertion order
     */
    void print() {
        for(auto x : order)
            if(dataInt.count(x))
                cout << setw(20) << x << ": " << dataInt[x] << endl;
            else if(dataFloat.count(x))
                cout << setw(20) << x << ": " << dataInt[x] << endl;
            else if(dataBool.count(x))
                cout << setw(20) << x << ": " << boolalpha << dataBool[x] << endl;
    }

};
emili::LocalSearch* ExamTTParser::search(prs::TokenManager& tm)
{
    static const std::string SA_BSU = "sa_bsu";
    static const std::string INTERACTIVE = "interactive";
    static const std::string TEST_KEMPE = "test_kempe";
    static const std::string TEST_KEMPE_RANDOM_VS_ITER = "test_kempe_random_vs_iter";
    static const std::string TEST_DELTA = "test_delta";
    static const std::string TEST_DELTA_REMOVE_ADD = "test_delta_remove_add";
    static const std::string INFO = "info";
    static const std::string BRUTE = "brute";
    static const std::string TEST_ITERATION_YIELD_VS_STATE = "test_iteration_yield_vs_state";
    static const std::string TEST_KEMPE_LARGE_COUNT = "test_kempe_large_count";
    static const std::string ITERATED_GREEDY = "iterated-greedy";
    static const std::string MY_ITERATED_GREEDY = "my-iterated-greedy";
    static const std::string MY_ITERATED_GREEDY_INIT = "my-iterated-greedy-init";

    static const std::string IDENTITY = "identity";

    static const std::string MY_ILS = "my-ils";

    static const std::vector<std::string> available = {
        ILS, MY_ILS, TABU, FIRST, BEST, SA_BSU, VND, ITERATED_GREEDY, MY_ITERATED_GREEDY,

        // below "test" or "info" => NoSearch
        BRUTE, INFO, INTERACTIVE, TEST_INIT, TEST_KEMPE, TEST_DELTA, TEST_DELTA_REMOVE_ADD,
        TEST_KEMPE_RANDOM_VS_ITER, TEST_ITERATION_YIELD_VS_STATE, TEST_KEMPE_LARGE_COUNT
    };

    auto rep = [](std::string s) -> string {
        for(char& c : s)
            if(c == '_')
                c = '-';
        return s;
    };

    auto is = [rep](std::string const& a, std::string const& b) {
        return rep(a) == rep(b);
    };

    prs::TabLevel level;

    if(checkTokenParams(tm, ILS, {"search", "term", "perturb", "acc"})) {
        auto a = search(tm);
        auto b = term(tm);
        auto c = per(tm);
        auto d = acc(tm);
        return new emili::IteratedLocalSearch(*a,*b,*c,*d);
    }
    else if(checkTokenParams(tm, MY_ILS, {"search", "term", "perturb", "acc"})) {
        auto a = search(tm);
        auto b = term(tm);
        auto c = per(tm);
        auto d = acc(tm);
        return new emili::MyIteratedLocalSearch(a,b,c,d);
    }
    else if(tm.checkToken(TABU))
    {
        printTab(TABU);
        return tparams(tm);
    }
    else if(checkTokenParams(tm, IDENTITY, {"init"}))
    {
        auto a = init(tm);
        return new emili::IdentityLocalSearch(*a);
    }
    else if(checkTokenParams(tm, FIRST, {"init", "term", "neigh"}))
    {
        auto a = init(tm);
        auto b = term(tm);
        auto c = neigh(tm);
        return new emili::FirstImprovementSearch(*a, *b, *c);
    }
    else if(checkTokenParams(tm, BEST, {"init", "term", "neigh"}))
    {
        auto a = init(tm);
        auto b = term(tm);
        auto c = neigh(tm);
        return new emili::BestImprovementSearch(*a, *b, *c);
    }
    else if(checkTokenParams(tm, SA, {"init", "nei", "inittemp", "acceptance", "cooling", "term", "templ"}))
    {
        auto initsol = init(tm);
        auto nei = neigh(tm);
        return sa.buildSA(tm, initsol, nei);
    }
    else if(is(tm.peek(), SA_BSU)){
        printTab(SA_BSU);
        tm.next();

        ArgParser p;
        p.addInt("factor", 1);
        p.addInt("factor-exponent", 0);
        p.addInt("freq", 0);
        p.addBool("percent", true);
        p.addInt("percent-kempe", 0);

        p.parse(tm);
        p.print();

        if(p.has("factor") && p.has("factor-exponent"))
            throw std::invalid_argument("factor and factor-exponent given");

        int factor =
            p.has("factor") ? p.Int("factor") :
            p.has("factor-exponent") ? std::pow(10, p.Int("factor-exponent")) : 1;

        int freq = p.Int("freq");
        bool percent = p.Bool("percent");

        int baselineitmax = 500000000; // 5e8
        int itmax = baselineitmax / factor;

        // baseline
        double t0 = 700.00;
        double tr = 1390.62;
        double alpha = 0.99;
        double rho = 0.14;
        double sr = 0.70;
        double wH = 21;

        int Cap = 0;
        for(emili::ExamTT::Room& r : instance.rooms)
            Cap += r.capacity;
        Cap *= instance.P();

        int PHC = instance.examsCoincidence.size() + instance.examsAfter.size() + instance.examsExclusion.size();
        int RHC = instance.examsRoomsExclusive.size();
        int FLP = instance.institutionalWeightings.frontload.time;

        int P = instance.P(),
            R = instance.R(),
            E = instance.E(),
            S = instance.students.size();

        double
            CD = instance.meta.ConflictDensity,
            PC = 0, // ?
            ExR = 1.0 * E / (P * R),
            SxE = [this, E, S](){ // sum(map(instance.numberOfExams, instance.students)) / E
                /*
                int s = 0;
                for(int student : instance.students)
                    s += instance.numberOfExamsOfStudent(student);
                return 1.0 * s / E;
                */

                std::map<int, int> numberOfExamsOfStudent;

                for(int student : instance.students)
                    numberOfExamsOfStudent[student] = 0;

                for(auto& e : instance.exams)
                    for(int student : e.students)
                        numberOfExamsOfStudent[student]++;

                int s = 0;
                for(auto p : numberOfExamsOfStudent)
                    s += p.second;

                return 1.0 * s / E;
            }(),
            SCap = 1.0 * S / Cap; // StudentCapacityRatio ?


        // *) Check validator
        // *) SxE ExR and PC
        // 1) run more isntances
        // 2) final temperature
        // 3) how many accepting moves / temperature
        // 4) cooling scheme
        //

        instance.hardWeight = 16.73100
                + 102.30000 * CD
                - 0.48330 * P
                - 0.17740 * SxE
                - 1.23900 * ExR
                - 0.00076 * S
                + 0.11590 * PHC
                + 0.66660 * FLP
                + 0.78080 * RHC
                + 32.46000 * SCap,
                + 0.10100 * PC;

        // derived
        // double tmin = 0.50;
        constexpr bool wantToUseIrace = false;
        constexpr bool useIrace = itmax == baselineitmax && wantToUseIrace;

        int nS = useIrace ? 764141 : (int)(- itmax / (std::log(tr) / std::log(alpha)));
        int nA = useIrace ? 106980 : rho * nS;

        cout << "Features as in BSU Table 2 : " << endl
             << "E S P R PHC RHC FLP CD ExR SxE S/Cap PC" << endl
             << E << ' ' << S << ' ' << P << ' ' << R << ' '<< PHC << ' '<< RHC << ' '
             << FLP << ' '<< CD << ' ' << ExR << ' ' << SxE << ' ' << SCap << ' ' << PC << ' '
             << endl
             << "nA " << nA << " " << "nS " << nS
             << endl;

        cout << setw(20) << "Hard weight " << instance.hardWeight << endl;

        // 251_629
        // 111_107

        // 206_134: assign = [[33, 6], [50, 3], [15, 3], [51, 3], [13, 0], [11, 2], [21, 5], [38, 4], [14, 2], [25, 6], [53, 3], [9, 2], [45, 3], [36, 2], [18, 2], [21, 6], [47, 0], [42, 1], [23, 1], [9, 0], [21, 5], [53, 2], [28, 2], [41, 0], [17, 3], [8, 3], [45, 3], [32, 3], [38, 6], [29, 3], [3, 2], [45, 0], [27, 2], [49, 3], [17, 3], [33, 2], [12, 0], [18, 6], [6, 4], [52, 5], [41, 4], [51, 6], [48, 6], [39, 3], [6, 2], [51, 1], [5, 5], [19, 0], [39, 6], [6, 4], [11, 5], [50, 2], [11, 1], [27, 0], [11, 6], [53, 0], [29, 3], [19, 0], [42, 0], [5, 6], [37, 4], [3, 5], [3, 0], [18, 6], [19, 4], [22, 6], [5, 3], [8, 6], [12, 4], [2, 1], [1, 3], [45, 5], [49, 3], [33, 1], [29, 5], [49, 4], [3, 4], [1, 4], [49, 2], [4, 0], [7, 6], [6, 2], [23, 1], [35, 2], [14, 5], [20, 2], [4, 0], [51, 0], [44, 2], [40, 4], [6, 1], [3, 0], [34, 2], [21, 3], [50, 4], [6, 2], [18, 1], [43, 6], [5, 5], [0, 2], [17, 6], [6, 4], [51, 5], [52, 3], [50, 3], [25, 6], [13, 4], [51, 2], [28, 4], [29, 5], [7, 1], [6, 3], [19, 3], [11, 3], [9, 4], [41, 5], [42, 0], [32, 5], [23, 5], [23, 4], [6, 3], [18, 3], [39, 6], [21, 2], [8, 4], [34, 4], [12, 4], [13, 5], [21, 6], [40, 1], [42, 3], [47, 0], [18, 6], [27, 2], [17, 1], [38, 1], [7, 1], [46, 3], [25, 2], [17, 0], [49, 4], [50, 1], [8, 5], [3, 2], [43, 5], [31, 6], [16, 2], [13, 5], [52, 0], [8, 0], [46, 1], [12, 1], [4, 4], [4, 2], [33, 5], [39, 3], [39, 3], [51, 2], [28, 3], [28, 3], [14, 4], [45, 0], [1, 6], [33, 4], [5, 2], [9, 0], [40, 0], [6, 1], [44, 4], [48, 3], [9, 4], [21, 2], [1, 5], [26, 6], [27, 5], [12, 4], [23, 1], [28, 2], [28, 2], [15, 5], [37, 3], [26, 2], [13, 1], [14, 6], [1, 3], [38, 6], [5, 4], [10, 5], [20, 0], [8, 3], [49, 1], [17, 2], [10, 1], [2, 1], [12, 1], [32, 3], [48, 4], [22, 0], [41, 0], [24, 2], [50, 4], [2, 4], [19, 1], [35, 2], [2, 1], [13, 5], [24, 4], [6, 6], [48, 0], [51, 4], [0, 3], [48, 3], [43, 6], [37, 5], [26, 4], [46, 5], [13, 6], [41, 1], [31, 3], [14, 0], [49, 6], [41, 2], [53, 6], [23, 1], [13, 3], [34, 2], [10, 2], [10, 3], [31, 2], [34, 3], [16, 1], [23, 1], [37, 4], [23, 3], [32, 3], [3, 5], [14, 4], [40, 3], [38, 1], [24, 0], [17, 2], [43, 3], [33, 2], [16, 0], [38, 0], [14, 6], [1, 4], [36, 2], [13, 4], [42, 6], [47, 0], [47, 6], [2, 5], [8, 5], [25, 3], [11, 4], [35, 3], [2, 2], [7, 0], [27, 5], [1, 6], [16, 4], [17, 1], [31, 6], [4, 6], [34, 0], [43, 0], [24, 4], [30, 4], [2, 6], [27, 0], [36, 6], [11, 5], [44, 1], [39, 5], [12, 1], [25, 1], [4, 5], [41, 3], [7, 0], [47, 3], [52, 5], [42, 2], [3, 2], [49, 5], [28, 4], [0, 6], [24, 1], [24, 6], [7, 0], [50, 4], [7, 2], [20, 4], [43, 2], [33, 3], [2, 3], [20, 5], [44, 4], [4, 2], [14, 0], [38, 5], [1, 3], [34, 0], [9, 3], [44, 2], [32, 2], [44, 4], [52, 4], [46, 3], [12, 6], [9, 6], [31, 3], [41, 4], [43, 6], [24, 6], [36, 3], [25, 6], [12, 0], [32, 5], [17, 5], [47, 4], [14, 2], [53, 6], [33, 2], [38, 5], [14, 6], [52, 2], [47, 3], [35, 1], [45, 3], [24, 5], [9, 3], [1, 6], [8, 2], [33, 6], [20, 1], [15, 6], [7, 1], [18, 6], [35, 1], [15, 5], [5, 6], [30, 2], [3, 2], [7, 1], [7, 6], [53, 0], [12, 2], [3, 5], [11, 4], [21, 1], [31, 6], [9, 4], [34, 1], [50, 0], [32, 6], [10, 4], [5, 0], [27, 2], [45, 3], [18, 2], [39, 0], [50, 1], [2, 4], [37, 2], [4, 1], [14, 1], [52, 4], [50, 4], [32, 0], [39, 5], [17, 0], [20, 1], [30, 2], [12, 6], [11, 2], [17, 0], [47, 6], [12, 4], [42, 0], [53, 2], [4, 5], [3, 6], [11, 2], [44, 6], [31, 5], [38, 5], [41, 1], [26, 4], [53, 5], [21, 0], [6, 5], [35, 3], [0, 6], [19, 1], [45, 4], [9, 1], [3, 3], [41, 2], [51, 5], [30, 6], [6, 0], [46, 4], [32, 6], [8, 2], [27, 1], [52, 0], [51, 5], [38, 1], [5, 4], [31, 3], [41, 2], [32, 4], [49, 4], [43, 6], [45, 4], [42, 2], [40, 2], [1, 2], [50, 1], [53, 3], [19, 2], [53, 6], [35, 6], [10, 5], [8, 6], [52, 0], [19, 6], [26, 0], [40, 6], [9, 6], [2, 5], [46, 5], [32, 1], [15, 2], [0, 6], [31, 0], [23, 6], [28, 1], [33, 3], [9, 4], [8, 0], [20, 1], [0, 1], [44, 1], [1, 6], [6, 2], [17, 5], [28, 4], [5, 3], [29, 3], [46, 5], [9, 3], [14, 6], [33, 2], [42, 5], [52, 4], [8, 5], [27, 2], [16, 3], [8, 1], [31, 0], [3, 4], [8, 4], [37, 3], [33, 4], [32, 0], [6, 5], [2, 0], [4, 6], [21, 4], [27, 4], [9, 0], [8, 5], [22, 6], [52, 4], [52, 0], [49, 6], [13, 5], [39, 0], [44, 4], [25, 0], [15, 6], [51, 2], [24, 6], [16, 5], [20, 1], [29, 4], [27, 6], [8, 1], [38, 2], [30, 1], [35, 3], [45, 3], [14, 6], [50, 2], [5, 1], [8, 4], [5, 3], [28, 4], [20, 2], [42, 5], [25, 1], [19, 5], [31, 2], [14, 0], [36, 1], [49, 0], [19, 5], [29, 5], [17, 0], [15, 4], [43, 4], [4, 6], [43, 6], [0, 1], [5, 4], [12, 4], [23, 5], [43, 6], [51, 0], [41, 3], [30, 6], [3, 6], [50, 5], [38, 5], [46, 3], [6, 2], [11, 0], [48, 4], [24, 1], [46, 0], [16, 4], [48, 0], [53, 1], [3, 0], [34, 1], [18, 4], [14, 2], [51, 4], [34, 0], [43, 5], [5, 4], [37, 6], [26, 3], [23, 1], [2, 6], [7, 5], [8, 3], [39, 3], [20, 2], [24, 6], [45, 6], [4, 0], [10, 4], [4, 5], [0, 4], [6, 4], [46, 3], [45, 4], [31, 5], [24, 6], [46, 6], [27, 0], [50, 2], [1, 5], [38, 5], [53, 1], [2, 0], [26, 2], [20, 1], [6, 2], [48, 0], [39, 4], [21, 3], [9, 3], [2, 0], [42, 4], [8, 0], [32, 4], [29, 6], [28, 4], [45, 1], [20, 5], [18, 1], [49, 0], [11, 1], [25, 3], [22, 4], [9, 0], [23, 6], [27, 6], [21, 3], [49, 6], [7, 4], [33, 0], [24, 0], [27, 0], [33, 0], [43, 5], [32, 6], [20, 6], [50, 0], [14, 2], [44, 2], [52, 4], [11, 3]]; hard = (455, 7843, 0, 6, 2, 0, 0); soft = (3059, 0, 9062, 9062, 8700, 265, 1560);

        auto init = new emili::ExamTT::RandomInitialSolution(instance);

        emili::Neighborhood* neigh = nullptr;
        if(p.Int("percent-kempe") == 0) {
            neigh = new emili::ExamTT::MixedMoveSwapNeighborhood(instance, sr);
        } else {
            double k = p.Int("percent-kempe") / 100.0;
            neigh = new emili::ExamTT::MixedRandomNeighborhoodProba(
                {
                    new emili::ExamTT::KempeChainNeighborhood(instance),
                    new emili::ExamTT::SwapNeighborhood(instance),
                    new emili::ExamTT::MoveNeighborhood(instance),
                },
                {
                    (float) k,
                    (float) (sr * (1 - k))
                } // (1 - sr) * (1 - k)
            );
        }

        auto acc = new SAMetropolisAcceptanceDebug(t0);

        auto term = percent ?
                    (SAMaxIterTermination*) new SAMaxIterTerminationDebug(itmax) :
                    (SAMaxIterTermination*) new SAMaxIterTermination(itmax);

        auto initialTemperature = new FixedInitTemp();
        initialTemperature->set(t0);

        auto tempLength = new CappedMaxAcceptedTempLength(nA, nS);
        auto tempRestart = new SANoRestart();

        auto cooling = new LinearCooling(alpha, initialTemperature);
        cooling->setTempLength(tempLength); // Shouldnt SimulatedAnnealing class do it ?
        cooling->setTempRestart(tempRestart); // Shouldnt SimulatedAnnealing class do it ?

        auto explo = new SARandomExplorationNoCopyDebug(neigh, acc, term);

        auto sa = new SimulatedAnnealing(init, initialTemperature, acc, cooling, tempRestart, term, tempLength, explo, neigh);

        if(percent) {
            ((SAMaxIterTerminationDebug*)term)->setOnEachPercent(freq, [sa,explo,acc]{
                cout << std::setprecision(10) << std::fixed <<
                    "&bestcost=" << sa->status->best_cost <<
                    "&curcost=" << sa->getBestSoFar()->getSolutionValue() <<
                    "&temp=" << sa->status->temp <<
                    "&nacc=" << explo->nacc <<
                    "&nnacc=" << explo->nnacc <<
                    "&tot=" << sa->status->total_counter <<
                    "&sumprob=" << acc->sumprob <<
                    "&nnotimproving=" << acc->nnotimproving <<
                endl;

                explo->nacc = 0;
                explo->nnacc = 0;
                acc->nnotimproving = 0;
                acc->sumprob = 0;
            });
        }

        return sa;
    }
    else if(tm.checkToken(VND))
    {
        return vparams(tm);
    }
    else if(checkTokenParams(tm, ITERATED_GREEDY, {"cons", "termin", "destr", "accept"}))
    {
        auto c = constructor(tm);
        auto t = term(tm);
        auto d = destructor(tm);
        auto ac = acc(tm);
        return new emili::IteratedGreedy(*c,*t,*d,*ac);
    }
    else if(checkTokenParams(tm, MY_ITERATED_GREEDY, {"search", "cons", "termin", "destr", "accept"}))
    {
        // use constructFull as init
        auto ls = search(tm);
        auto c = constructor(tm);
        auto t = term(tm);
        auto d = destructor(tm);
        auto ac = acc(tm);
        return new emili::MyIteratedGreedy(c,t,d,ac,ls);
    }
    else if(checkTokenParams(tm, MY_ITERATED_GREEDY_INIT, {"init", "search", "cons", "termin", "destr", "accept"})) {
        auto a = init(tm);
        auto b = search(tm);
        auto c = constructor(tm);
        auto t = term(tm);
        auto d = destructor(tm);
        auto e = acc(tm);
        return new emili::MyIteratedGreedyInit(a,c,t,d,e,b);
    }
    else if(tm.checkToken(INFO)) {
        instance.presentation(std::cout);

        throw NoSearch();
    }
    else if(tm.peek() == TEST_DELTA || tm.peek() == TEST_DELTA_REMOVE_ADD) {
        bool isTestDelta = tm.peek() == TEST_DELTA;
        bool isTestDeltaRemoveAdd = tm.peek() == TEST_DELTA_REMOVE_ADD;
        tm.next();

        ArgParser p;

        p.addInt("N", 1000);
        p.addBool("each-move", true);

        if(isTestDeltaRemoveAdd) {
            p.addInt("G", 1);
            p.addBool("each-move-partial", false);
        }

        p.parse(tm);
        p.print();

        int N = p.Int("N");
        int G = p.Int("G");
        bool checkEachMove = p.Bool("each-move");
        bool checkEachMovePartial = p.Bool("each-move-partial");

        if(isTestDelta)
            emili::ExamTT::test::delta(instance, N, checkEachMove);
        else if(isTestDeltaRemoveAdd)
            emili::ExamTT::test::deltaRemoveAdd(instance, N, checkEachMove, G, checkEachMovePartial);

        throw NoSearch();
    }
    else if(tm.checkToken(INTERACTIVE)) {
        emili::ExamTT::test::interactive(instance);
        throw NoSearch();
    }
    else if(tm.checkToken(TEST_KEMPE_LARGE_COUNT)) {
        emili::ExamTT::stats::kempe_compare_size_fast_iter(instance);
        throw NoSearch();
    }
    else if(tm.peek() == TEST_KEMPE || tm.peek() == TEST_KEMPE_RANDOM_VS_ITER) {
        bool isTestKempe = tm.peek() == TEST_KEMPE;
        bool isTestKempeRandomVsIter = tm.peek() == TEST_KEMPE_RANDOM_VS_ITER;
        tm.next();

        std::vector<int> x = {-1};
        while(tm.checkInteger(x.back()))
            x.push_back(-1);
        x.pop_back();

        if(x.empty())
            genericError("ints... are expected");

        const int P = instance.P();
        for(int y : x)
            if(!(0 <= y && y < P))
                genericError(cerr << "int " << y << " is not a valid period in [0," << P << "[");

        ArgParser p;
        p.addBool("use-fast-iter", false);
        p.addBool("use-iterate", false);

        if(isTestKempeRandomVsIter)
            p.addInt("N", 1000);

        p.parse(tm);
        p.print();

        if(isTestKempe) {
            emili::ExamTT::stats::kempe_print_iteration(instance, x, p.Bool("use-fast-iter"), p.Bool("use-iterate"));
        } else if(isTestKempeRandomVsIter) {
            emili::ExamTT::test::kempe_iteration_vs_random(instance, x, p.Int("N"), p.Bool("use-fast-iter"), p.Bool("use-iterate"));
        }

        throw NoSearch();
    }
    else if(tm.checkToken(TEST_ITERATION_YIELD_VS_STATE)) {
        ArgParser p;
        p.addBool("use-fast-iter", false);

        p.parse(tm);
        p.print();

        emili::ExamTT::ExamTTSolution sol;
        sol.initRandom(instance);

        emili::Neighborhood* n = p.Bool("use-fast-iter")
            ? new emili::ExamTT::KempeChainNeighborhoodFastIter(instance)
            : new emili::ExamTT::KempeChainNeighborhood(instance);

        try {
            emili::ExamTT::test::iterateVsComputeStep(&sol, n);
        } catch(std::exception& e) {
            delete n;
            throw e;
        }

        delete n;

        throw NoSearch();
    }
    else if(tm.checkToken(TEST_INIT))
    {
        emili::InitialSolution* ini = init(tm);
        emili::Solution* s = ini->generateSolution();
        std::cout << s->getSolutionRepresentation() << std::endl;
        std::cout << s->getSolutionValue() << std::endl;
        std::cerr << s->getSolutionValue() << std::endl;

        throw NoSearch();
    }
    else if(tm.checkToken(BRUTE)) {
        int n = -1;
        if(tm.checkToken("n")) {
            n = tm.getInteger();
        }

        emili::ExamTT::BruteForce* ret = nullptr;
        if(tm.checkToken("[")) {
            std::vector<int> ints;
            int x;
            while(tm.checkInteger(x))
                ints.push_back(x);

            if(ints.size() % 2 == 0) {
                if(tm.checkToken("]")) {
                    std::vector<std::pair<int,int>> assignements;
                    for(int i = 0; i < ints.size(); i += 2)
                        assignements.push_back({ints[i], ints[i+1]});
                    ret = new emili::ExamTT::BruteForce(instance, assignements);
                } else {
                    errorExpected(tm, "]", {"]"});
                }
            } else {
                cerr << "(BRUTE) Even number of integers expected, got " << ints.size() << endl;
                exit(-1);
            }

        } else {
            ret = new emili::ExamTT::BruteForce(instance);
        }

        if(n >= 0)
            ret->setSizeDebug(n);

        return ret;
    }

    errorExpected(tm, "SEARCH", available);

    return nullptr;
}

SimulatedAnnealing* ExamTTParser::SAParser::buildSA(prs::TokenManager& tm, emili::InitialSolution* initsol, emili::Neighborhood* nei) {
    auto inittemp = INITTEMP(tm, initsol, nei);
    auto acceptance = ACCEPTANCE(tm, inittemp);
    auto cooling = COOL(tm, inittemp, nei);
    auto temprestart = TEMPRESTART(tm, inittemp, nei);
    cooling->setTempRestart(temprestart);
    auto term = TERMINATION(tm, nei);
    auto templ = TEMPLENGTH(tm, nei);
    cooling->setTempLength(templ);
    auto explo = EXPLORATION(tm, nei, acceptance, term);

    return new SimulatedAnnealing(initsol, inittemp, acceptance, cooling, temprestart, term, templ, explo, nei);
}

emili::Perturbation* ExamTTParser::per(prs::TokenManager& tm)
{
    printTab(*tm);

    static const std::string NO = "no";
    static const std::string NOPERTURB = "nopertub";
    static const std::string RandomMovePerturbationInPlace = "random-move";
    static const std::string VNRandomMovePerturbationInPlace = "variable-neigh-random-move";
    static const std::string constructperturb = "construct-perturb";

    static const std::vector<std::string> available = {
        NO, RandomMovePerturbationInPlace,
    };

    prs::TabLevel level;

    if(tm.checkToken(NOPERTURB) || tm.checkToken(NO)) {
        return new emili::NoPerturbation();
    } else if(checkTokenParams(tm, RandomMovePerturbationInPlace, {"rneigh", "steps"})) {
        auto n = neigh(tm);
        int s = getIntParam(tm, "s");
        return new emili::RandomMovePerturbationInPlace(*n, s);
    } else if(checkTokenParams(tm, VNRandomMovePerturbationInPlace, {"steps", "iter", "rneigh..."})) {
        TabLevel l;
        int steps = getIntParam(tm, "steps");
        int iterations = getIntParam(tm, "iter");
        auto neigh = neighs(tm);
        return new emili::VNRandomMovePerturbationInPlace(neigh, steps, iterations);
    } else if(checkTokenParams(tm, constructperturb, {"construct", "destruct"})) {
        auto a = constructor(tm);
        auto b = destructor(tm);
        return new emili::ConstructDestructPertub(a,b);
    }

    errorExpected(tm, "PERTURBATION", available);
    return nullptr;
}

emili::Acceptance* ExamTTParser::acc(prs::TokenManager& tm)
{
    prs::TabLevel level;

    if(tm.checkToken(ACCEPTANCE_METRO))
    {
        float n = tm.getDecimal();
        printTab("metropolis acceptance. temperature : " + to_string(n));
        return new emili::MetropolisAcceptance(n);
    }
    else if(tm.checkToken(ACCEPTANCE_ALWAYS))
    {
        emili::accept_candidates accc;
        string t1 = *tm;

        if(tm.checkToken(ACCEPTANCE_INTENSIFY))
            accc = emili::ACC_INTENSIFICATION;

        else if(tm.checkToken(ACCEPTANCE_DIVERSIFY))
            accc = emili::ACC_DIVERSIFICATION;

        else
            errorExpected(tm, ACCEPTANCE_INTENSIFY, {ACCEPTANCE_INTENSIFY, ACCEPTANCE_DIVERSIFY});

        printTab("Acceptance always " + t1);
        return new emili::AlwaysAccept(accc);
    }
    else if(tm.checkToken(ACCEPTANCE_IMPROVE)) {

        printTab("Improve acceptance");

        return new  emili::ImproveAccept();
    }
    else if(tm.checkToken(ACCEPTANCE_SA_METRO))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        printTab("Metropolis acceptance. start, end, ratio : " + to_string(start) + ", " + to_string(end)+ "," + to_string(ratio));
        return new emili::Metropolis(start,end,ratio);
    }
    else if(tm.checkToken(ACCEPTANCE_PMETRO))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        int iterations = tm.getInteger();
        ostringstream oss; oss << "metropolis acceptance. start, end, ratio, frequence : "<< start << ", "<< end << "," << ratio <<"," << iterations;
        printTab(oss.str());
        return new emili::Metropolis(start,end,ratio,iterations);
    }
    else if(tm.checkToken(ACCEPTANCE_SA))
    {
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        int iterations = tm.getInteger();
        float alpha =tm.getDecimal();
        ostringstream oss; oss << "metropolis acceptance. start ,end , ratio, frequence, alpha : "<< start << ", "<< end << "," << ratio <<","<< iterations << "," << alpha;
        printTab(oss.str());
        return new emili::Metropolis(start,end,ratio,iterations,alpha);
    }
    else if(tm.checkToken(ACCEPTANCE_IMPROVE_PLATEAU))
    {
        int plateau_steps = tm.getInteger();
        int threshold = tm.getInteger();
        ostringstream oss; oss << "Accept a diversification solution if it improves on the intensification otherwise it will accept "<< plateau_steps << " non improving steps once it reaches the threshold of " << threshold;
        printTab(oss.str());
        return new emili::AcceptPlateau(plateau_steps,threshold);
    }

    errorExpected(tm, "ACCEPTANCE_CRITERIA", {ACCEPTANCE_METRO, ACCEPTANCE_ALWAYS, ACCEPTANCE_INTENSIFY, ACCEPTANCE_DIVERSIFY, ACCEPTANCE_IMPROVE, ACCEPTANCE_SA_METRO, ACCEPTANCE_PMETRO, ACCEPTANCE_SA, ACCEPTANCE_IMPROVE_PLATEAU});
    return nullptr;
}

emili::BestTabuSearch* ExamTTParser::tparams(prs::TokenManager& tm)
{
    if(tm.checkToken(BEST))
    {
        auto a = init(tm);
        auto b = term(tm);
        auto c = neigh(tm);
        auto tmem = tmemory(c, tm);
        return new emili::BestTabuSearch(*a, *b, *c, *tmem);
    }
    else if(tm.checkToken(FIRST))
    {
        auto a = init(tm);
        auto b = term(tm);
        auto c = neigh(tm);
        auto tmem = tmemory(c, tm);
        return new emili::FirstTabuSearch(*a, *b, *c, *tmem);
    }

    errorExpected(tm, "PIVOTAL_RULE", {BEST,FIRST});
    return nullptr;
}

emili::TabuMemory* ExamTTParser::tmemory(emili::Neighborhood* n,prs::TokenManager& tm)
{
    prs::TabLevel level;

    emili::TabuMemory* tmem = nullptr;

    /*TODO HERE GOES THE CODE TO INSTATIATE A TABUTENURE IMPLEMENTATION*/

    return tmem;
}

emili::LocalSearch* ExamTTParser::vparams(prs::TokenManager& tm)
{
    printTab("VND SEARCH");

    prs::TabLevel level;

    if(checkTokenParams(tm, FIRST, {"init", "term", "neigh"}))
    {
        auto in = init(tm);
        auto te = term(tm);
        auto nes = neighs(tm);
        return new emili::VNDSearch<emili::FirstImprovementSearch>(*in,*te,nes);
    }
    else if(checkTokenParams(tm, BEST, {"init", "term", "neigh"}))
    {
        auto in = init(tm);
        auto te = term(tm);
        auto nes = neighs(tm);
        return new emili::VNDSearch<emili::BestImprovementSearch>(*in,*te,nes);
    }

    errorExpected(tm, "VND_SEARCH", {FIRST, BEST});
    return nullptr;
}

emili::ExamTT::InsertHeuristic* ExamTTParser::insertHeuristic(prs::TokenManager& tm) {
    static const std::string BestInsertHeuristic = "best";

    TabLevel l;

    const std::vector<std::string> available = {
        BestInsertHeuristic
    };

    if(checkTokenParams(tm, BestInsertHeuristic, {})) {
        return new emili::ExamTT::BestInsertHeuristic(instance);
    }

    errorExpected(tm, "INSERT_HEURISTIC", available);
    return nullptr;
}

emili::Constructor* ExamTTParser::constructor(prs::TokenManager& tm) {
    static const std::string RandomOrderInserter = "random-order";
    static const std::string DegreeInserter = "degree";

    const std::vector<std::string> available = {
        RandomOrderInserter, DegreeInserter,
    };

    TabLevel l;

    if(checkTokenParams(tm, RandomOrderInserter, {"heur"})) {
        auto heur = insertHeuristic(tm);
        return new emili::ExamTT::RandomOrderInserter(instance, heur);
    } else if(checkTokenParams(tm, DegreeInserter, {"heur"})) {
        auto heur = insertHeuristic(tm);
        return new emili::ExamTT::DegreeInserter(instance, heur);
    }

    errorExpected(tm, "CONSTRUCTOR", available);
    return nullptr;
}

emili::Destructor* ExamTTParser::destructor(prs::TokenManager& tm) {
    static const std::string FixedRandomDestructor = "fixed-random";
    static const std::string NRandomDaysDestructor = "random-days";
    static const std::string NGroupedDaysDestructor = "grouped-days";
    static const std::string NGroupedPeriodsDestructor = "grouped-periods";

    TabLevel l;

    const std::vector<std::string> available = {
        FixedRandomDestructor, NRandomDaysDestructor, NGroupedDaysDestructor,
        NGroupedPeriodsDestructor,
    };

    if(checkTokenParams(tm, FixedRandomDestructor, {"N"})) {
        TabLevel l;
        int N = getIntPercentParam(tm, "N", instance.E());
        return new emili::ExamTT::FixedRandomDestructor(instance, N);
    } else if(checkTokenParams(tm, NRandomDaysDestructor, {"N"})) {
        TabLevel l;
        int N = getIntPercentParam(tm, "N", instance.Days());
        return new emili::ExamTT::NRandomDaysDestructor(instance, N);
    } else if(checkTokenParams(tm, NGroupedDaysDestructor, {"N"})) {
        TabLevel l;
        int N = getIntPercentParam(tm, "N", instance.Days());
        return new emili::ExamTT::NGroupedDaysDestructor(instance, N);
    } else if(checkTokenParams(tm, NGroupedPeriodsDestructor, {"N"})) {
        TabLevel l;
        int N = getIntPercentParam(tm, "N", instance.P());
        return new emili::ExamTT::NGroupedPeriodsDestructor(instance, N);
    }

    errorExpected(tm, "DESTRUCTOR", available);
    return nullptr;
}

emili::InitialSolution* ExamTTParser::init(prs::TokenManager& tm)
{
    prs::TabLevel level;
    static const std::string INITIAL_CONSTRUCTOR = "initial-constructor";

    static const std::vector<std::string> available = {
        INITIAL_RANDOM, INITIAL_CONSTRUCTOR
    };

    printTab(*tm);

    if(tm.checkToken(INITIAL_RANDOM)) {
        return new emili::ExamTT::RandomInitialSolution(instance);
    } else if(tm.checkToken(INITIAL_CONSTRUCTOR)){
        auto cons = constructor(tm);
        return new emili::ExamTT::ConstructorInitialSolution(instance, cons);
    }

    errorExpected(tm, "INITIAL_SOLUTION", available);
    return nullptr;
}

emili::Termination* ExamTTParser::term(prs::TokenManager& tm)
{
    prs::TabLevel level;

    static const std::vector<std::string> available = {
        TERMINATION_LOCMIN, TERMINATION_WTRUE, TERMINATION_TIME, TERMINATION_MAXSTEPS
    };

    if(tm.checkToken(TERMINATION_LOCMIN))
    {
        printTab(TERMINATION_LOCMIN);
        return new emili::LocalMinimaTermination();
    }
    else if(tm.checkToken(TERMINATION_WTRUE))
    {
        printTab(TERMINATION_WTRUE);
        return new emili::WhileTrueTermination();
    }
    else if(tm.checkToken(TERMINATION_TIME))
    {
        float time = tm.getDecimal();
        if(time == 0)
            time = 1;

        printTab("Timed termination. ratio: " + to_string(time));
        return new emili::TimedTermination(time);
    }
    else if(tm.checkToken(TERMINATION_MAXSTEPS))
    {
        int steps = tm.getInteger();
        printTab("Max Steps termination. # steps: " + to_string(steps));
        return new emili::MaxStepsTermination(steps);
    }

    errorExpected(tm, "TERMINATION_CRITERIA", available);
    return nullptr;
}

emili::Neighborhood* ExamTTParser::neigh(prs::TokenManager& tm, bool errorIfNotFound)
{
    const static std::string MoveNeighborhood = "move";
    const static std::string SwapNeighborhood = "swap";
    const static std::string KempeChainNeighborhoodFastIter = "kempe";

    const static std::vector<std::string> available = {
        MoveNeighborhood, SwapNeighborhood, KempeChainNeighborhoodFastIter
    };

    prs::TabLevel level;

    printTab(*tm);

    if(tm.checkToken(MoveNeighborhood)) {
        return new emili::ExamTT::MoveNeighborhood(instance);
    }
    else if(tm.checkToken(SwapNeighborhood)) {
        return new emili::ExamTT::SwapNeighborhood(instance);
    }
    else if(tm.checkToken(KempeChainNeighborhoodFastIter)) {
        return new emili::ExamTT::KempeChainNeighborhoodFastIter(instance);
    }

    if(errorIfNotFound)
        errorExpected(tm, "NEIGHBORHOOD", available);

    return nullptr;
}

std::vector<emili::Neighborhood*> ExamTTParser::neighs(prs::TokenManager& tm)
{
    std::vector<emili::Neighborhood*> vnds;

    if(tm.checkToken("[")) {
        while(! tm.checkToken("]"))
            vnds.push_back(neigh(tm, true));
    } else {
        vnds.push_back(neigh(tm, true));

        while(vnds.back() != nullptr)
            vnds.push_back(neigh(tm, false));

        vnds.pop_back();
    }

    if(vnds.empty())
        genericError("neighs+");

    return vnds;
}

void ExamTTParser::problem(prs::TokenManager& tm)
{
    tm.nextToken();
    emili::ExamTT::InstanceParser parser(tm.tokenAt(1));
    parser.parse(instance);
}

emili::LocalSearch* ExamTTParser::buildAlgo(prs::TokenManager& tm)
{
    problem(tm);
    emili::LocalSearch* local = search(tm);
    std::cout << "------" << std::endl;
    return local;
}

bool ExamTTParser::isParsable(std::string& problem)
{
    return PROBLEMS_DEF.count(problem); // return problem in PROBLEMS_DEF
}

std::string ExamTTParser::availableProblems() const
{
    ostringstream oss;

    for(string problem : PROBLEMS_DEF)
        oss << problem << " ";

    return oss.str();
}

} // namespace ExamTT
} // namespace prs
