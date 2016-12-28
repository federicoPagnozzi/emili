#include "examttparser.h"
#include <ctime>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <set>
#include <string>
#include <iomanip>
#include <map>
#include <functional>

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
    TERMINATION_MAXSTEPS_DEBUG = "maxstep-debug",
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
    oss << "EMILI INSTANCE_FILE_PATH ExamTT SEARCH [rnds seed]";
    return oss.str();
}

namespace {
bool checkTokenParams(prs::TokenManager& tm, std::string val, std::vector<std::string> const& params, std::string prefix = "") {
    if(tm.checkToken(val)) {
        std::ostringstream oss;
        if(! prefix.empty())
            oss << prefix << ": ";
        oss << val;
        /*
        if(params.size())
            oss << " (" << params.size() << ") ";
        for(auto const& p : params)
            oss << p << " ";
        */
        printTab(oss.str());
        return true;
    } else {
        return false;
    }
}

int getIntParam(prs::TokenManager& tm, std::string val) {
    TabLevel l;
    int i = tm.getInteger();
    std::ostringstream oss;
    oss << val << ": " << i;
    printTab(oss.str());
    return i;
}

float getDecimalParam(prs::TokenManager& tm, std::string val) {
    float i = tm.getDecimal();
    std::ostringstream oss;
    oss << val << ": " << i;
    printTab(oss.str());
    return i;
}

int getIntPercentParam(prs::TokenManager& tm, int max) {
    if(tm.checkToken("percent") || tm.checkToken("%")) {
        return tm.getInteger() * max / 100;
    } else {
        return tm.getInteger();
    }
}

float getFloatPercentParam(prs::TokenManager& tm, float max) {
    if(tm.checkToken("percent") || tm.checkToken("%")) {
        return tm.getDecimal() * max / 100;
    } else {
        return tm.getDecimal();
    }
}

}

double ExamTTParser::getHardweightFromFeatures(bool USE_FORMULA_1) {
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

    std::string canonicalInstanceName = instanceFilename.rfind('/') == std::string::npos
                ? instanceFilename
                : instanceFilename.substr(instanceFilename.rfind('/') + 1); // set1.exam

    std::cout << "instance name : " << canonicalInstanceName << endl;

    double
        CD = instance.meta.ConflictDensity,
        ExR = 1.0 * E / (P * R),
        SxE = [this, E, S](){ // sum(map(instance.numberOfExamsOfStudent, instance.students)) / E
            // int s = 0;
            // for(int student : instance.students)
            //     s += instance.numberOfExamsOfStudent(student);
            // return 1.0 * s / E;

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
        SCap = 1.0 * S / Cap, // StudentCapacityRatio ?
        PC = 0; // ?

    static map<string, tuple<double,double>> InstanceDataMap = {
        {"set1.exam", make_tuple(0.75, 15.93)},
        {"set2.exam", make_tuple(0.23, 26.45)},
        {"set3.exam", make_tuple(0.33, 34.11)},
        {"set4.exam", make_tuple(0.86, 19.05)},
        {"set5.exam", make_tuple(0.34, 72.62)},
        {"set6.exam", make_tuple(0.56, 35.00)},
        {"set7.exam", make_tuple(0.22, 43.63)},
        {"set8.exam", make_tuple(0.43, 177.00)},
        {"set9.exam", make_tuple(0.60, 32.80)},
        {"set10.exam", make_tuple(0.13, 89.38)},
        {"set11.exam", make_tuple(0.48, 51.08)},
        {"set12.exam", make_tuple(0.20, 36.67)},
    };

    if(InstanceDataMap.count(canonicalInstanceName))
        tie(SCap, PC) = InstanceDataMap[canonicalInstanceName];

    // *) Check validator
    // *) SxE ExR and PC
    // 1) run more isntances
    // 2) final temperature
    // 3) how many accepting moves / temperature
    // 4) cooling scheme

    double hardWeight;
    if(USE_FORMULA_1) {
        hardWeight = 23.45040
            - 0.01924   * E
            + 29.10456  * SCap
            - 191.44258 * CD
            - 0.05139   * PC
            - 0.11727   * R;
    } else {
        hardWeight = 16.73100
            + 102.30000 * CD
            - 0.48330 * P
            - 0.17740 * SxE
            - 1.23900 * ExR
            - 0.00076 * S
            + 0.11590 * PHC
            + 0.66660 * FLP
            + 0.78080 * RHC
            + 32.46000 * SCap
            + 0.10100 * PC;
    }

    if(hardWeight < 1) // also negative ones...
        hardWeight = 1;

    printTab("Features as in BSU Table 2 : ");
    printTab("E S P R PHC RHC FLP CD ExR SxE S/Cap PC");

    ostringstream oss;
    oss << E << ' ' << S << ' ' << P << ' ' << R << ' '<< PHC << ' '<< RHC << ' '
        << FLP << ' '<< CD << ' ' << ExR << ' ' << SxE << ' ' << SCap << ' ' << PC;
    printTab(oss.str());
    oss.str("");

    return hardWeight;
}

emili::LocalSearch* ExamTTParser::search(prs::TokenManager& tm, bool mustHaveInit, std::string prefix)
{
    static const std::string
        SA_BSU = "sa-bsu",
        INTERACTIVE = "interactive",
        SHOW_PARAMS = "show-params",
        TEST_KEMPE = "test_kempe",
        TEST_KEMPE_RANDOM_VS_ITER = "test_kempe_random_vs_iter",
        TEST_DELTA = "test_delta",
        TEST_DELTA_REMOVE_ADD = "test_delta_remove_add",
        INFO = "info",
        INFO_INIT = "info-init",
        BRUTE = "brute",
        TEST_ITERATION_YIELD_VS_STATE = "test_iteration_yield_vs_state",
        TEST_KEMPE_LARGE_COUNT = "test_kempe_large_count",
        ITERATED_GREEDY = "iterated-greedy",
        MY_ITERATED_GREEDY = "my-iterated-greedy",
        INIT_AND_SEARCH = "init-search",
        NO_SEARCH = "no-search",

        IDENTITY = "identity",
        MY_ILS = "my-ils",

        __DUMMY__;

    static const std::vector<std::string> available = {
        ILS, MY_ILS, TABU, FIRST, BEST, SA_BSU, VND, ITERATED_GREEDY, MY_ITERATED_GREEDY, INIT_AND_SEARCH, IDENTITY,

        // below "test" or "info" => NoSearch
        NO_SEARCH, BRUTE, INFO, INFO_INIT, INTERACTIVE, TEST_INIT, TEST_KEMPE, TEST_DELTA, TEST_DELTA_REMOVE_ADD,
        TEST_KEMPE_RANDOM_VS_ITER, TEST_ITERATION_YIELD_VS_STATE, TEST_KEMPE_LARGE_COUNT
    };

    prs::TabLevel level;

    if(checkTokenParams(tm, ILS, {"search", "termination", "perturbation", "acceptance"})) {
        auto a = search(tm, mustHaveInit); // the init will be the init of the search
        auto b = termination(tm);
        auto c = perturbation(tm);
        auto d = acceptance(tm);
        return new emili::IteratedLocalSearch(*a,*b,*c,*d);
    }
    else if(checkTokenParams(tm, MY_ILS, {"search", "term", "perturbation", "acc"})) {
        emili::LocalSearch* a = nullptr;
        emili::Termination* b = nullptr;
        emili::Perturbation* c = nullptr;
        emili::Acceptance* d = nullptr;

        // what about the hard weight ?
        if(tm.peekIs("{")) {
            ArgParser p;

            p.addFuncRequired("search", [&a, this, mustHaveInit](TokenManager& tm){
                a = search(tm, mustHaveInit); });
            p.addFuncRequired("termination", [&b, this](TokenManager& tm){
                b = termination(tm); });
            p.addFuncRequired("perturbation", [&c, this](TokenManager& tm){
                c = perturbation(tm); });
            p.addFuncRequired("acceptance", [&d, this](TokenManager& tm){
                d = acceptance(tm); });

            p.parse(tm);
            p.print();
        } else {
            CheckBrac br(tm);
            a = search(tm, mustHaveInit); // the init will be the init of the search
            b = termination(tm);
            c = perturbation(tm);
            d = acceptance(tm);
        }

        return new emili::MyIteratedLocalSearch(a,b,c,d);
    }
    else if(tm.checkToken(TABU))
    {
        printTab(TABU);
        return tabu(tm);
    }
    else if(tm.checkToken(NO_SEARCH))
    {
        search(tm, mustHaveInit, prefix);
        throw NoSearch();
    }
    else if(mustHaveInit && checkTokenParams(tm, IDENTITY, {"init"}, prefix))
    {
        CheckBrac br(tm);
        auto a = initializer(tm);
        return new emili::IdentityLocalSearch(*a);
    }
    else if(!mustHaveInit && checkTokenParams(tm, IDENTITY, {}, prefix))
    {
        CheckBrac br(tm);
        return new emili::IdentityLocalSearch();
    }
    else if(mustHaveInit && checkTokenParams(tm, FIRST, {"init", "term", "neigh"}, prefix))
    {
        CheckBrac br(tm);
        auto a = initializer(tm);
        auto b = termination(tm);
        auto c = neigh(tm);
        return new emili::FirstImprovementSearch(*a, *b, *c);
    }
    else if(!mustHaveInit && checkTokenParams(tm, FIRST, {"term", "neigh"}, prefix))
    {
        CheckBrac br(tm);
        auto b = termination(tm);
        auto c = neigh(tm);
        return new emili::FirstImprovementSearch(*b, *c);
    }
    else if(mustHaveInit && checkTokenParams(tm, BEST, {"init", "term", "neigh"}, prefix))
    {
        CheckBrac br(tm);
        auto a = initializer(tm);
        auto b = termination(tm);
        auto c = neigh(tm);
        return new emili::BestImprovementSearch(*a, *b, *c);
    }
    else if(!mustHaveInit && checkTokenParams(tm, BEST, {"term", "neigh"}, prefix))
    {
        CheckBrac br(tm);
        auto b = termination(tm);
        auto c = neigh(tm);
        return new emili::BestImprovementSearch(*b, *c);
    }
    else if(mustHaveInit && checkTokenParams(tm, SA, {"init", "neigh", "inittemp", "acceptance", "cooling", "term", "templ"}, prefix))
    {
        CheckBrac br(tm);
        auto initsol = initializer(tm);
        auto nei = neigh(tm);
        return sa.buildSA(tm, initsol, nei);
    }
    else if(!mustHaveInit && checkTokenParams(tm, SA, {"nei", "inittemp", "acceptance", "cooling", "term", "templ"}, prefix))
    {
        CheckBrac br(tm);
        auto initsol = new emili::ExamTT::RandomInitialSolution(instance);
        auto nei = neigh(tm);
        return sa.buildSA(tm, initsol, nei); // initsol will not be used
    }
    else if(mustHaveInit && checkTokenParams(tm, INIT_AND_SEARCH, {"init", "search"}))
    {
        CheckBrac br(tm);
        auto a = initializer(tm);
        auto b = search(tm, false);
        return new emili::InitAndSearch(a,b);
    }
    else if(tm.peek() == SA_BSU){
        printTab(prefix + ": " + SA_BSU);
        tm.next();

        ArgParser p;
        p.addInt("factor", 1);
        p.addInt("factor-exponent", 0);
        p.addInt("freq", 0);
        p.addBool("percent", false);
        p.addFloat("percent-kempe", 0);

        p.parse(tm);
        p.print();
        TabLevel lvl;

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

        // derived
        // double tmin = 0.50;
        constexpr bool wantToUseIrace = false;
        bool useIrace = itmax == baselineitmax && wantToUseIrace;

        int nS = useIrace ? 764141 : (int)(- itmax / (std::log(tr) / std::log(alpha)));
        int nA = useIrace ? 106980 : rho * nS;

        ostringstream oss;

        oss << "nA " << nA << " " << "nS " << nS;
        printTab(oss.str());
        oss.str("");

        // 251_629
        // 111_107

        // 206_134: assign = [[33, 6], [50, 3], [15, 3], [51, 3], [13, 0], [11, 2], [21, 5], [38, 4], [14, 2], [25, 6], [53, 3], [9, 2], [45, 3], [36, 2], [18, 2], [21, 6], [47, 0], [42, 1], [23, 1], [9, 0], [21, 5], [53, 2], [28, 2], [41, 0], [17, 3], [8, 3], [45, 3], [32, 3], [38, 6], [29, 3], [3, 2], [45, 0], [27, 2], [49, 3], [17, 3], [33, 2], [12, 0], [18, 6], [6, 4], [52, 5], [41, 4], [51, 6], [48, 6], [39, 3], [6, 2], [51, 1], [5, 5], [19, 0], [39, 6], [6, 4], [11, 5], [50, 2], [11, 1], [27, 0], [11, 6], [53, 0], [29, 3], [19, 0], [42, 0], [5, 6], [37, 4], [3, 5], [3, 0], [18, 6], [19, 4], [22, 6], [5, 3], [8, 6], [12, 4], [2, 1], [1, 3], [45, 5], [49, 3], [33, 1], [29, 5], [49, 4], [3, 4], [1, 4], [49, 2], [4, 0], [7, 6], [6, 2], [23, 1], [35, 2], [14, 5], [20, 2], [4, 0], [51, 0], [44, 2], [40, 4], [6, 1], [3, 0], [34, 2], [21, 3], [50, 4], [6, 2], [18, 1], [43, 6], [5, 5], [0, 2], [17, 6], [6, 4], [51, 5], [52, 3], [50, 3], [25, 6], [13, 4], [51, 2], [28, 4], [29, 5], [7, 1], [6, 3], [19, 3], [11, 3], [9, 4], [41, 5], [42, 0], [32, 5], [23, 5], [23, 4], [6, 3], [18, 3], [39, 6], [21, 2], [8, 4], [34, 4], [12, 4], [13, 5], [21, 6], [40, 1], [42, 3], [47, 0], [18, 6], [27, 2], [17, 1], [38, 1], [7, 1], [46, 3], [25, 2], [17, 0], [49, 4], [50, 1], [8, 5], [3, 2], [43, 5], [31, 6], [16, 2], [13, 5], [52, 0], [8, 0], [46, 1], [12, 1], [4, 4], [4, 2], [33, 5], [39, 3], [39, 3], [51, 2], [28, 3], [28, 3], [14, 4], [45, 0], [1, 6], [33, 4], [5, 2], [9, 0], [40, 0], [6, 1], [44, 4], [48, 3], [9, 4], [21, 2], [1, 5], [26, 6], [27, 5], [12, 4], [23, 1], [28, 2], [28, 2], [15, 5], [37, 3], [26, 2], [13, 1], [14, 6], [1, 3], [38, 6], [5, 4], [10, 5], [20, 0], [8, 3], [49, 1], [17, 2], [10, 1], [2, 1], [12, 1], [32, 3], [48, 4], [22, 0], [41, 0], [24, 2], [50, 4], [2, 4], [19, 1], [35, 2], [2, 1], [13, 5], [24, 4], [6, 6], [48, 0], [51, 4], [0, 3], [48, 3], [43, 6], [37, 5], [26, 4], [46, 5], [13, 6], [41, 1], [31, 3], [14, 0], [49, 6], [41, 2], [53, 6], [23, 1], [13, 3], [34, 2], [10, 2], [10, 3], [31, 2], [34, 3], [16, 1], [23, 1], [37, 4], [23, 3], [32, 3], [3, 5], [14, 4], [40, 3], [38, 1], [24, 0], [17, 2], [43, 3], [33, 2], [16, 0], [38, 0], [14, 6], [1, 4], [36, 2], [13, 4], [42, 6], [47, 0], [47, 6], [2, 5], [8, 5], [25, 3], [11, 4], [35, 3], [2, 2], [7, 0], [27, 5], [1, 6], [16, 4], [17, 1], [31, 6], [4, 6], [34, 0], [43, 0], [24, 4], [30, 4], [2, 6], [27, 0], [36, 6], [11, 5], [44, 1], [39, 5], [12, 1], [25, 1], [4, 5], [41, 3], [7, 0], [47, 3], [52, 5], [42, 2], [3, 2], [49, 5], [28, 4], [0, 6], [24, 1], [24, 6], [7, 0], [50, 4], [7, 2], [20, 4], [43, 2], [33, 3], [2, 3], [20, 5], [44, 4], [4, 2], [14, 0], [38, 5], [1, 3], [34, 0], [9, 3], [44, 2], [32, 2], [44, 4], [52, 4], [46, 3], [12, 6], [9, 6], [31, 3], [41, 4], [43, 6], [24, 6], [36, 3], [25, 6], [12, 0], [32, 5], [17, 5], [47, 4], [14, 2], [53, 6], [33, 2], [38, 5], [14, 6], [52, 2], [47, 3], [35, 1], [45, 3], [24, 5], [9, 3], [1, 6], [8, 2], [33, 6], [20, 1], [15, 6], [7, 1], [18, 6], [35, 1], [15, 5], [5, 6], [30, 2], [3, 2], [7, 1], [7, 6], [53, 0], [12, 2], [3, 5], [11, 4], [21, 1], [31, 6], [9, 4], [34, 1], [50, 0], [32, 6], [10, 4], [5, 0], [27, 2], [45, 3], [18, 2], [39, 0], [50, 1], [2, 4], [37, 2], [4, 1], [14, 1], [52, 4], [50, 4], [32, 0], [39, 5], [17, 0], [20, 1], [30, 2], [12, 6], [11, 2], [17, 0], [47, 6], [12, 4], [42, 0], [53, 2], [4, 5], [3, 6], [11, 2], [44, 6], [31, 5], [38, 5], [41, 1], [26, 4], [53, 5], [21, 0], [6, 5], [35, 3], [0, 6], [19, 1], [45, 4], [9, 1], [3, 3], [41, 2], [51, 5], [30, 6], [6, 0], [46, 4], [32, 6], [8, 2], [27, 1], [52, 0], [51, 5], [38, 1], [5, 4], [31, 3], [41, 2], [32, 4], [49, 4], [43, 6], [45, 4], [42, 2], [40, 2], [1, 2], [50, 1], [53, 3], [19, 2], [53, 6], [35, 6], [10, 5], [8, 6], [52, 0], [19, 6], [26, 0], [40, 6], [9, 6], [2, 5], [46, 5], [32, 1], [15, 2], [0, 6], [31, 0], [23, 6], [28, 1], [33, 3], [9, 4], [8, 0], [20, 1], [0, 1], [44, 1], [1, 6], [6, 2], [17, 5], [28, 4], [5, 3], [29, 3], [46, 5], [9, 3], [14, 6], [33, 2], [42, 5], [52, 4], [8, 5], [27, 2], [16, 3], [8, 1], [31, 0], [3, 4], [8, 4], [37, 3], [33, 4], [32, 0], [6, 5], [2, 0], [4, 6], [21, 4], [27, 4], [9, 0], [8, 5], [22, 6], [52, 4], [52, 0], [49, 6], [13, 5], [39, 0], [44, 4], [25, 0], [15, 6], [51, 2], [24, 6], [16, 5], [20, 1], [29, 4], [27, 6], [8, 1], [38, 2], [30, 1], [35, 3], [45, 3], [14, 6], [50, 2], [5, 1], [8, 4], [5, 3], [28, 4], [20, 2], [42, 5], [25, 1], [19, 5], [31, 2], [14, 0], [36, 1], [49, 0], [19, 5], [29, 5], [17, 0], [15, 4], [43, 4], [4, 6], [43, 6], [0, 1], [5, 4], [12, 4], [23, 5], [43, 6], [51, 0], [41, 3], [30, 6], [3, 6], [50, 5], [38, 5], [46, 3], [6, 2], [11, 0], [48, 4], [24, 1], [46, 0], [16, 4], [48, 0], [53, 1], [3, 0], [34, 1], [18, 4], [14, 2], [51, 4], [34, 0], [43, 5], [5, 4], [37, 6], [26, 3], [23, 1], [2, 6], [7, 5], [8, 3], [39, 3], [20, 2], [24, 6], [45, 6], [4, 0], [10, 4], [4, 5], [0, 4], [6, 4], [46, 3], [45, 4], [31, 5], [24, 6], [46, 6], [27, 0], [50, 2], [1, 5], [38, 5], [53, 1], [2, 0], [26, 2], [20, 1], [6, 2], [48, 0], [39, 4], [21, 3], [9, 3], [2, 0], [42, 4], [8, 0], [32, 4], [29, 6], [28, 4], [45, 1], [20, 5], [18, 1], [49, 0], [11, 1], [25, 3], [22, 4], [9, 0], [23, 6], [27, 6], [21, 3], [49, 6], [7, 4], [33, 0], [24, 0], [27, 0], [33, 0], [43, 5], [32, 6], [20, 6], [50, 0], [14, 2], [44, 2], [52, 4], [11, 3]]; hard = (455, 7843, 0, 6, 2, 0, 0); soft = (3059, 0, 9062, 9062, 8700, 265, 1560);

        auto init = new emili::ExamTT::RandomInitialSolution(instance);

        emili::Neighborhood* neigh = nullptr;
        if(p.Float("percent-kempe") == 0) {
            neigh = new emili::ExamTT::MixedMoveSwapNeighborhood(instance, sr);
        } else {
            double k = p.Float("percent-kempe") / 100.0;
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
    else if(checkTokenParams(tm, VND, {}, prefix))
    {
        return vnd(tm, mustHaveInit);
    }
    else if(checkTokenParams(tm, ITERATED_GREEDY, {"cons", "termin", "destr", "accept"}, prefix))
    {
        CheckBrac br(tm);
        auto c = constructor(tm); // if must have init, will use constructFull of constructor
        auto t = termination(tm);
        auto d = destructor(tm);
        auto ac = acceptance(tm);
        return new emili::IteratedGreedy(*c,*t,*d,*ac);
    }
    else if(checkTokenParams(tm, MY_ITERATED_GREEDY, {"search", "cons", "termin", "destr", "accept"}, prefix))
    {
        emili::LocalSearch* ls = nullptr;
        emili::Constructor* c = nullptr;
        emili::Termination* t = nullptr;
        emili::Destructor* d = nullptr;
        emili::Acceptance* ac = nullptr;

        if(tm.peekIs("{")) {
            ArgParser p;

            p.addFuncRequired("search", [this, &ls](TokenManager& tm) {
                 ls = search(tm, false); });
            p.addFuncRequired("constructor", [this, &c](TokenManager& tm) {
                c = constructor(tm); });
            p.addFuncRequired("termination", [this, &t](TokenManager& tm) {
                t = termination(tm); });
            p.addFuncRequired("destructor", [this, &d](TokenManager& tm) {
                d = destructor(tm); });
            p.addFuncRequired("acceptance", [this, &ac](TokenManager& tm) {
                ac = acceptance(tm); });

            p.parse(tm);
            // p.print();
        } else {
            CheckBrac br(tm);
            ls = search(tm, false);
            c = constructor(tm); // if mustHaveInit, use constructFull as init
            t = termination(tm);
            d = destructor(tm);
            ac = acceptance(tm);
        }

        return new emili::MyIteratedGreedy(c,t,d,ac,ls);
    }
    else if(tm.checkToken(INFO)) {
        instance.presentation(std::cout);

        throw NoSearch();
    }
    else if(tm.peek() == TEST_DELTA || tm.peek() == TEST_DELTA_REMOVE_ADD) {
        bool isTestDelta = tm.peek() == TEST_DELTA;
        bool isTestDeltaRemoveAdd = tm.peek() == TEST_DELTA_REMOVE_ADD;
        printTab(tm.peek());
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
    else if(tm.checkToken(SHOW_PARAMS)) {
        search(tm, mustHaveInit);
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
            throw ParsingError("ints... are expected");

        const int P = instance.P();
        for(int y : x)
            if(!(0 <= y && y < P))
                throw ParsingError("int " + to_string(y) + " is not a valid period in [0," + to_string(P) + "[");

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
        emili::InitialSolution* ini = initializer(tm);
        emili::Solution* s = ini->generateSolution();
        std::cout << s->getSolutionRepresentation() << std::endl;
        std::cout << s->getSolutionValue() << std::endl;
        std::cerr << s->getSolutionValue() << std::endl;

        throw NoSearch();
    }
    else if(tm.checkToken(INFO_INIT))
    {
        auto ini = initializer(tm);
        auto raw = ini->generateSolution();
        auto s = (emili::ExamTT::ExamTTSolution*) raw;
        s->printTo(instance, std::cout);

        throw NoSearch();
    }
    else if(tm.checkToken(BRUTE)) {
        printTab(prefix + ": " + BRUTE);
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
                    throw ErrorExpected(tm, "]", {"]"});
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

    throw ErrorExpected(tm, "SEARCH", available);
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

emili::Perturbation* ExamTTParser::perturbation(prs::TokenManager& tm, std::string prefix)
{
    static const std::string NO = "no";
    static const std::string NOPERTURB = "nopertub";
    static const std::string RandomMovePerturbationInPlace = "random-move";
    static const std::string VNRandomMovePerturbationInPlace = "variable-neigh-random-move";
    static const std::string constructperturb = "construct-perturb";

    static const std::vector<std::string> available = {
        NO, RandomMovePerturbationInPlace, VNRandomMovePerturbationInPlace,
        constructperturb
    };

    prs::TabLevel level;

    if(tm.checkToken(NOPERTURB) || tm.checkToken(NO)) {
        printTab(prefix + ": no");
        return new emili::NoPerturbation();
    } else if(checkTokenParams(tm, RandomMovePerturbationInPlace, {"rneigh", "steps"}, prefix)) {
        auto n = neigh(tm);
        int s = getIntParam(tm, "s");
        return new emili::RandomMovePerturbationInPlace(*n, s);
    } else if(checkTokenParams(tm, VNRandomMovePerturbationInPlace, {"steps", "iter", "rneigh..."}, prefix)) {
        int steps = getIntParam(tm, "steps");
        int iterations = getIntParam(tm, "iter");
        auto neigh = neighs(tm);
        return new emili::VNRandomMovePerturbationInPlace(neigh, steps, iterations);
    } else if(checkTokenParams(tm, constructperturb, {"construct", "destruct"}, prefix)) {
        auto a = constructor(tm);
        auto b = destructor(tm);
        return new emili::ConstructDestructPertub(a,b);
    }

    throw ErrorExpected(tm, "PERTURBATION", available);
}

emili::Acceptance* ExamTTParser::acceptance(prs::TokenManager& tm, std::string prefix)
{
    prs::TabLevel level;

    ostringstream oss; oss << prefix << ": ";

    if(tm.checkToken(ACCEPTANCE_METRO))
    {
        CheckBrac br(tm);
        float n;

        if(! br.has && tm.checkToken("{")) {
            ArgParser p;
            p.addFloatRequired("temperature");
            p.parse(tm);
            n = p.Float("temperature");
            if(! tm.checkToken("}"))
                throw ErrorExpected(tm, "}", {"}"});
        } else {
            n = tm.getDecimal();
        }

        oss << "metropolis(temperature = " << n << ")";
        printTab(oss.str());
        return new emili::MetropolisAcceptance(n);
    }
    else if(tm.checkToken("metropolis-proba"))
    {
        CheckBrac br(tm);
        float value, proba;
        if(! br.has && tm.checkToken("{")) {
            ArgParser p;
            p.addFloatRequired("value");
            p.addFloat("proba", 0.5f);
            p.parse(tm);
            value = p.Float("value");
            proba = p.Float("proba");
            if(!(0 <= proba && proba <= 1))
                throw ParsingError("proba must be between 0 and 1");
            if(! tm.checkToken("}"))
                throw ErrorExpected(tm, "}", {"}"});
        } else {
            value = tm.getDecimal();
            proba = 0.5f;
        }
        // when the difference in cost = value, the proba is <proba>
        double temperature = value / log(proba);
        oss << "metropolis(temperature = " << temperature << ")";
        return new emili::MetropolisAcceptance(temperature);
    }
    else if(tm.checkToken(ACCEPTANCE_ALWAYS))
    {
        CheckBrac br(tm);
        emili::accept_candidates accc;
        string t1 = *tm;

        if(tm.checkToken(ACCEPTANCE_INTENSIFY))
            accc = emili::ACC_INTENSIFICATION;

        else if(tm.checkToken(ACCEPTANCE_DIVERSIFY))
            accc = emili::ACC_DIVERSIFICATION;

        else
            throw ErrorExpected(tm, ACCEPTANCE_INTENSIFY, {ACCEPTANCE_INTENSIFY, ACCEPTANCE_DIVERSIFY});

        oss << "always(" << t1 << ")";
        printTab(oss.str());
        return new emili::AlwaysAccept(accc);
    }
    else if(tm.checkToken(ACCEPTANCE_IMPROVE))
    {
        CheckBrac br(tm);
        oss << "improve";
        printTab(oss.str());
        return new emili::ImproveAccept();
    }
    else if(tm.checkToken(ACCEPTANCE_SA_METRO))
    {
        CheckBrac br(tm);
        float start =tm.getDecimal();
        float end =tm.getDecimal();
        float ratio =tm.getDecimal();
        oss << "metropolis(start = " << start << ", end = " << end << ", ratio = " << ratio << ")";
        printTab(oss.str());
        return new emili::Metropolis(start,end,ratio);
    }
    else if(tm.checkToken(ACCEPTANCE_PMETRO))
    {
        CheckBrac br(tm);
        float start = tm.getDecimal();
        float end = tm.getDecimal();
        float ratio = tm.getDecimal();
        int iterations = tm.getInteger();
        oss << "metropolis(start = " << start << ", end = " << end << ", ratio = " << ratio << ", frequence = " << iterations << ")";
        printTab(oss.str());
        return new emili::Metropolis(start,end,ratio,iterations);
    }
    else if(tm.checkToken(ACCEPTANCE_SA))
    {
        CheckBrac br(tm);
        float start = tm.getDecimal();
        float end = tm.getDecimal();
        float ratio = tm.getDecimal();
        int iterations = tm.getInteger();
        float alpha = tm.getDecimal();
        oss << "metropolis(start = " << start << ", end = " << end << ", ratio = " << ratio << ", frequence = " << iterations << ", alpha = " << alpha << ")";
        printTab(oss.str());
        return new emili::Metropolis(start,end,ratio,iterations,alpha);
    }
    else if(tm.checkToken(ACCEPTANCE_IMPROVE_PLATEAU))
    {
        CheckBrac br(tm);
        int plateau_steps = tm.getInteger();
        int threshold = tm.getInteger();
        oss << "Accept a diversification solution if it improves on the intensification otherwise it will accept " << plateau_steps << " non improving steps once it reaches the threshold of " << threshold;
        printTab(oss.str());
        return new emili::AcceptPlateau(plateau_steps,threshold);
    }

    throw ErrorExpected(tm, "ACCEPTANCE_CRITERIA", {ACCEPTANCE_METRO, ACCEPTANCE_ALWAYS, ACCEPTANCE_INTENSIFY, ACCEPTANCE_DIVERSIFY, ACCEPTANCE_IMPROVE, ACCEPTANCE_SA_METRO, ACCEPTANCE_PMETRO, ACCEPTANCE_SA, ACCEPTANCE_IMPROVE_PLATEAU});
}

emili::BestTabuSearch* ExamTTParser::tabu(prs::TokenManager& tm)
{
    if(tm.checkToken(BEST))
    {
        auto a = initializer(tm);
        auto b = termination(tm);
        auto c = neigh(tm);
        auto tmem = tabuMemory(c, tm);
        return new emili::BestTabuSearch(*a, *b, *c, *tmem);
    }
    else if(tm.checkToken(FIRST))
    {
        auto a = initializer(tm);
        auto b = termination(tm);
        auto c = neigh(tm);
        auto tmem = tabuMemory(c, tm);
        return new emili::FirstTabuSearch(*a, *b, *c, *tmem);
    }

    throw ErrorExpected(tm, "PIVOTAL_RULE", {BEST,FIRST});
}

emili::TabuMemory* ExamTTParser::tabuMemory(emili::Neighborhood* n,prs::TokenManager& tm)
{
    prs::TabLevel level;

    emili::TabuMemory* tmem = nullptr;

    /*TODO HERE GOES THE CODE TO INSTATIATE A TABUTENURE IMPLEMENTATION*/

    return tmem;
}

emili::LocalSearch* ExamTTParser::vnd(prs::TokenManager& tm, bool mustHaveInit, std::string prefix)
{
    prs::TabLevel level;

    if(mustHaveInit && checkTokenParams(tm, FIRST, {"init", "term", "neigh"}, prefix))
    {
        CheckBrac br(tm);
        auto in = initializer(tm);
        auto te = termination(tm);
        auto nes = neighs(tm);
        return new emili::VNDSearch<emili::FirstImprovementSearch>(*in,*te,nes);
    }
    else if(!mustHaveInit && checkTokenParams(tm, FIRST, {"term", "neigh"}, prefix))
    {
        CheckBrac br(tm);
        auto in = new emili::ExamTT::RandomInitialSolution(instance);
        auto te = termination(tm);
        auto nes = neighs(tm);
        return new emili::VNDSearch<emili::FirstImprovementSearch>(*in,*te,nes);
    }
    else if(mustHaveInit && checkTokenParams(tm, BEST, {"init", "term", "neigh"}, prefix))
    {
        CheckBrac br(tm);
        auto in = initializer(tm);
        auto te = termination(tm);
        auto nes = neighs(tm);
        return new emili::VNDSearch<emili::BestImprovementSearch>(*in,*te,nes);
    }
    else if(!mustHaveInit && checkTokenParams(tm, BEST, {"term", "neigh"}, prefix))
    {
        CheckBrac br(tm);
        auto in = new emili::ExamTT::RandomInitialSolution(instance);
        auto te = termination(tm);
        auto nes = neighs(tm);
        return new emili::VNDSearch<emili::BestImprovementSearch>(*in,*te,nes);
    }

    throw ErrorExpected(tm, "VND_SEARCH", {FIRST, BEST});
}

emili::ExamTT::InsertHeuristic* ExamTTParser::insertHeuristic(prs::TokenManager& tm, std::string prefix) {
    static const std::string BestInsertHeuristic = "best";

    TabLevel l;

    const std::vector<std::string> available = {
        BestInsertHeuristic
    };

    if(checkTokenParams(tm, BestInsertHeuristic, {}, prefix)) {
        CheckBrac br(tm);
        return new emili::ExamTT::BestInsertHeuristic(instance);
    }

    throw ErrorExpected(tm, "INSERT_HEURISTIC", available);
}

emili::Constructor* ExamTTParser::constructor(prs::TokenManager& tm, std::string prefix) {
    static const std::string RandomOrderInserter = "random-order";
    static const std::string DegreeInserter = "degree";

    const std::vector<std::string> available = {
        RandomOrderInserter, DegreeInserter,
    };

    TabLevel l;

    if(checkTokenParams(tm, RandomOrderInserter, {"heur"}, prefix)) {
        CheckBrac br(tm);
        auto heur = insertHeuristic(tm, "heur");
        return new emili::ExamTT::RandomOrderInserter(instance, heur);
    } else if(checkTokenParams(tm, DegreeInserter, {"heur"}, prefix)) {
        CheckBrac br(tm);
        auto heur = insertHeuristic(tm, "heur");
        return new emili::ExamTT::DegreeInserter(instance, heur);
    }

    throw ErrorExpected(tm, "CONSTRUCTOR", available);
}

emili::Destructor* ExamTTParser::destructor(prs::TokenManager& tm, std::string prefix) {
    static const std::string FixedRandomDestructor = "fixed-random";
    static const std::string NRandomDaysDestructor = "random-days";
    static const std::string NGroupedDaysDestructor = "grouped-days";
    static const std::string NGroupedPeriodsDestructor = "grouped-periods";
    static const std::string NBiggestWeightedDestructor = "greedy-biggest";

    TabLevel l;
    ostringstream oss;
    oss << prefix << ": ";

    static const std::vector<std::string> available = {
        FixedRandomDestructor, NRandomDaysDestructor, NGroupedDaysDestructor,
        NGroupedPeriodsDestructor, NBiggestWeightedDestructor
    };

    if(tm.checkToken(FixedRandomDestructor)) {
        CheckBrac br(tm);
        int N = getIntPercentParam(tm, instance.E());
        oss << FixedRandomDestructor << "(N = " << N << ")";
        printTab(oss.str());
        return new emili::ExamTT::FixedRandomDestructor(instance, N);
    } else if(tm.checkToken(NRandomDaysDestructor)) {
        CheckBrac br(tm);
        int N = getIntPercentParam(tm, instance.Days());
        oss << NRandomDaysDestructor << "(N = " << N << ")";
        printTab(oss.str());
        return new emili::ExamTT::NRandomDaysDestructor(instance, N);
    } else if(tm.checkToken(NGroupedDaysDestructor)) {
        CheckBrac br(tm);
        int N = getIntPercentParam(tm, instance.Days());
        oss << NGroupedDaysDestructor << "(N = " << N << ")";
        printTab(oss.str());
        return new emili::ExamTT::NGroupedDaysDestructor(instance, N);
    } else if(tm.checkToken(NGroupedPeriodsDestructor)) {
        CheckBrac br(tm);
        int N = getIntPercentParam(tm, instance.P());
        oss << NGroupedPeriodsDestructor << "(N = " << N << ")";
        printTab(oss.str());
        return new emili::ExamTT::NGroupedPeriodsDestructor(instance, N);
    } else if(tm.checkToken(NBiggestWeightedDestructor)) {
        CheckBrac br(tm);
        int N = getIntPercentParam(tm, instance.E());
        oss << NBiggestWeightedDestructor << "(N = " << N << ")";
        printTab(oss.str());
        return new emili::ExamTT::NBiggestWeightedDestructor(instance, N);
    }

    throw ErrorExpected(tm, "DESTRUCTOR", available);
}

emili::InitialSolution* ExamTTParser::initializer(prs::TokenManager& tm, std::string prefix)
{
    prs::TabLevel level;
    static const std::string INITIAL_CONSTRUCTOR = "initial-constructor";
    static const std::string given_solution = "given-solution";

    static const std::vector<std::string> available = {
        INITIAL_RANDOM, INITIAL_CONSTRUCTOR, given_solution
    };

    if(checkTokenParams(tm, INITIAL_RANDOM, {}, prefix)){
        CheckBrac br(tm);
        return new emili::ExamTT::RandomInitialSolution(instance);
    } else if(checkTokenParams(tm, INITIAL_CONSTRUCTOR, {"constructor"}, prefix)){
        CheckBrac br(tm);
        auto cons = constructor(tm);
        return new emili::ExamTT::ConstructorInitialSolution(instance, cons);
    } else if(checkTokenParams(tm, given_solution, {}, prefix)) {
        CheckBrac br(tm);
        std::vector<int> vec;

        if(tm.checkToken("list-exact")) {
            // [ 1 2 3 4 5 6 ]
            // 1 2 3 4 5 6
            CheckBrac br(tm);
            while(vec.size() != instance.E() * 2)
                vec.push_back(tm.getInteger());
        }
        else if(tm.checkToken("string-like")) {
            // One String containing E*2 numbers
            auto x = tm.nextToken();
            if(!x)
                throw ParsingError("expected string");
            string s = x;
            for(char & c : s)
                if(!('0' <= c && c <= '9'))
                    c = ' ';
            istringstream iss(s);
            vec.push_back(-1);
            while(iss >> vec.back())
                vec.push_back(-1);
            vec.pop_back();
        }
        else if(tm.checkToken("list-like")) {
            // valid: [[1, 2], [3, 4], [5, 6]]
            // valid: [1, 2, 3, 4, 5, 6]
            // valid: [[1 ,2 [3,4] 5 6,,
            // last token must have at least one int
            while(vec.size() < instance.E() * 2) {
                auto x = tm.nextToken();
                if(!x)
                    throw ParsingError("int...");
                std::string s = x;
                for(char & c : s)
                    if(!('0' <= c && c <= '9'))
                        c = ' ';
                istringstream iss(s);
                vec.push_back(-1);
                while(iss >> vec.back())
                    vec.push_back(-1);
                vec.pop_back();
            }
        // } else if(tm.checkToken("file")) {
        } else {
            throw ErrorExpected(tm, "given-solution-type", {"list-like", "list-exact", "string-like"});
        }

        if(vec.size() != instance.E() * 2)
            throw ParsingError("expeccted E * 2 ints, got " + to_string(vec.size()));
        std::vector<std::pair<int,int>> firstAssign;
        for(int i = 0; i < vec.size(); i += 2)
            firstAssign.push_back(make_pair(vec[i], vec[i+1]));

        return new emili::ExamTT::ZeroInitialSolution(instance, firstAssign);
    }

    throw ErrorExpected(tm, "INITIAL_SOLUTION", available);
}

emili::Termination* ExamTTParser::termination(prs::TokenManager& tm, std::string prefix)
{
    prs::TabLevel level;

    ostringstream oss;
    oss << prefix << ": ";

    static const std::vector<std::string> available = {
        TERMINATION_LOCMIN, TERMINATION_WTRUE, TERMINATION_TIME, TERMINATION_MAXSTEPS
    };

    if(tm.checkToken(TERMINATION_LOCMIN))
    {
        CheckBrac br(tm);
        oss << TERMINATION_LOCMIN;
        printTab(oss.str());
        return new emili::LocalMinimaTermination();
    }
    else if(tm.checkToken(TERMINATION_WTRUE))
    {
        CheckBrac br(tm);
        oss << TERMINATION_WTRUE;
        printTab(oss.str());
        return new emili::WhileTrueTermination();
    }
    else if(tm.checkToken(TERMINATION_TIME))
    {
        CheckBrac br(tm);
        float time = tm.getDecimal();
        if(time == 0)
            time = 1;

        oss << "timed(ratio = " << time << ")";
        printTab(oss.str());
        return new emili::TimedTermination(time);
    }
    else if(tm.checkToken(TERMINATION_MAXSTEPS))
    {
        CheckBrac br(tm);
        int steps = tm.getInteger();
        oss << "maxSteps(steps = " << steps << ")";
        printTab(oss.str());
        return new emili::MaxStepsTermination(steps);
    }
    else if(tm.checkToken(TERMINATION_MAXSTEPS_DEBUG))
    {
        int steps;
        double percent = 10;
        std::string prefix;
        if(tm.peekIs("{")) {
            ArgParser p;
            p.addIntRequired("steps");
            p.addInt("percent", 10);
            p.addString("prefix", "");
            p.parse(tm);
            oss << "maxSteps";
            p.printInline(oss.str());
            steps = p.Int("steps");
            percent = p.Int("percent");
            prefix = p.String("prefix");
        } else {
            CheckBrac br(tm);
            steps = tm.getInteger();
            percent = 10;
            prefix = "";
            oss << "maxSteps(steps = " << steps << ", debugPercent = " << percent << ")";
            printTab(oss.str());
        }
        auto x = new emili::MaxStepsTerminationDebug(steps, percent);
        x->setPrefix(prefix);
        return x;
    }

    throw ErrorExpected(tm, "TERMINATION_CRITERIA", available);
}

emili::Neighborhood* ExamTTParser::neigh(prs::TokenManager& tm, bool errorIfNotFound, string prefix)
{
    const static std::string MoveNeighborhood = "move";
    const static std::string SwapNeighborhood = "swap";
    const static std::string KempeChainNeighborhoodFastIter = "kempe";

    const static std::vector<std::string> available = {
        MoveNeighborhood, SwapNeighborhood, KempeChainNeighborhoodFastIter
    };

    prs::TabLevel level;

    ostringstream oss;
    oss << prefix << ": ";

    if(tm.checkToken(MoveNeighborhood)) {
        oss << MoveNeighborhood;
        printTab(oss.str());
        return new emili::ExamTT::MoveNeighborhood(instance);
    }
    else if(tm.checkToken(SwapNeighborhood)) {
        oss << SwapNeighborhood;
        printTab(oss.str());
        return new emili::ExamTT::SwapNeighborhood(instance);
    }
    else if(tm.checkToken(KempeChainNeighborhoodFastIter)) {
        oss << KempeChainNeighborhoodFastIter;
        printTab(oss.str());
        return new emili::ExamTT::KempeChainNeighborhoodFastIter(instance);
    }

    if(errorIfNotFound)
        throw ErrorExpected(tm, "NEIGHBORHOOD", available);

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
        throw ParsingError("neigh+");

    return vnds;
}

void ExamTTParser::problem(prs::TokenManager& tm)
{
    tm.next(); // skip the problem "ExamTT" without reading it
    instanceFilename = tm.tokenAt(1);
    emili::ExamTT::InstanceParser parser(instanceFilename);
    parser.parse(instance);
    globalParams(tm);
}

void ExamTTParser::globalParams(TokenManager &tm) {
    ArgParser p;
    p.addString("hard-weight", "sa-bsu-1");
    p.addString("output-hard-weight", "same");

    // p.addBool("search", true);
    // p.addInt("time-seconds", 0);
    // p.addInt("time-ratio-size", 0);
    // p.addInt("random-seed", 0);
    //      error if p.has("time-seconds") && p.has("time-ratio-size")
    // else parse if p.has("time-seconds") || p.has("time-ratio-size")
    // else do nothing

    p.parseSequence(tm);
    p.print();
    auto x = p.String("hard-weight");
    auto y = p.String("output-hard-weight");
    if(y == "same")
        y = x;

    auto checkHardWeightParam = [](std::string s) -> bool {
        return s == "sa-bsu-1" || s == "sa-bsu-2" || all_of(
            s.begin(), s.end(), [](char c){ return ::isdigit(c) || c == '.'; }
        );
    };

    if(! checkHardWeightParam(x))
        throw ParsingError("hard-weight expected, got '" + x + "'");

    if(! checkHardWeightParam(y))
        throw ParsingError("hard-weight expected, got '" + y + "'");

    // compute hardWeight
    instance.hardWeight =
        x == "sa-bsu-1" ? getHardweightFromFeatures(true) :
        x == "sa-bsu-2" ? getHardweightFromFeatures(false) :
        atof(x.c_str()); // float

    // compute hardWeightFinal
    instance.hardWeightFinal =
        y == "sa-bsu-1" ? getHardweightFromFeatures(true) :
        y == "sa-bsu-2" ? getHardweightFromFeatures(false) :
        atof(y.c_str()); // float

    printTab("Hard weight: " + to_string(instance.hardWeight));
}


emili::LocalSearch* ExamTTParser::buildAlgo(prs::TokenManager& tm)
{
    problem(tm);
    emili::LocalSearch* local = search(tm, true);
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

std::vector<string> ExamTTParser::availableProblemsList() const {
    return std::vector<string>(PROBLEMS_DEF.begin(), PROBLEMS_DEF.end());
}

} // namespace ExamTT
} // namespace prs
