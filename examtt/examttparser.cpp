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
        return startingSolution->clone(); // do we have to clone ?
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

emili::LocalSearch* ExamTTParser::eparams(prs::TokenManager& tm)
{
    prs::TabLevel level;

    if(tm.checkToken(ILS))
        return ils(tm);
    else
        return search(tm);
}


emili::LocalSearch* ExamTTParser::search(prs::TokenManager& tm)
{
    const std::string SA_BSU = "sa_bsu";
    const std::string INTERACTIVE = "interactive";
    const std::string TEST_KEMPE = "test_kempe";
    const std::string TEST_DELTA = "test_delta";
    const std::string TEST_DELTA_REMOVE_ADD = "test_delta_remove_add";
    const std::string INFO = "info";
    const std::string BRUTE = "brute";

    const std::vector<std::string> available = {
        ILS, TABU, FIRST, BEST, SA_BSU, VND,
        BRUTE, INFO, INTERACTIVE, TEST_INIT, TEST_KEMPE, TEST_DELTA, TEST_DELTA_REMOVE_ADD
    };

    prs::TabLevel level;

    if(tm.checkToken(ILS))
        return ils(tm);

    else if(tm.checkToken(TABU))
    {
        printTab(TABU);
        return tparams(tm);
    }
    else if(tm.checkToken(FIRST))
    {
        printTab(FIRST + " (3) init term neigh");
        auto p = params(tm);
        return new emili::FirstImprovementSearch(*get<0>(p), *get<1>(p), *get<2>(p));
    }
    else if(tm.checkToken(BEST))
    {
        printTab(BEST + " (3) init term neigh");
        auto p = params(tm);
        return new emili::BestImprovementSearch(*get<0>(p), *get<1>(p), *get<2>(p));
    }
    else if(tm.checkToken(SA)) {
        printTab(SA + " (7) init nei inittemp acceptance cooling term templ");

        auto initsol = init(tm);
        auto nei = neigh(tm);
        return sa.buildSA(tm, initsol, nei);
    }
    else if(tm.peek() == SA_BSU){
        printTab(SA_BSU);
        tm.next();

        int factor = 1;
        int freq = 0;
        bool percent = true;

        for(;;) {
            if(tm.checkToken("factor"))
                factor = tm.getInteger();
            else if(tm.checkToken("freq"))
                freq = tm.getInteger();
            else if(tm.checkToken("no-percent"))
                percent = false;
            else
                break;
        }

        cout << setw(20) << "factor " << factor << endl;
        cout << setw(20) << "freq " << freq << endl;
        cout << setw(20) << "percent " << boolalpha << percent << endl;

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
        auto neigh = new emili::ExamTT::MixedMoveSwapNeighborhood(instance, sr);
        auto acc = new SAMetropolisAcceptance(t0);

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

        auto explo = new SARandomExplorationNoCopy(neigh, acc, term);

        auto sa = new SimulatedAnnealing(init, initialTemperature, acc, cooling, tempRestart, term, tempLength, explo, neigh);

        if(percent) {
            ((SAMaxIterTerminationDebug*)term)->setOnEachPercent(freq, [sa,explo,acc]{
                cout << std::setprecision(10) << std::fixed <<
                    "&bestcost=" << sa->status->best_cost <<
                    "&curcost=" << sa->bestSoFar->getSolutionValue() <<
                    "&temp=" << sa->status->temp <<
                    // "&nacc=" << explo->nacc <<
                    // "&nnacc=" << explo->nnacc <<
                    "&tot=" << sa->status->total_counter <<
                    // "&sumprob=" << acc->sumprob <<
                    // "&nnotimproving=" << acc->nnotimproving <<
                endl;

                // explo->nacc = 0;
                // explo->nnacc = 0;
                // acc->nnotimproving = 0;
                // acc->sumprob = 0;
            });
        }

        return sa;
    }
    else if(tm.checkToken(VND))
    {
        return vparams(tm);
    }
    else if(tm.checkToken(INFO)) {
        instance.presentation(std::cout);

        throw NoSearch();
    }
    else if(tm.peek() == TEST_DELTA || tm.peek() == TEST_DELTA_REMOVE_ADD) {
        bool isTestDelta = tm.peek() == TEST_DELTA;
        bool isTestDeltaRemoveAdd = tm.peek() == TEST_DELTA_REMOVE_ADD;
        tm.next();

        int N = 1000;
        int G = 1;
        bool checkEachMove = true;

        for(;;)
            if(tm.checkToken("N"))
                N = tm.getInteger();
            else if(tm.checkToken("no-each-move"))
                checkEachMove = false;
            else if(tm.checkToken("each-move"))
                checkEachMove = true;
            else if(isTestDeltaRemoveAdd && tm.checkToken("G"))
                G = tm.getInteger();
            else
                break;

        cout
            << setw(20) << "N " << N << endl
            << setw(20) << "each-move " << boolalpha << checkEachMove << endl
        ;

        if(isTestDeltaRemoveAdd)
            cout
                << setw(20) << "G " << G << endl
            ;

        if(isTestDelta)
            emili::ExamTT::test::delta(instance, N, checkEachMove);
        else if(isTestDeltaRemoveAdd)
            emili::ExamTT::test::deltaRemoveAdd(instance, N, checkEachMove, G);

        throw NoSearch();
    }
    else if(tm.checkToken(INTERACTIVE)) {
        emili::ExamTT::test::interactive(instance);
        throw NoSearch();
    }
    else if(tm.checkToken(TEST_KEMPE)) {

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

        emili::ExamTT::test::kempe(instance, x);

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
    else
    {
        errorExpected(tm, "SEARCH", available);
        return nullptr;
    }
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

emili::LocalSearch* ExamTTParser::ils(prs::TokenManager& tm)
{
    printTab("ILS (4) search term perturb acc");
    auto a = search(tm);
    auto b = term(tm);
    auto c = per(tm);
    auto d = acc(tm);
    return new emili::IteratedLocalSearch(*a,*b,*c,*d);
}

emili::Perturbation* ExamTTParser::per(prs::TokenManager& tm)
{
    printTab(*tm);

    prs::TabLevel level;

    if(tm.checkToken("noperturb"))
        return new emili::NoPerturbation();

    errorExpected(tm, "PERTURBATION", {"noperturb"});
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
    else
    {
        errorExpected(tm, "ACCEPTANCE_CRITERIA", {ACCEPTANCE_METRO, ACCEPTANCE_ALWAYS, ACCEPTANCE_INTENSIFY, ACCEPTANCE_DIVERSIFY, ACCEPTANCE_IMPROVE, ACCEPTANCE_SA_METRO, ACCEPTANCE_PMETRO, ACCEPTANCE_SA, ACCEPTANCE_IMPROVE_PLATEAU});
        return nullptr;
    }
}

emili::BestTabuSearch* ExamTTParser::tparams(prs::TokenManager& tm)
{
    if(tm.checkToken(BEST))
    {
        auto p = params(tm);
        auto tmem = tmemory(get<2>(p), tm);
        return new emili::BestTabuSearch(*get<0>(p), *get<1>(p), *get<2>(p), *tmem);
    }
    else if(tm.checkToken(FIRST))
    {
        auto p = params(tm);
        auto tmem = tmemory(get<2>(p), tm);
        return new emili::FirstTabuSearch(*get<0>(p), *get<1>(p), *get<2>(p), *tmem);
    }
    else
    {
        errorExpected(tm, "PIVOTAL_RULE", {BEST,FIRST});
        return nullptr;
    }
}

emili::TabuMemory* ExamTTParser::tmemory(emili::Neighborhood* n,prs::TokenManager& tm)
{
    prs::TabLevel level;

    emili::TabuMemory* tmem = nullptr;

    /*TODO HERE GOES THE CODE TO INSTATIATE A TABUTENURE IMPLEMENTATION*/

    return tmem;
}

std::tuple<emili::InitialSolution*, emili::Termination*, emili::Neighborhood*> ExamTTParser::params(prs::TokenManager& tm)
{
    auto in = init(tm);
    auto te = term(tm);
    auto ne = neigh(tm);
    return std::make_tuple(in, te, ne);
}

emili::LocalSearch* ExamTTParser::vparams(prs::TokenManager& tm)
{
    printTab("VND SEARCH");

    prs::TabLevel level;

    if(tm.checkToken(FIRST))
    {
        printTab(FIRST + " (3) init term neigh");
        auto in = init(tm);
        auto te = term(tm);
        auto nes = neighs(tm);
        return new emili::VNDSearch<emili::FirstImprovementSearch>(*in,*te,nes);
    }
    else if(tm.checkToken(BEST))
    {
        printTab(BEST + " (3) init term neigh");
        auto in = init(tm);
        auto te = term(tm);
        auto nes = neighs(tm);
        return new emili::VNDSearch<emili::BestImprovementSearch>(*in,*te,nes);
    }
    else
    {
        errorExpected(tm, "VND_SEARCH", {FIRST, BEST});
        return nullptr;
    }
}

emili::InitialSolution* ExamTTParser::init(prs::TokenManager& tm)
{
    prs::TabLevel level;

    printTab(*tm);

    if(tm.checkToken(INITIAL_RANDOM)) {
        return new emili::ExamTT::RandomInitialSolution(instance);
    } else {
        errorExpected(tm, "INITIAL_SOLUTION", {INITIAL_RANDOM});
        return nullptr;
    }
}

emili::Termination* ExamTTParser::term(prs::TokenManager& tm)
{
    prs::TabLevel level;

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
    else
    {
        errorExpected(tm, "TERMINATION_CRITERIA", {TERMINATION_LOCMIN, TERMINATION_WTRUE, TERMINATION_TIME, TERMINATION_MAXSTEPS});
        return nullptr;
    }
}

emili::Neighborhood* ExamTTParser::neigh(prs::TokenManager& tm)
{
    prs::TabLevel level;

    printTab(*tm);

    if(tm.checkToken("move")) {
        return new emili::ExamTT::MoveNeighborhood(instance);
    }
    else if(tm.checkToken("swap")) {
        return new emili::ExamTT::SwapNeighborhood(instance);
    }
    else if(tm.checkToken("kempe")) {
        return new emili::ExamTT::KempeChainNeighborhood(instance);
    }
    else {
        errorExpected(tm, "NEIGHBORHOOD", {"move", "swap", "kempe"});
        return nullptr;
    }
}

emili::Neighborhood* ExamTTParser::neighV(prs::TokenManager& tm)
{
    prs::TabLevel level;

    if(tm.checkToken("NEIGH")){
        return nullptr;
    }
    else {
        return nullptr;
    }
}

std::vector<emili::Neighborhood*> ExamTTParser::neighs(prs::TokenManager& tm)
{
    std::vector<emili::Neighborhood*> vnds = {neigh(tm)};
    neighs1(tm, vnds);
    return vnds;
}

void ExamTTParser::neighs1(prs::TokenManager& tm, std::vector<emili::Neighborhood*> & nes)
{
    emili::Neighborhood* n = neighV(tm);
    if(n != nullptr)
    {
        nes.push_back(n);
        neighs1(tm, nes);
    }
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
    emili::LocalSearch* local = eparams(tm);
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
