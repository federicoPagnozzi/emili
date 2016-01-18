#include "examttparser.h"
#include <ctime>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <set>
#include <string>

namespace prs {
namespace ExamTT {

const std::string

/* Algos */
    IG = "ig",
    ILS = "ils",
    TABU = "tabu",
    FIRST = "first",
    BEST = "best",
    VND = "vnd",
    TEST_INIT = "stin",

/* initial solution heuristics */
    INITIAL_RANDOM = "random",

/* Termination criteria*/
    TERMINATION_MAXSTEPS = "maxstep",
    TERMINATION_TIME = "time",
    TERMINATION_LOCMIN = "locmin",
    TERMINATION_ITERA = "iteration",
    TERMINATION_WTRUE = "true",
    TERMINATION_SOA = "soater",

/* acceptance criteria*/
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

using namespace std;

std::string ExamTTParser::info()
{
    ostringstream oss;
    oss << "Usage:\n\n";
    oss << "EMILI INSTANCE_FILE_PATH ExamTT <LOCAL_SEARCH | ITERATED_LOCAL_SEARCH | TABU_SEARCH | VND_SEARCH> [rnds seed]\n\n";
    return oss.str();
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

emili::InitialSolution* in = nullptr;
emili::Termination* te = nullptr;
emili::TabuMemory* tmem = nullptr;
emili::Termination* ilt = nullptr;
emili::Neighborhood* ne = nullptr;
std::vector< emili::Neighborhood*> nes;

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
        params(tm);
        return new emili::FirstImprovementSearch(*in,*te,*ne);
    }
    else if(tm.checkToken(BEST))
    {
        printTab(BEST + " (3) init term neigh");
        params(tm);
        return new emili::BestImprovementSearch(*in,*te,*ne);
    }
    else if(tm.checkToken(VND))
    {
        return vparams(tm);
    }
    else if(tm.checkToken(TEST_INIT))
    {
        emili::InitialSolution* ini = init(tm);
        emili::Solution* s = ini->generateSolution();
        std::cout << s->getSolutionRepresentation() << std::endl;
        std::cout << s->getSolutionValue() << std::endl;
        std::cerr << s->getSolutionValue() << std::endl;
        exit(-1);
        return nullptr;
    }
    else
    {
        errorExpected(tm, "SEARCH", {ILS, TABU, FIRST, BEST, VND, TEST_INIT});
        return nullptr;
    }
}

emili::LocalSearch* ExamTTParser::ils(prs::TokenManager& tm)
{
    printTab("ILS (3) search term perturb");
    /*Parse the inner local search*/
    emili::LocalSearch* ls = search(tm);
    /*Parse the termination*/
    emili::Termination* pft = term(tm);
    /*Parse the perturbation*/
    emili::Perturbation* prsp = per(tm);
    /*Parse the Acceptance*/
    emili::Acceptance* tac = acc(tm);
    /*Build the ils and returns it*/
    emili::LocalSearch* iils = new emili::IteratedLocalSearch(*ls,*pft,*prsp,*tac);   
    return iils;
}

emili::Perturbation* ExamTTParser::per(prs::TokenManager& tm)
{
    printTab(*tm);

    prs::TabLevel level;

    if(tm.checkToken("noperturb"))
        return new emili::NoPertubation();

    errorExpected(tm, "PERTURBATION", {"noperturb"});
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
        params(tm);
        tmem = tmemory(ne,tm);
        return new emili::BestTabuSearch(*in,*te,*ne,*tmem);
    }
    else if(tm.checkToken(FIRST))
    {
        params(tm);
        tmem = tmemory(ne,tm);
        return new emili::FirstTabuSearch(*in,*te,*ne,*tmem);
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

void ExamTTParser::params(prs::TokenManager& tm)
{
    in = init(tm);
    te = term(tm);
    ne = neigh(tm);
}

emili::LocalSearch* ExamTTParser::vparams(prs::TokenManager& tm)
{
    printTab("VND SEARCH");

    prs::TabLevel level;

    if(tm.checkToken(FIRST))
    {
        printTab(FIRST + " (3) init term neigh");
        in = init(tm);
        te = term(tm);
        neighs(tm);
        return new emili::VNDSearch<emili::FirstImprovementSearch>(*in,*te,nes);
    }
    else if(tm.checkToken(BEST))
    {
        printTab(BEST + " (3) init term neigh");
        in = init(tm);
        te = term(tm);
        neighs(tm);
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

    if(tm.checkToken("move")){
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
        /*TODO HERE INSERT THE CODE FOR PARSING THE NEIGHBORHOODS FOR VND */
        return nullptr;
    }
    else
    {  
        return nullptr;
    }
}

void ExamTTParser::neighs(prs::TokenManager& tm)
{
    std::vector<emili::Neighborhood*> vnds;
    vnds.push_back(neigh(tm));
    nes = vnds;
    neighs1(tm);
}

void ExamTTParser::neighs1(prs::TokenManager& tm)
{
    emili::Neighborhood* n = neighV(tm);
    if(n != nullptr)
    {
        nes.push_back(n);
        neighs1(tm);
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
