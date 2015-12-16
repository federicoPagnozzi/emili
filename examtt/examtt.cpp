#include "examtt.h"

#include <iomanip>
#include <sstream>

using namespace std;

namespace emili
{
namespace ExamTT
{

template<typename Type, unsigned N, unsigned Last>
struct tuple_printer {

    static void print(std::ostream& out, const Type& value) {
        out << std::get<N>(value) << ", ";
        tuple_printer<Type, N + 1, Last>::print(out, value);
    }
};

template<typename Type, unsigned N>
struct tuple_printer<Type, N, N> {

    static void print(std::ostream& out, const Type& value) {
        out << std::get<N>(value);
    }

};

template<typename... Types>
std::ostream& operator<<(std::ostream& out, const std::tuple<Types...>& value) {
    out << "(";
    tuple_printer<std::tuple<Types...>, 0, sizeof...(Types) - 1>::print(out, value);
    out << ")";
    return out;
}

template <typename U, typename V>
std::ostream& operator <<(std::ostream& out, const std::pair<U,V>& value) {
    return out << "(" << value.first << ", " << value.second << ")";
}

template <typename T, typename U>
struct KeyValue {
    T key;
    U value;
};

template <typename T, typename U>
bool operator <(KeyValue<T,U> const& a, KeyValue<T,U> const& b) {
    return make_pair(a.key,a.value) < make_pair(b.key, b.value);
}

template <typename T, typename U>
std::ostream& operator<<(std::ostream& out, KeyValue<T,U> const& t) {
    return out << t.key << "(" << t.value << ")";
}

template <typename T>
struct SetPrinter {
    std::set<T> const& self;
    int unsigned flags;
};

template <typename T>
struct VectorPrinter {
    std::vector<T> const& self;
    unsigned int flags;
};

// default : [1, 3, 4, 5]
constexpr unsigned int
    nocomma = 1 << 0,  // [1 3 4 5]
    nobraces = 1 << 1, // 1, 3, 4, 5
    nospaces = 1 << 2; // [1,3,4,5]

template <typename T>
SetPrinter<T> print(std::set<T> const& t) {
    return {t, 0};
}

template <typename T>
SetPrinter<T> print(std::set<T> const& t, unsigned int flags) {
    return {t, flags};
}

template <typename T>
VectorPrinter<T> print(std::vector<T> const& t) {
    return {t, 0};
}

template <typename T>
VectorPrinter<T> print(std::vector<T> const& t, unsigned int flags) {
    return {t, flags};
}

template <typename T>
std::ostream& operator<<(std::ostream& out, SetPrinter<T> printer) {
    auto const& t = printer.self;
    auto F = printer.flags;

    char const* sep = F & nocomma ? (F & nospaces ? "" : " ") :
                                    (F & nospaces ? "," : ", ");
    char const* beg = F & nobraces ? "" : "{";
    char const* end = F & nobraces ? " " : "}";

    out << beg;
    auto it = t.begin();
    if(it != t.end()) {
        out << *it++;
        while(it != t.end())
            out << sep << *it++;
    }
    return out << end;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, VectorPrinter<T> printer) {
    auto const& t = printer.self;
    auto F = printer.flags;

    char const* sep = F & nocomma ? (F & nospaces ? "" : " ") :
                                    (F & nospaces ? "," : ", ");
    char const* beg = F & nobraces ? "" : "[";
    char const* end = F & nobraces ? " " : "]";

    out << beg;
    auto it = t.begin();
    if(it != t.end()) {
        out << *it++;
        while(it != t.end())
            out << sep << *it++;
    }
    return out << end;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, std::set<T> const& t) {
    return out << print(t);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> const& t) {
    return out << print(t);
}


bool isDigit(char c) {
    return '0' <= c && c <= '9';
}

bool isDigit(std::string const& s) {
    if(s.empty())
        return false;

    for(char c : s)
        if(! isDigit(c))
            return false;

    return true;
}

struct CommaSpaceSeparated {
    const std::string s;
    size_t i = 0;
    CommaSpaceSeparated(const std::string s_) : s(s_) {
        skip();
    }

    void skip() {
        while(i < s.size() && (s[i] == ' ' || s[i] == ','))
            i++;
    }

    bool hasNext() {
        return i < s.size();
    }

    std::string next() {
        if(!hasNext())
            throw invalid_argument("no char to read");
        int a = i;
        while(i < s.size() && !(s[i] == ' ' || s[i] == ','))
            i++;
        string r = s.substr(a, i-a);
        skip();
        return r;
    }

    int nextInt() {
        return stoi(next());
    }
};


// 31:12:2020 => 2020 December 31 == (2020, 12, 31)
Date makeDateFromITCFormat(std::string s) {
    if(s.size() == 10) {
        if(s[2] == ':' && s[5] == ':')
            if(isDigit(s.substr(0,2)) && isDigit(s.substr(3,2)) && isDigit(s.substr(6,4)))
                return Date(stoi(s.substr(6,4)), stoi(s.substr(3,2)), stoi(s.substr(0,2)) );
    }
    throw invalid_argument("date has no known format");
}

// 23:05:02 => 23h05 and 2 seconds == (23, 5, 2)
// 23:05 => 23h05 and 0 seconds == (23, 5, 0)
Time makeTimeFromITCFormat(std::string s) {
    if(s.size() == 8) {
        if(s[2] == ':' && s[5] == ':')
            if(isDigit(s.substr(0,2)) && isDigit(s.substr(3,2)) && isDigit(s.substr(6,2)))
                return Time(stoi(s.substr(0,2)), stoi(s.substr(3,2)), stoi(s.substr(6,2)) );
    }
    else if(s.size() == 6) {
        if(s[2] == ':')
            if(isDigit(s.substr(0,2)) && isDigit(s.substr(3,2)))
                return Time(stoi(s.substr(0,2)), stoi(s.substr(3,2)), 0);
    }

    throw invalid_argument("time has no known format");
}

constexpr unsigned int
    nopadding = 1 << 0;

struct Prefix {
    const std::string S;
    const int padding; // < 0 for no padding

    std::string operator ()(int i, unsigned int F) {
        if(F & nopadding)
            return S + to_string(i);
        return operator ()(i);
    }

    std::string operator()(int i) const {
        if(padding < 0)
            return S + to_string(i);

        string base = S + to_string(i);
        return string(max(0, padding - (int)base.size()), ' ') + base;
    }
};

std::ostream& operator <<(std::ostream& out, Prefix const& p) {
    return out << p.S;
}

const Prefix EPrefix{"e", 4};
const Prefix PPrefix{"p", 3};
const Prefix RPrefix{"r", 3};

bool ExamTT::periodInTheEnd(PeriodId p) const {
    return ! correctPeriod(p + institutionalWeightings.frontload.time);
}

int ExamTT::sizeOfExam(ExamId e) const {
    return exams[e].students.size();
}

void ExamTT::buildStructures()  {
    int E = exams.size();

    examsRelated.resize(E);

    for(Two<ExamId> p : examsCoincidence) {
        coincidenceOfExam(p.first).insert(p.second);
        coincidenceOfExam(p.second).insert(p.first);
    }

    for(Two<ExamId> p : examsExclusion) {
        exclusionOfExam(p.first).insert(p.second);
        exclusionOfExam(p.second).insert(p.first);
    }

    for(Two<ExamId> p : examsAfter) {
        afterOfExam(p.first).insert(p.second);
        beforeOfExam(p.second).insert(p.first);
    }

    for(ExamId i = 0; i < E; i++)
    for(ExamId j = i + 1; j < E; j++) {
        if(numberStudentsInCommon[i][j]) {
            hasStudentsInCommonOfExam(i).insert(j);
            hasStudentsInCommonOfExam(j).insert(i);
        }
    }

    vector<ExamId> examsBySize(E);
    for(ExamId i = 0; i < E; i++)
        examsBySize[i] = i;

    std::stable_sort(examsBySize.begin(), examsBySize.end(), [this](ExamId i, ExamId j){
        return sizeOfExam(i) > sizeOfExam(j);
    });

    int M = min(E, institutionalWeightings.frontload.largestExams);
    frontLoadExams.assign(examsBySize.begin(), examsBySize.begin() + M);

    // isInFrontLoad
    isInFrontLoad.resize(E, false);
    for(ExamId i : frontLoadExams)
        isInFrontLoad[i] = true;

    // examIsRoomExclusive
    examIsRoomExclusive.resize(E, false);
    for(ExamId i : examsRoomsExclusive)
        examIsRoomExclusive[i] = true;
}

bool ExamTT::correctExam(int i) const {
    return 0 <= i && i < (int)exams.size();
}

bool ExamTT::correctPeriod(int i) const {
    return 0 <= i && i < (int)periods.size();
}

bool ExamTT::correctRoom(int i) const {
    return 0 <= i && i < (int)rooms.size();
}

bool ExamTT::correctIndexes() const {
    for(std::vector<Two<ExamId>> const* vec : {&examsCoincidence, &examsExclusion, &examsAfter})
        for(Two<ExamId>p : *vec)
            if(! correctExam(p.first) || ! correctExam(p.second))
                return false;

    for(ExamId i : examsRoomsExclusive)
        if(! correctExam(i))
            return false;

    return true;
}

bool ExamTT::nonNegativePenalties() {
    InstitutionalWeightings& w = institutionalWeightings;

    return (
        w.twoInARow >= 0 &&
        w.twoInADay >= 0 &&
        w.periodSpread >= 0 &&
        w.nonMixedDurations >= 0 &&
        w.frontload.penalty >= 0
    ) && all_of(rooms.begin(), rooms.end(), [](Room& x){
        return x.penalty >= 0;
    }) && all_of(periods.begin(), periods.end(), [](Period& x){
        return x.penalty >= 0;
    });
}

void ExamTT::compute() {
    int E = exams.size();

    for(Exam& exam : exams)
        students.insert(exam.students.begin(), exam.students.end());
    meta.Students = students.size();

    numberStudentsInCommon.assign(E, vector<ExamId>(E,0));

    for(ExamId i = 0; i < E; i++)
    for(ExamId j = i + 1; j < E; ++j)
        numberStudentsInCommon[i][j] = numberStudentsInCommon[j][i] = set_intersection(
            exams[i].students.begin(), exams[i].students.end(),
            exams[j].students.begin(), exams[j].students.end(),
            CountIterator()
        );


    int c = 0;
    for(ExamId i = 0; i < E; i++)
    for(ExamId j = i + 1; j < E; ++j)
    if(numberStudentsInCommon[i][j])
        c++;

    meta.ConflictDensity = 2.0 * c / (E * E);

    // number of components of graph of student in common
    std::set<ExamId> toDo;
    for(ExamId i = 0; i < E; i++)
        toDo.insert(i);
    std::set<ExamId> currentComponent;

    function<void(ExamId)> removeNeighbors = [this, &toDo, E, &removeNeighbors, &currentComponent](ExamId i){
        auto it = toDo.find(i);
        if(it != toDo.end()) {
            toDo.erase(it);
            currentComponent.insert(i);
            for(ExamId j = 0; j < E; j++)
            if(numberStudentsInCommon[i][j])
                removeNeighbors(j);
        }
    };

    meta.conflictsComponents = 0;
    while(toDo.size()) {
        currentComponent.clear();
        removeNeighbors(*toDo.begin());
        meta.conflictsComponents++;

        int c = 0;
        for(int i : currentComponent)
        for(int j : currentComponent)
        if(i != j && numberStudentsInCommon[i][j])
            c++;

        meta.conflictsComponentsDensities.push_back(float(c) / (currentComponent.size() * currentComponent.size()));
    }
}

bool ExamTT::periodsSorted() const {
    return is_sorted(periods.begin(), periods.end(), &Period::compareDateTime);
}

bool ExamTT::allPeriodsConsecutives() const {
    for(size_t i = 0; i < periods.size() - 1; i++)
        if(periods[i].date == periods[i+1].date)
            if(seconds(periods[i].time) + periods[i].duration * 60 != seconds(periods[i+1].time))
                return false;

    return true;
}

Two<Room&> ExamTT::roomsOf(Two<RoomId> p) {
    return {rooms[p.first], rooms[p.second]};
}

Two<Exam&> ExamTT::examsOf(Two<ExamId> p) {
    return {exams[p.first], exams[p.second]};
}

Two<Period&> ExamTT::periodsOf(Two<PeriodId> p) {
    return {periods[p.first], periods[p.second]};
}

Two<Room&> ExamTT::roomsOf(RoomId a, RoomId b) {
    return {rooms[a], rooms[b]};
}

Two<Exam&> ExamTT::examsOf(ExamId a, ExamId b) {
    return {exams[a], exams[b]};
}

Two<Period&> ExamTT::periodsOf(PeriodId a, PeriodId b) {
    return {periods[a], periods[b]};
}

Two<Room const&> ExamTT::roomsOf(Two<RoomId> p) const {
    return {rooms[p.first], rooms[p.second]};
}

Two<Exam const&> ExamTT::examsOf(Two<ExamId> p) const {
    return {exams[p.first], exams[p.second]};
}

Two<Period const&> ExamTT::periodsOf(Two<PeriodId> p) const {
    return {periods[p.first], periods[p.second]};
}

Two<Room const&> ExamTT::roomsOf(RoomId a, RoomId b) const {
    return {rooms[a], rooms[b]};
}

Two<Exam const&> ExamTT::examsOf(ExamId a, ExamId b) const {
    return {exams[a], exams[b]};
}

Two<Period const&> ExamTT::periodsOf(PeriodId a, PeriodId b) const {
    return {periods[a], periods[b]};
}

ExamTT::ExamTT(char* instance_path) {
    InstanceParser parser(instance_path);
    parser.parse(*this);
}

double ExamTT::evaluateSolution(Solution & raw_solution) {
    double hardWeight = 1.0;

    ExamTTSolution & sol = (ExamTTSolution&) raw_solution;
    sol.computeCost(*this);
    sol.setSolutionValue(sol.costs.total(hardWeight));

    return sol.getSolutionValue();
}

/**
 * Solution components
 */

std::ostream& operator <<(std::ostream& out, HardCostComponents::Printer printer)
{
    auto const& hard = printer.self;

    return out << endl << "* Hard " << hard.sum() << endl
                << "     roomConstraint " << setw(10) << hard.roomConstraint << " exams" << endl
                << "  excessiveDuration " << setw(10) << hard.excessiveDuration << " exams" << endl
                << "   periodConstraint " << setw(10) << hard.periodConstraint() <<
                   " = AF(" << hard.periodConstraintAfter << ")"
                   " + EX(" << hard.periodConstraintExclusion << ")"
                   " + CO(" << hard.periodConstraintCoincidence << ")" << " exams" << endl
                << "  simultaneousExams " << setw(10) << hard.simultaneousExams << " students" << endl
                << "       overCapacity " << setw(10) << hard.overCapacity << " students" << endl;
}

std::ostream& operator <<(std::ostream& out, SoftCostComponents::Printer printer)
{
    auto const& soft = printer.self;
    auto const& instance = printer.instance;

    return out << endl << "* Soft " << soft.sum() << endl
               << "    periodsPenalty " << setw(10) << soft.periodsPenalty << endl
               << "      roomsPenalty " << setw(10) << soft.roomsPenalty << endl
               << "    twoExamsInARow " << setw(10) << soft.twoExamsInARow << " = "
                    << setw(5) << soft.twoExamsInARow / instance.institutionalWeightings.twoInARow << "x"
                    << setw(5) << instance.institutionalWeightings.twoInARow << " soft cost / student" << endl
               << "    twoExamsInADay " << setw(10) << soft.twoExamsInADay << " = "
                    << setw(5) << soft.twoExamsInADay / instance.institutionalWeightings.twoInADay << "x"
                    << setw(5) << instance.institutionalWeightings.twoInADay << " soft cost / student" << endl
               << "      periodSpread " << setw(10) << soft.periodSpread << " = "
                    << setw(5) << soft.periodSpread << "x    1 soft cost / student" << endl
               << "         frontload " << setw(10) << soft.frontload << " = "
                    << setw(5) << soft.frontload / instance.institutionalWeightings.frontload.penalty << "x"
                    << setw(5) << instance.institutionalWeightings.frontload.penalty << " soft cost / big exam too late" << endl
               << "     mixedDuration " << setw(10) << soft.mixedDuration << " = "
                    << setw(5) << soft.mixedDuration / instance.institutionalWeightings.nonMixedDurations << "x"
                    << setw(5) << instance.institutionalWeightings.nonMixedDurations << endl;
}

std::ostream& operator <<(std::ostream& out, CostComponents::Printer printer) {
    auto const& costs = printer.self;
    auto const& instance = printer.instance;

    return out << costs.hard.print(instance) << endl
               << costs.soft.print(instance) << endl;
}

/**
  STUFF
  */

std::ostream& operator <<(std::ostream& out, ExamTTSolution::Printer print) {
    print.self.printTo(print.instance, out);
    return out;
}

std::ostream& operator <<(std::ostream& out, ExamTTSolution::Writer obj) {
    obj.self.writeTo(out);
    return out;
}

void ExamTTSolution::writeTo(std::ostream& out) const {
    for(size_t i = 0; i < periods.size(); i++)
        out << periods[i] << ", " << rooms[i] << endl;
}

void ExamTTSolution::printTimelineTo(ExamTT const& instance, std::ostream& out) const {
    int E = instance.exams.size();
    int P = instance.periods.size();
    int R = instance.rooms.size();

    ExamTTSolution const& sol = *this;

    out << "* Timeline" << endl;
    int day = 0;
    for(PeriodId p = 0; p < P; p++) {
        if(p == 0 || instance.periods[p].date != instance.periods[p-1].date)
            out << string(5, '-') << "Day " << setw(3) << (day++) << string(5, '-') << endl;

        Period const& period = instance.periods[p];

        out << "   " << string(5, '_')
            << " " << PPrefix(p) << " "
            << (instance.periodInTheEnd(p) ? "[end]" : "")
            << "[duration " << instance.periods[p].duration << "]"
            << " " << string(5, '_')
            << endl;

        for(RoomId r = 0; r < R; r++) {
            Room const& room = instance.rooms[r];
            std::set<ExamId> exams = intersection(sol.examsByPeriods[p], sol.examsByRooms[r]);

            out << setw(5) << " " << RPrefix(r)
                << " : " << exams;

            std::vector<NumberOfStudents> sizes;
            for(ExamId i : exams)
                sizes.push_back(instance.exams[i].students.size());

            NumberOfStudents sumsizes = std::accumulate(sizes.begin(), sizes.end(), 0);

            std::set<Minutes> durations;
            std::vector<Minutes> excessiveDurations;

            for(ExamId i : exams) {
                Minutes duration = instance.exams[i].duration;
                durations.insert(duration);
                if(duration > instance.periods[p].duration)
                    excessiveDurations.push_back(duration);
            }

            out << " [capacity " << sumsizes << "/" << instance.rooms[r].capacity << "] " << sizes << (sumsizes > room.capacity ? "[over]" : "")
                << " [durations " << "[" << durations.size() << "]" << ":" << durations << (durations.size() > 1 ? "[mixed]" : "")
                << " [excessive " << "[" << excessiveDurations.size() << "]" << ":" << excessiveDurations << "]"
                << " [penalties " << (exams.size() * (instance.periods[p].penalty + instance.rooms[r].penalty)) << "]"
                    << " " << "[" << exams.size() << "x" << "(" << instance.periods[p].penalty << "+" << instance.rooms[r].penalty << ")" << "]"
                << endl;

            for(ExamId i : exams) {
                for(ExamId j : instance.afterOfExam(i))
                    if(periods[i] <= periods[j])
                        out << string(12, ' ') << "(" << i << "," << j << ")[after]" << endl;

                for(ExamId j : instance.exclusionOfExam(i))
                    if(periods[i] == periods[j])
                        out << string(12, ' ') << "(" << i << "," << j << ")[exclu]" << endl;

                for(ExamId j : instance.coincidenceOfExam(i))
                    if(periods[i] != periods[j])
                        out << string(12, ' ') << "(" << i << "," << j << ")[coinc]" << endl;
            }

            for(ExamId i : exams) {
                for(ExamId j : instance.hasStudentsInCommonOfExam(i)) {
                    Two<PeriodId> periodsId = periodsOf(i,j);
                    Two<Period const&> periodObjects = instance.periodsOf(periodsId);

                    Cost cost = instance.numberStudentsInCommon[i][j];
                    if(periodsId.first == periodsId.second) {
                        out << string(12, ' ') << "(" << i << "," << j << ")[simul:" << cost << "]" << endl;
                    } else {
                        NumberOfPeriods diff = abs(periodsId.first - periodsId.second);

                        if(periodObjects.first.date == periodObjects.second.date) {
                            if(diff == 1) {
                                out << string(12, ' ') << "(" << i << "," << j << ")[twoInARow:" << cost << "]" << endl;
                            } else {
                                out << string(12, ' ') << "(" << i << "," << j << ")[twoInADay:" << cost << "]" << endl;
                            }
                        }

                        if(diff <= instance.institutionalWeightings.periodSpread) {
                            out << string(12, ' ') << "(" << i << "," << j << ")[periodSpread:" << cost << "]" << endl;
                        }
                    }
                }
            }

            if(instance.periodInTheEnd(p))
                for(ExamId i : exams)
                    if(instance.isInFrontLoad[i])
                        out << string(12, ' ') << "(" << i << ")" << " [frontload]" << endl;

            if(exams.size() > 1)
                for(ExamId i : exams)
                    if(instance.examIsRoomExclusive[i])
                        out << string(12, ' ') << "(" << i << ")" << " [room exclusive]" << endl;

            out << endl;
        }

        out << endl;
    }
}

void ExamTTSolution::printTo(ExamTT const& instance, std::ostream& out) const {
    int E = instance.exams.size();
    int P = instance.periods.size();
    int R = instance.rooms.size();
    ExamTTSolution const& sol = *this;

    vector<tuple<RoomId, PeriodId, int>> byRoom(E);
    vector<tuple<PeriodId, RoomId, int>> byPeriod(E);

    for(int i = 0; i < E; i++) {
        byRoom[i]   = make_tuple(sol.rooms[i], sol.periods[i], i);
        byPeriod[i] = make_tuple(sol.periods[i], sol.rooms[i], i);
    }

    std::sort(byRoom.begin(), byRoom.end());
    std::sort(byPeriod.begin(), byPeriod.end());

    vector<bool> costExclusive(E, false);
    for(ExamId i : instance.examsRoomsExclusive) {
        auto a = assignement(i);

        for(ExamId j = 0; j < E; j++) {
            if(a == assignement(j) && i != j) {
                costExclusive[i] = true;
                break;
            }
        }
    }

    out << endl << "* Assignments " << endl
         << "R$|P$ = Room|Period soft cost" << endl
         << "rE = 1HC for roomExclusive" << endl
         << "eD = 1 HC for excessiveDuration)" << endl;

    for(ExamId i = 0; i < E; i++) {
        PeriodId p = sol.periods[i];
        RoomId r = sol.rooms[i];

        out << EPrefix(i)
             << PPrefix(sol.periods[i])
             << RPrefix(sol.rooms[i])
             << (instance.periods[p].penalty ? " " + to_string(instance.periods[p].penalty) + ".P$" : string(""))
             << (instance.rooms[r].penalty ? " " + to_string(instance.rooms[r].penalty) + ".R$" : string(""))
             << (costExclusive[i] ? " [rE]" : "")
             << (instance.exams[i].duration > instance.periods[p].duration ? " [eD]" : "")
             << endl;
    }

    out << endl << "* Sorted Assignements" << endl;
    for(ExamId i = 0; i < E; i++)
        out << "" << RPrefix(get<0>(byRoom[i])) << PPrefix(get<1>(byRoom[i])) << EPrefix(get<2>(byRoom[i]))
            << " |" << PPrefix(get<0>(byPeriod[i])) << RPrefix(get<1>(byPeriod[i])) << EPrefix(get<2>(byPeriod[i])) << endl;

    printTimelineTo(instance, out);

    out << endl << "* ExamsByPeriods " << endl;
    for(PeriodId i = 0; i < P; i++)
        out << PPrefix(i) << setw(7) << "[" + to_string(sol.examsByPeriods[i].size()) + "]"
             << " " << print(sol.examsByPeriods[i], nocomma) << endl;

    out << endl << "* ExamsByRooms" << endl;
    for(RoomId i = 0; i < R; i++)
        out << RPrefix(i) << setw(7) << "[" + to_string(sol.examsByRooms[i].size()) + "]"
             << " " << print(sol.examsByRooms[i], nocomma) << endl;

    out << endl << costs.print(instance) << endl;
}

int InstanceParser::readHeaderInt(std::string name){
    std::string line = readline();

    size_t i = 0;
    if(line.substr(i,1) != "[")
        throw invalid_argument("readHeader no opening [");
    i += 1;

    if(line.substr(i, name.size()) != name)
        throw invalid_argument("readHeader not correct name");
    i += name.size();

    if(line.substr(i, 1) != ":")
        throw invalid_argument("readHeader no :");
    i += 1;

    size_t a = i;
    while(i < line.size() && isDigit(line[i]))
        i++;

    if(a == i)
        throw invalid_argument("readHeader no number");

    int r = stoi(line.substr(a, i-a));
    if(line.substr(i) != "]")
        throw invalid_argument("readHeader no closing ]");

    return r;
}

void InstanceParser::readHeaderVoid(std::string name){
    std::string line = readline();

    int i = 0;
    if(line.substr(i,1) != "[")
        throw invalid_argument("readHeader no opening [");
    i += 1;

    if(line.substr(i, name.size()) != name)
        throw invalid_argument("readHeader not correct name");
    i += name.size();

    if(line.substr(i) != "]")
        throw invalid_argument("readHeader no closing ]");
}

std::string InstanceParser::readHeaderName() {
    std::string line = readline();

    if(line.substr(0,1) != "[")
        throw invalid_argument("readHeader no opening [");

    if(line.substr(line.size() - 1, 1) != "]")
        throw invalid_argument("readHeader no closing ]");

    return line.substr(0, line.size() - 1);
}

std::string InstanceParser::readline() {
    std::string line;
    readline(line);
    return line;
}

void InstanceParser::readline(std::string &line) {
    if(! std::getline(file, line))
        throw invalid_argument("no lines left");
    while(line.size() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' '))
        line.pop_back();
}

void InstanceParser::parse(ExamTT &i) {
    string line;

    i.exams.resize(readHeaderInt("Exams"));
    for(Exam& exam : i.exams) {
        readline(line);
        CommaSpaceSeparated tok(line);
        exam.duration = tok.nextInt();
        while(tok.hasNext()) {
            int student = tok.nextInt();
            bool inserted = exam.students.insert(student).second;
            if(! inserted)
                throw invalid_argument("one exam has a repeated student : " + to_string(student));
        }
    }

    i.periods.resize(readHeaderInt("Periods"));

    for(Period& period : i.periods){
        readline(line);
        CommaSpaceSeparated tok(line);
        period.date = makeDateFromITCFormat(tok.next());
        period.time = makeTimeFromITCFormat(tok.next());
        period.duration = tok.nextInt();
        period.penalty = tok.nextInt();
    }

    Minutes maxExams = max_element(i.exams.begin(), i.exams.end(), &Exam::compareDuration)->duration;
    Minutes maxPeriods = max_element(i.periods.begin(), i.periods.end(), &Period::compareDuration)->duration;

    if(maxExams > maxPeriods)
        throw invalid_argument("The longest exam is longer than the longest period " + to_string(maxExams) + " > " + to_string(maxPeriods));

    if(! i.periodsSorted())
        throw invalid_argument("Periods are not sorted");

    i.rooms.resize(readHeaderInt("Rooms"));

    for(Room& room : i.rooms) {
        readline(line);
        CommaSpaceSeparated tok(line);
        room.capacity = tok.nextInt();
        room.penalty = tok.nextInt();
    }

    enum class State {
        PeriodHardConstraints,
        RoomHardConstraints,
        InstitutionalWeightings,
        Null,
    };

    vector<string> nameOfState = {
        "PeriodHardConstraints",
        "RoomHardConstraints",
        "InstitutionalWeightings"
    };

    State state = State::Null;
    while(getline(file, line)) {
        while(line.size() && (line.back() == '\n' || line.back() == '\r' || line.back() == ' '))
            line.pop_back();
        // if line match \[.*\] try to change state
        if(line.substr(0,1) == "[" && line.substr(line.size() - 1, 1) == "]") {
            string name = line.substr(1, line.size() - 2);
            state = (State)(std::find(nameOfState.begin(), nameOfState.end(), name) - nameOfState.begin());
            // if not found state will be NoState
        } else switch(state) {
        case State::PeriodHardConstraints: {
            CommaSpaceSeparated tok(line);
            Two<int> pair;

            pair.first = tok.nextInt();
            string name = tok.next();
            pair.second = tok.nextInt();

            if(name == "AFTER")
                i.examsAfter.push_back(pair);
            else if(name == "EXCLUSION")
                i.examsExclusion.push_back(pair);
            else if(name == "EXAM_COINCIDENCE")
                i.examsCoincidence.push_back(pair);
            else
                throw invalid_argument("unknown PeriodHardConstraints : " + name);
        }
        break;

        case State::RoomHardConstraints: {
            CommaSpaceSeparated tok(line);
            int room = tok.nextInt();
            string name = tok.next();

            if(name == "ROOM_EXCLUSIVE")
                i.examsRoomsExclusive.insert(room);
            else
                throw invalid_argument("unknown PeriodHardConstraints : " + name);
        }
        break;

        case State::InstitutionalWeightings: {
            CommaSpaceSeparated tok(line);
            string name = tok.next();

            InstitutionalWeightings& w = i.institutionalWeightings;

            if(name == "TWOINAROW") {
                w.twoInARow = tok.nextInt();
            } else if(name == "TWOINADAY") {
                w.twoInADay = tok.nextInt();
            } else if(name == "PERIODSPREAD") {
                w.periodSpread = tok.nextInt();
            } else if(name == "NONMIXEDDURATIONS") {
                w.nonMixedDurations = tok.nextInt();
            } else if(name == "FRONTLOAD") {

                w.frontload.largestExams = tok.nextInt();
                w.frontload.time = tok.nextInt();
                w.frontload.penalty = tok.nextInt();

                if(w.frontload.largestExams >= (int)i.exams.size())
                    throw invalid_argument("Too many exams for the Frontload Constraint");
                if(w.frontload.time >= (int)i.periods.size())
                    throw invalid_argument("Too many periods for the Frontload Constraint");
                if(w.frontload.penalty < 0)
                    throw invalid_argument(name + " has a negative penalty");

            } else {
                throw invalid_argument("unknow InstitutionalWeightings : " + name);
            }
        }
        break;
        case State::Null:
            break;
        }
    }

    // check duplicate or reverse duplicate for AFTER, EXCLUSIVE, EXAM_COINCIDENCE

    for(size_t a = 0; a < i.examsAfter.size(); a++)
    for(size_t b = a + 1; b < i.examsAfter.size(); b++){
        Two<ExamId> p = i.examsAfter[a];
        Two<ExamId> q = i.examsAfter[b];
        if(p.first == q.first && p.second == q.second || p.second == q.first && p.first == q.second) {
            ostringstream oss;
            oss << "Exams after, duplicate #" << a << p << " and #" << b << q;
            throw invalid_argument(oss.str());
        }
    }

    for(std::vector<Two<ExamId>>* v : {&i.examsCoincidence, &i.examsExclusion}) {
        for(size_t a = 0; a < v->size(); a++)
        for(size_t b = a+1; b < v->size(); b++){
            Two<ExamId> p = (*v)[a];
            Two<ExamId> q = (*v)[b];
            if(p.first == q.first && p.second == q.second) {
                ostringstream oss;
                oss << "Exams coincidence or exclusion, duplicate #" << a << p << " and #" << b << q;
                throw invalid_argument(oss.str());
            }
        }
    }

    if(i.institutionalWeightings.periodSpread > i.periods.size())
        throw invalid_argument("Period Spread longer then the number of periods");
}

void InstanceParser::parse(ExamTT const& instance, ExamTTSolution &sol)
{
    string line;
    while(getline(file, line)) {
        CommaSpaceSeparated tok(line);
        tok.hasNext();
        sol.periods.push_back( tok.nextInt() );
        sol.rooms.push_back( tok.nextInt() );
    }

    if(instance.exams.size() != sol.periods.size())
        throw invalid_argument("Not correct number of assignments: " + to_string(sol.periods.size()));

    for(int period : sol.periods)
        if(! instance.correctPeriod(period))
            throw invalid_argument("Period " + to_string(period) + " does not exist");

    for(int room : sol.rooms)
        if(! instance.correctRoom(room))
            throw invalid_argument("Room " + to_string(room) + " does not exist");
}

// access
std::pair<PeriodId&, RoomId&> ExamTTSolution::assignement(ExamId exam) {
    return {periods[exam], rooms[exam]};
}

Two<std::pair<PeriodId&, RoomId&>> ExamTTSolution::assignementOf(Two<ExamId> exams) {
    return {assignement(exams.first), assignement(exams.second)};
}

Two<PeriodId &> ExamTTSolution::periodsOf(Two<ExamId> exams) {
    return {periods[exams.first], periods[exams.second]};
}

Two<RoomId &> ExamTTSolution::roomsOf(Two<ExamId> exams) {
    return {rooms[exams.first], rooms[exams.second]};
}

std::pair<std::pair<int &, int &>, std::pair<int &, int &> > ExamTTSolution::assignementOf(int a, int b) {
    return {assignement(a), assignement(b)};
}

std::pair<int &, int &> ExamTTSolution::periodsOf(int a, int b) {
    return {periods[a], periods[b]};
}

std::pair<int &, int &> ExamTTSolution::roomsOf(int a, int b) {
    return {rooms[a], rooms[b]};
}

// access const
std::pair<PeriodId const &, RoomId const &> ExamTTSolution::assignement(ExamId exam) const {
    return {periods[exam], rooms[exam]};
}

Two<std::pair<PeriodId const &, RoomId const &>> ExamTTSolution::assignementOf(Two<ExamId> exams) const {
    return {assignement(exams.first), assignement(exams.second)};
}

Two<PeriodId const &> ExamTTSolution::periodsOf(Two<ExamId> exams) const {
    return {periods[exams.first], periods[exams.second]};
}

Two<RoomId const &> ExamTTSolution::roomsOf(Two<ExamId> exams) const {
    return {rooms[exams.first], rooms[exams.second]};
}

std::pair<std::pair<int const &, int const &>, std::pair<int const &, int const &> > ExamTTSolution::assignementOf(int a, int b) const {
    return {assignement(a), assignement(b)};
}

std::pair<int const &, int const &> ExamTTSolution::periodsOf(int a, int b) const {
    return {periods[a], periods[b]};
}

std::pair<int const &, int const &> ExamTTSolution::roomsOf(int a, int b) const {
    return {rooms[a], rooms[b]};
}

void ExamTTSolution::initRandom(ExamTT const& instance, Random & r) {
    int E = instance.exams.size();
    int P = instance.periods.size();
    int R = instance.rooms.size();

    periods.resize(E);
    rooms.resize(E);

    for(ExamId i = 0; i < E; i++) {
        periods[i] = r.randrange(P);
        rooms[i] = r.randrange(R);
    }

    computeCost(instance);
    buildStructures(instance);
}

void ExamTTSolution::movePeriod(ExamTT const& instance, ExamId e, PeriodId p) {
    move(instance, e, p, rooms[e]);
}

void ExamTTSolution::moveRoom(ExamTT const& instance, ExamId e, RoomId r) {
    move(instance, e, periods[e], r);
}

void ExamTTSolution::move(ExamTT const& instance, ExamId e, PeriodId nextP, RoomId nextR) {
    PeriodId prevP = periods[e];
    RoomId prevR = rooms[e];

    updateMove(instance, e, nextP, nextR, this->costs);

    examsByPeriods[prevP].erase(examsByPeriodsIterators[e]);
    examsByPeriodsIterators[e] = examsByPeriods[nextP].insert(e).first;
    examsByRooms[prevR].erase(examsByRoomsIterators[e]);
    examsByRoomsIterators[e] = examsByRooms[nextR].insert(e).first;

    periods[e] = nextP;
    rooms[e] = nextR;
}

void ExamTTSolution::swap(ExamTT const& instance, ExamId e1, ExamId e2) {
    PeriodId p1 = periods[e1], p2 = periods[e2];
    RoomId r1 = rooms[e1], r2 = rooms[e2];
    move(instance, e1, p2, r2);
    move(instance, e2, p1, r1);
}

CostComponents ExamTTSolution::computeAndGetCost(ExamTT const& instance) const {
    CostComponents real;
    computeCost(instance, real);
    return real;
}

void ExamTTSolution::computeCost(ExamTT const& instance) {
    computeCost(instance, costs);
}

void ExamTTSolution::computeCost(ExamTT const& instance, CostComponents& costs) const {
    int E = instance.exams.size();
    int R = instance.rooms.size();
    int P = instance.periods.size();

    InstitutionalWeightings const& institutionalWeightings = instance.institutionalWeightings;
    FrontLoadParams const& frontload = institutionalWeightings.frontload;

    // period constraints : pair of exams

    costs.hard.periodConstraintAfter = 0;
    costs.hard.periodConstraintExclusion = 0;
    costs.hard.periodConstraintCoincidence = 0;

    for(Two<ExamId> exams : instance.examsAfter)
        if(periods[exams.first] <= periods[exams.second])
            costs.hard.periodConstraintAfter++;

    for(Two<ExamId> exams : instance.examsExclusion)
        if(periods[exams.first] == periods[exams.second])
            costs.hard.periodConstraintExclusion++;

    for(Two<ExamId> exams : instance.examsCoincidence)
        if(periods[exams.first] != periods[exams.second])
            costs.hard.periodConstraintCoincidence++;

    // pairs of exams with students in common

    costs.hard.simultaneousExams = 0;

    costs.soft.twoExamsInARow = 0;
    costs.soft.twoExamsInADay = 0;
    costs.soft.periodSpread = 0;

    for(ExamId i = 0; i < E; i++)
    for(ExamId j = i + 1; j < E; j++)
    if(instance.numberStudentsInCommon[i][j]) {
        Two<PeriodId> periodsId = periodsOf(i,j);
        Two<Period const&> periodObjects = instance.periodsOf(periodsId);

        Cost cost = instance.numberStudentsInCommon[i][j];
        if(periodsId.first == periodsId.second) {
            costs.hard.simultaneousExams += cost;
        } else {
            int diff = abs(periodsId.first - periodsId.second);

            if(periodObjects.first.date == periodObjects.second.date) {
                if(diff == 1) {
                    costs.soft.twoExamsInARow += cost;
                } else {
                    costs.soft.twoExamsInADay += cost;
                }
            }

            if(diff <= institutionalWeightings.periodSpread) {
                costs.soft.periodSpread += cost;
            }
        }
    }

    costs.soft.twoExamsInARow *= institutionalWeightings.twoInARow;
    costs.soft.twoExamsInADay *= institutionalWeightings.twoInADay;

    // frontload

    costs.soft.frontload = 0;

    // select the frontload.largestExams biggest exams
    // and give cost that the period is too late
    // aka period + frontload.time is out of range

    for(ExamId e : instance.frontLoadExams)
        if(instance.periodInTheEnd(periods[e]))
            costs.soft.frontload++;

    costs.soft.frontload *= frontload.penalty;

    // directly from assignement

    costs.hard.excessiveDuration = 0;

    costs.soft.periodsPenalty = 0;
    costs.soft.roomsPenalty = 0;

    MapVec<PeriodId, MapVec<RoomId, int>> numberStudentsInRoom(P, MapVec<RoomId, int>(R,0));
    MapVec<PeriodId, MapVec<RoomId, set<int>>> mixedDurations(P, MapVec<RoomId, set<int>>(R));

    PeriodId p;
    RoomId r;
    for(ExamId i = 0; i < E; i++){
        tie(p,r) = assignement(i);
        numberStudentsInRoom[p][r] += instance.exams[i].students.size();
        mixedDurations[p][r].insert(instance.exams[i].duration);

        costs.soft.periodsPenalty += instance.periods[p].penalty;
        costs.soft.roomsPenalty += instance.rooms[r].penalty;
        if(instance.exams[i].duration > instance.periods[p].duration)
            costs.hard.excessiveDuration++;
    }

    // linked between period and room : overCapacity and mixedDurations

    costs.hard.overCapacity = 0;

    costs.soft.mixedDuration = 0;

    for(PeriodId p = 0; p < P; p++)
    for(RoomId r = 0; r < R; r++) {
        int neededSeats = numberStudentsInRoom[p][r] - instance.rooms[r].capacity;
        if(neededSeats > 0)
            costs.hard.overCapacity += neededSeats;

        if(mixedDurations[p][r].size() > 1)
            costs.soft.mixedDuration += mixedDurations[p][r].size() - 1;
    }

    costs.soft.mixedDuration *= institutionalWeightings.nonMixedDurations;

    // exams that are room exclusive

    costs.hard.roomConstraint = 0;

    for(ExamId i : instance.examsRoomsExclusive) {
        auto a = assignement(i);

        for(ExamId j = 0; j < E; j++) {
            if(a == assignement(j) && i != j) {
                costs.hard.roomConstraint++;
                break;
            }
        }
    }
}

string lower(string s) {
    // first time we Want a copy...
    for(char& c : s)
        if('A' <= c && c <= 'Z')
            c += 'a' - 'A';
    return s;
}

void test() {
    string filename = "../simple1.exam"; // ../itc2007-exam-instances/exam_comp_set1.exam"; // args[0];
    ExamTT inst;
    ExamTT& instance = inst;
    InstanceParser parser(filename);

    if(! parser.file) {
        cerr << "NO FILE" << endl;
        return;
    }

    parser.parse(inst);
    inst.compute();
    inst.buildStructures();

    int E = inst.exams.size(),
        P = inst.periods.size(),
        R = inst.rooms.size();

    ofstream log("log");

    log << endl << "* Stats" << endl
         << "             Exams    " << inst.exams.size() << endl
         << "             Rooms    " << inst.rooms.size() << endl
         << "           Periods    " << inst.periods.size() << endl
         << "          Students    " << inst.students.size() << endl
         << "         Frontload    " << inst.institutionalWeightings.frontload << endl
         << " Conflict Density     " << inst.meta.ConflictDensity << endl
         << " conflictsComponents  " << inst.meta.conflictsComponents << " " << inst.meta.conflictsComponentsDensities << endl
         << " Period spread        " << inst.institutionalWeightings.periodSpread << "p" << endl
         << " Period contraints    " << inst.examsAfter.size() + inst.examsExclusion.size() + inst.examsCoincidence.size()
            << " = AF[" << inst.examsAfter.size() << "]"
            << " + EX[" << inst.examsExclusion.size() << "]"
            << " + CO[" << inst.examsCoincidence.size() << "]" << endl
         << " Room exclusive exams " << inst.examsRoomsExclusive.size() << endl;

    std::set<Minutes> allDurations;
    for(ExamId i = 0; i < E; i++)
        allDurations.insert(inst.exams[i].duration);
    log << endl << "* All mixed durations (" << allDurations.size() << ") : " << allDurations << endl;

    log << endl << "* Period constraints" << endl;

    log << endl << "** AF After" << endl;
    for(ExamId i = 0; i < E; i++)
        if(inst.afterOfExam(i).size())
            log << EPrefix(i) << " AF " << print(inst.afterOfExam(i), nocomma) << endl;

    log << endl << "** BE Before" << endl;
    for(ExamId i = 0; i < E; i++)
        if(inst.beforeOfExam(i).size())
            log << EPrefix(i) << " BE " << print(inst.beforeOfExam(i), nocomma) << endl;

    log << endl << "** EX Exclusion" << endl;
    for(ExamId i = 0; i < E; i++)
        if(inst.exclusionOfExam(i).size())
            log << EPrefix(i) << " EX " << print(inst.exclusionOfExam(i), nocomma) << endl;

    log << endl << "** CO Coincidence" << endl;
    for(ExamId i = 0; i < E; i++)
        if(inst.coincidenceOfExam(i).size())
            log << EPrefix(i) << " CO " << print(inst.coincidenceOfExam(i), nocomma) << endl;

    log << "* Exams | minutes | id | size | connected exams : student in commons : {exam1 exam2} [studentsWith1 studentsWith2]" << endl;
    for(ExamId i = 0; i < E; i++) {
        vector<KeyValue<int,int>> M;
        vector<int> M2;
        M.reserve(instance.hasStudentsInCommonOfExam(i).size());
        M.reserve(M.size());
        for(ExamId j : instance.hasStudentsInCommonOfExam(i)) {
            M.push_back({j, instance.numberStudentsInCommon[i][j]});
            M2.push_back(instance.numberStudentsInCommon[i][j]);
        }

        log << setw(5) << to_string(instance.exams[i].duration) << "m "
            << EPrefix(i)
            << setw(5) << "|" + to_string(instance.exams[i].students.size()) << "|" << " "
            << setw(3) << "[" + to_string(instance.hasStudentsInCommonOfExam(i).size()) << "] : "
            << print(M, nocomma | nobraces) << print(instance.hasStudentsInCommonOfExam(i), nocomma) << " " << print(M2, nocomma) << endl;
    }

    std::set<ExamId> setFrontloadExams(instance.frontLoadExams.begin(), instance.frontLoadExams.end());
    log << "* Frontload exams (" << instance.frontLoadExams.size() << ") : " << print(setFrontloadExams, nocomma) << " : ";
    for(ExamId i : instance.frontLoadExams)
        log << EPrefix(i) << "(" << inst.exams[i].students.size() << ") ";
    log << endl;

    log << "* Frontload periods (" << inst.institutionalWeightings.frontload.time << ") : "
        << PPrefix(P - inst.institutionalWeightings.frontload.time) << " - " << PPrefix(P-1) << endl;

    log << "* ROOM_EX Exams that must be room exclusive" << endl;
    for(ExamId e : inst.examsRoomsExclusive)
        log << e << endl;

    log << "* InstitutionalWeightings" << endl
         << inst.institutionalWeightings << endl;

    ExamTTSolution sol;

    Random r;
    r.seed(89);
    sol.initRandom(inst, r);

    /*
    InstanceParser parser2("../my-solutions/exam_comp_set1-seed-89-crand.sol"); // ("../my-solutions/art000-seed-89-crand.sol");
    parser2.parse(sol);
    */

    sol.computeCost(inst);
    sol.buildStructures(inst);

    log << endl << "* Solution" << endl;

    sol.printTo(inst, log);
    log.close();

    // hard after test
    // test on simple with 4 after 6
    vector<pair<int,int>> data;
    data.push_back({4,2});
    data.push_back({6,2});
    for(int i = 0; i < P; i++) {
        data.push_back({6, i});
        for(int j = 0; j < P; j++) {
            data.push_back({4, j});
            data.push_back({4, 2});
        }
        data.push_back({6, 2});
    }

    cout << sol.printer(inst);
    cout << "Going to move" << endl;
    if(data.size())
        cout << "Next " << EPrefix(data[0].first) << " to " << PPrefix(data[0].second) << endl;
    /*std::string s;
    getline(cin, s);*/

    for(size_t moveId = 0; moveId < data.size(); moveId++) {
        auto move = data[moveId];

        sol.movePeriod(inst, move.first, move.second);

        cout << sol.printer(inst);
        cout << "Moved " << EPrefix(move.first) << " to " << PPrefix(move.second) << endl;
        if(moveId < data.size() + 1)
            cout << "Next " << EPrefix(data[moveId+1].first) << " to " << PPrefix(data[moveId+1].second) << endl;

        if(sol.costs.hard.periodConstraintAfter != sol.computeAndGetCost(inst).hard.periodConstraintAfter) {
            cout << "Error " << " Real " << sol.computeAndGetCost(inst).hard.periodConstraintAfter << " vs diffed " << sol.costs.hard.periodConstraintAfter << endl;
            sol.computeCost(inst);
            break;
        }

        if(! sol.costs.exactlyEqual(sol.computeAndGetCost(inst))) {
            cout << "Error " << " Real " << sol.computeAndGetCost(inst).print(inst) << " vs diffed " << sol.costs.print(inst) << endl;
            sol.computeCost(inst);
            break;
        }

        /*std::string s;
        getline(cin, s);*/
    }

    bool isHelp = true;
    string str;
    for(;;) {
        if(isHelp)
            cout
             << "Change solution, C in front means compute and print difference in cost" << endl
             << "[C]Period <exam> <period>      : move <exam> to period <period>" << endl
             << "[C]Room <exam> <room>          : move <exam> to room <room>" << endl
             << "[C]Move <exam> <period> <room> : move <exam> to <period> <room>" << endl
             << "[C]Swap <exam1> <exam2>        : swap <exam> and <exam>" << endl
             << "Display                     : display current solution" << endl
             << "Write <filename>            : write current solution to <file>" << endl
             << "WriteVerbose <filename>     : write current solution verbosly <file>" << endl
             << "DCost                       : display current solution cost" << endl
             << "ReCalculate                 : recalculate cost of current solution" << endl
             << "Help                        : display help" << endl
            ;
        getline(cin, str);
        if(! cin)
            break;
        istringstream iss(str);
        string name;
        iss >> name;
        name = lower(name);

        bool isPeriod       = name == "p" || name == "period";
        bool isRoom         = name == "r" || name == "room";
        bool isMove         = name == "m" || name == "move";
        bool isSwap         = name == "s" || name == "swap";

        bool isPeriodC      = name == "cp" || name == "cperiod";
        bool isRoomC        = name == "cr" || name == "croom";
        bool isMoveC        = name == "cm" || name == "cmove";
        bool isSwapC        = name == "cs" || name == "cswap";

        bool isDisplay      = name == "d" || name == "disp" || name == "display";
        bool isWrite        = name == "w" || name == "write";
        bool isWriteVerbose = name == "wv" || name == "writeverbose";
        bool isDisplayCost  = name == "dc" || name == "dcost";
        bool isRecalculate  = name == "rc" || name == "recal" || name == "recalculate";
        isHelp              = name == "h" || name == "help" || name == "?";

        if(isPeriod || isRoom || isMove || isSwap || isPeriodC || isRoomC || isMoveC || isSwapC) {
            CostComponents realBefore = sol.computeAndGetCost(inst);

            CostComponents afterReal, diffAnnounced;

            if(isPeriod || isPeriodC) {
                int exam, period;
                iss >> exam >> period;
                if(iss && inst.correctExam(exam) && inst.correctPeriod(period)) {
                    if(isPeriod)
                        sol.move(inst, exam, period, sol.rooms[exam]);
                    else {
                        diffAnnounced = sol.differenceCostMove(inst, exam, period, sol.rooms[exam]);

                        int periodBefore = sol.periods[exam];
                        sol.movePeriod(inst, exam, period);
                        sol.computeCost(inst, afterReal);
                        sol.movePeriod(inst, exam, periodBefore);
                    }
                } else
                    cout << "Wront exam id or period id !";
            }
            else if(isRoom || isRoomC) {
                int exam, room;
                iss >> exam >> room;
                if(iss && inst.correctExam(exam) && inst.correctRoom(room)) {
                    if(isRoom)
                        sol.move(inst, exam, sol.periods[exam], room);
                    else {
                        diffAnnounced = sol.differenceCostMove(inst, exam, sol.periods[exam], room);

                        int roomBefore = sol.rooms[exam];
                        sol.moveRoom(inst, exam, room);
                        sol.computeCost(inst, afterReal);
                        sol.moveRoom(inst, exam, roomBefore);
                    }
                } else
                    cout << "Wrong id !" << endl;
            }
            else if(isMove || isMoveC) {
                int exam, period, room;
                iss >> exam >> period >> room;
                if(iss && inst.correctExam(exam) && inst.correctPeriod(period) && inst.correctRoom(room)) {
                    if(isMove)
                        sol.move(inst, exam, period, room);
                    else {
                        diffAnnounced = sol.differenceCostMove(inst, exam, period, room);

                        int roomBefore = sol.rooms[exam];
                        int periodBefore = sol.periods[exam];
                        sol.move(inst, exam, period, room);
                        sol.computeCost(inst, afterReal);
                        sol.move(inst, exam, periodBefore, roomBefore);
                    }
                } else
                    cout << "Wrong id !" << endl;
            }
            else if(isSwap || isSwapC){
                int exam1, exam2;
                iss >> exam1 >> exam2;
                if(iss && inst.correctExam(exam1) && inst.correctExam(exam2)) {
                    if(isSwap)
                        sol.swap(inst, exam1, exam2);
                    else {
                        diffAnnounced = sol.differenceCostSwap(inst, exam1, exam2);
                        sol.swap(inst, exam1, exam2);
                        sol.computeCost(inst, afterReal);
                        sol.swap(inst, exam1, exam2);
                    }
                } else
                    cout << "Wrong id !" << endl;
            }

            if(isPeriodC || isRoomC || isMoveC || isSwapC) {
                CostComponents realDiff = afterReal - realBefore;
                cout << "Announced difference" << diffAnnounced.print(instance) << endl;
                cout << "Real difference" << realDiff.print(instance) << endl;
                cout << (diffAnnounced.exactlyEqual(realDiff) ? "EQUAL" : "DIFFERENT") << endl;
            } else {
                CostComponents real = sol.computeAndGetCost(inst);
                if(! real.exactlyEqual(sol.costs)) {
                    cout << "Wrong cost difference : costs - real = " << (sol.costs - real).print(instance) << endl
                         << "costs - costsBefore = " << (sol.costs - realBefore).print(instance) << endl;
                } else
                    cout << "Diffed cost correctly calculated" << endl;

                cout << "Real changed cost : real - realBefore"
                     << (real - realBefore).print(instance) << endl;
            }
        } else if(isDisplay) {
            CostComponents recal = sol.computeAndGetCost(inst);
            std::cout << sol.printer(inst) << endl
                      << "* Real cost " << recal.print(instance) << endl
                      << "* Diff = real - costs" << (recal - sol.costs).print(instance) << endl;

        } else if(isWrite || isWriteVerbose) {
            string filename;
            iss >> filename;
            ofstream out(filename);
            if(isWrite)
                out << sol.writer();
            else if(isWriteVerbose) {
                CostComponents recal = sol.computeAndGetCost(inst);
                out << sol.printer(inst) << endl
                          << "* Real cost " << recal.print(instance) << endl
                          << "* Diff = real - costs" << (recal - sol.costs).print(instance) << endl;

            }
        } else if(isDisplayCost){
            CostComponents recal = sol.computeAndGetCost(inst);
            std::cout << sol.costs.print(instance) << endl
                      << "* Real cost " << recal.print(instance) << endl
                      << "* Diff = real - costs" << (recal - sol.costs).print(instance) << endl;
        } else if(isRecalculate){
            CostComponents before = sol.costs;
            sol.computeCost(inst);

            std::cout << "Real cost" << sol.costs.print(instance) << endl;
            if(sol.costs.exactlyEqual(before))
                cout << "No Error" << endl;
            else
                cout << "Error " << (sol.costs - before).print(instance);
        } else if(isHelp){

        } else {
            cout << "Unknown command '" << name << "'" << endl;
        }
    }

    /*

    // solve 1 after
    sol.move(106, 50, sol.rooms[106]);
    sol.print();

    CostComponents recalculate;
    sol.computeCost(recalculate);
    if(! recalculate.exactlyEqual(sol.costs))
        cout << "Wrong cost diffed !" << endl << recalculate.print(inst) << endl;

    // break 1 after
    sol.move(106, 14, sol.rooms[106]);
    sol.print();
    sol.computeCost(recalculate);
    if(! recalculate.exactlyEqual(sol.costs))
        cout << "Wrong cost diffed !" << endl << recalculate.print(inst) << endl;

    // solve 1 exclusive
    sol.move(175, 45, sol.rooms[175]);
    sol.print();
    sol.computeCost(recalculate);
    if(! recalculate.exactlyEqual(sol.costs))
        cout << "Wrong cost diffed !" << endl << recalculate.print(inst) << endl;

    // break 1 exclusive
    sol.move(175, 44, sol.rooms[175]);
    sol.print();
    sol.computeCost(recalculate);
    if(! recalculate.exactlyEqual(sol.costs))
        cout << "Wrong cost diffed !" << endl << recalculate.print(inst) << endl;

    for(int n = 0; n < 100; n++)
        cout << endl;
        */
}

/*
 * SOLUTION
 */

const void* ExamTTSolution::getRawData() const {
    return (const void*) this;
}

void ExamTTSolution::setRawData(const void *data) {
    ExamTTSolution const * other = (ExamTTSolution const *) data;
    // TODO
}

Solution* ExamTTSolution::clone() {
    ExamTTSolution* other = new ExamTTSolution(getSolutionValue());

    other->costs = costs;
    other->periods = periods;
    other->rooms = rooms;

    other->examsByPeriods.resize(periods.size());
    other->examsByRooms.resize(rooms.size());
    other->examsByPeriodsIterators.resize(periods.size());
    other->examsByRoomsIterators.resize(rooms.size());

    for(size_t e = 0; e < periods.size(); e++) {
        PeriodId p = periods[e];
        RoomId r = rooms[e];
        other->examsByPeriodsIterators[e] = other->examsByPeriods[p].insert(e).first;
        other->examsByRoomsIterators[e] = other->examsByRooms[r].insert(e).first;
    }

    return other;
}

/*
 * Methods
 */

void ExamTTSolution::buildStructures(InstanceRef instance)
{
   // structures are empty

   int P = instance.periods.size();
   int E = instance.exams.size();
   int R = instance.rooms.size();

   examsByPeriodsIterators.resize(E);
   examsByRoomsIterators.resize(E);
   examsByPeriods.resize(P);
   examsByRooms.resize(R);

   for(ExamId e = 0; e < E; e++) {
       PeriodId p = periods[e];
       RoomId r = rooms[e];
       examsByPeriodsIterators[e] = examsByPeriods[p].insert(e).first;
       examsByRoomsIterators[e] = examsByRooms[r].insert(e).first;
   }
}

CostComponents ExamTTSolution::differenceCostMove(Instance const& instance, ExamId ex, PeriodId nextP, RoomId nextR) const
{
    CostComponents diff = CostComponents::zero();
    updateMove(instance, ex, nextP, nextR, diff);
    return diff;
}

CostComponents ExamTTSolution::differenceCostSwap(Instance const& instance, ExamId e1, ExamId e2) const
{
    PeriodId p1 = periods[e1], p2 = periods[e2];
    RoomId r1 = rooms[e1], r2 = rooms[e2];

    ExamTTSolution& self = const_cast<ExamTTSolution&>(*this);

    CostComponents before = self.costs;
    self.move(instance, e1, p2, r2);
    self.move(instance, e2, p1, r1);
    CostComponents diff = self.costs - before;
    self.move(instance, e1, p1, r1);
    self.move(instance, e2, p2, r2);
    return diff;
}

void ExamTTSolution::updateMove(Instance const& instance, ExamId ex, PeriodId nextP, RoomId nextR, CostComponents& costs) const
{
    PeriodId prevP = periods[ex];
    RoomId prevR = rooms[ex];

    if(make_pair(prevP, prevR) == make_pair(nextP, nextR))
        return;

    cout << endl << "* Move E " << ex << endl;
    cout << endl << "** Periods" << endl;

    InstitutionalWeightings const& weightings = instance.institutionalWeightings;

    for(int n = 0; n < 2; n++) {
        PeriodId p = n == 0 ? prevP : nextP;
        cout << (n == 0 ? "From P" : "To   P") << setw(4) << p << " ";

        Instance::ExamsRelated const& related = instance.examsRelated[ex];

        cout << " AF=" << print(intersection(examsByPeriods[p], related.afters), nocomma);
        cout << " BE=" << print(intersection(examsByPeriods[p], related.befores), nocomma);
        cout << " CO=" << print(intersection(examsByPeriods[p], related.coincidences), nocomma);
        cout << " EX=" << print(intersection(examsByPeriods[p], related.exclusions), nocomma);

        set<ExamId> I = intersection(examsByPeriods[p], related.hasStudentsInCommon);

        std::vector<KeyValue<ExamId,int>> M;
        M.reserve(I.size());
        for(ExamId j : I)
            M.push_back({j, instance.numberStudentsInCommon[ex][j]});

        cout << " STU=" << print(M);

        cout << " -- " << print(examsByPeriods[p], nocomma) << endl;
    }

    Instance::ExamsRelated const& related = instance.examsRelated[ex];

    int s = nextP > prevP ? +1 :
            nextP < prevP ? -1 : 0;

    // if the period has changed
    if(s) {

        // simple period penalty
        costs.soft.periodsPenalty += instance.periods[nextP].penalty - instance.periods[prevP].penalty;

        /*
        costs.hard.periodConstraintCoincidence += countIntersection(examsByPeriods[prevP], instance.coincidenceOfExam[e])
                                                - countIntersection(examsByPeriods[nextP], instance.coincidenceOfExam[e]);

        costs.hard.periodConstraintExclusion -= countIntersection(examsByPeriods[prevP], instance.exclusionOfExam[e])
                                              - countIntersection(examsByPeriods[nextP], instance.exclusionOfExam[e]);
        */

        for(ExamId j : related.coincidences) {
            int d = periods[j] == prevP ? +1 :
                    periods[j] == nextP ? -1 : 0;
            costs.hard.periodConstraintCoincidence += d;
            costs.hard.periodConstraintExclusion -= d;
        }

        for(ExamId j : related.exclusions) {
            int d = periods[j] == prevP ? +1 :
                    periods[j] == nextP ? -1 : 0;
            costs.hard.periodConstraintCoincidence -= d;
            costs.hard.periodConstraintExclusion += d;
        }


        /*
        cout << "Passing through" << endl;

        for(int i = prevP + s; (nextP - i) * s >= 0; i += s) {
            cout << setw(4) << i << " ";

            cout << (s == 1 ? " -" : " +") << countIntersection(examsByPeriods[i], instance.afterOfExam[e])
                 << (s == 1 ? " +" : " -") << countIntersection(examsByPeriods[i], instance.beforeEqualOfExam[e]) << " ";

            set<int> I;
            makeIntersection(examsByPeriods[i], instance.afterOfExam[e], I);
            cout << "AF" << I;

            I.clear();
            makeIntersection(examsByPeriods[i], instance.beforeEqualOfExam[e], I);
            cout << "BE" << I;

            I.clear();
            makeIntersection(examsByPeriods[i], instance.coincidenceOfExam[e], I);
            cout << "CO" << I;

            I.clear();
            makeIntersection(examsByPeriods[i], instance.exclusionOfExam[e], I);
            cout << "EX" << I;

            cout << " -- " << examsByPeriods[i] << " ";
            cout << endl;

            costs.hard.periodConstraintAfter -= s * countIntersection(examsByPeriods[i], instance.afterOfExam[e]);
            costs.hard.periodConstraintAfter += s * countIntersection(examsByPeriods[i], instance.beforeEqualOfExam[e]);
        }
        */
        int m = min(prevP, nextP);
        int M = max(prevP, nextP);

        for(ExamId j : instance.afterOfExam(ex))
            if(m <= periods[j] && periods[j] < M)
                costs.hard.periodConstraintAfter -= s;

        for(ExamId j : instance.beforeOfExam(ex))
            if(m < periods[j] && periods[j] <= M)
                costs.hard.periodConstraintAfter += s;
    }

    // room penalty
    costs.soft.roomsPenalty += instance.rooms[nextR].penalty - instance.rooms[prevR].penalty;

    // about students in common

    for(ExamId j : related.hasStudentsInCommon) {
        Cost cost = instance.numberStudentsInCommon[ex][j];

        int relatedP = periods[j];
        int prevDiff = abs(relatedP - prevP);
        int nextDiff = abs(relatedP - nextP);

        Two<Period const&> before = {instance.periods[relatedP], instance.periods[prevP]};
        Two<Period const&> after  = {instance.periods[relatedP], instance.periods[nextP]};

        if(1) {
            if(nextP == relatedP)
                cout << EPrefix(j) << "[simul " << " @" << PPrefix(periods[j]) << " +" << instance.numberStudentsInCommon[ex][j] << endl;

            if(prevP == relatedP)
                cout << EPrefix(j) << "[simul " << " @" << PPrefix(periods[j]) << " -" << instance.numberStudentsInCommon[ex][j] << endl;


            if(before.first.date == before.second.date) {
                if(prevDiff == 1)
                    cout << EPrefix(j) << "[twoExamsInARow -" << cost * weightings.twoInARow << endl;
                else
                    cout << EPrefix(j) << "[twoExamsInADay -" << cost * weightings.twoInADay << endl;
            }

            if(after.first.date == after.second.date) {
                if(nextDiff == 1)
                    cout << EPrefix(j) << "[twoExamsInARow +" << cost * weightings.twoInARow << endl;
                else
                    cout << EPrefix(j) << "[twoExamsInADay +" << cost * weightings.twoInADay << endl;
            }
        }

        // remove the weight from before
        if(prevP == relatedP)
            costs.hard.simultaneousExams -= cost;

        if(before.first.date == before.second.date) {
            if(prevDiff == 1)
                costs.soft.twoExamsInARow -= cost * weightings.twoInARow;
            else
                costs.soft.twoExamsInADay -= cost * weightings.twoInADay;
        }

        if(prevDiff <= weightings.periodSpread)
            costs.soft.periodSpread -= cost;

        // add the new
        if(nextP == relatedP)
            costs.hard.simultaneousExams += cost;

        if(after.first.date == after.second.date) {
            if(nextDiff == 1)
                costs.soft.twoExamsInARow += cost * weightings.twoInARow;
            else
                costs.soft.twoExamsInADay += cost * weightings.twoInADay;
        }

        if(nextDiff <= weightings.periodSpread)
            costs.soft.periodSpread += cost;
    }

    // frontload

    if(instance.isInFrontLoad[ex]) {
        if(instance.periodInTheEnd(prevP))
            costs.soft.frontload -= weightings.frontload.penalty;
        if(instance.periodInTheEnd(nextP))
            costs.soft.frontload += weightings.frontload.penalty;
    }

    // excessiveDuration
    if(instance.exams[ex].duration > instance.periods[prevP].duration)
        costs.hard.excessiveDuration--;
    if(instance.exams[ex].duration > instance.periods[nextP].duration)
        costs.hard.excessiveDuration++;

    // comparing assignement before and after

    set<ExamId> leaving; // exams that I am leaving
    makeIntersection(examsByRooms[prevR], examsByPeriods[prevP], leaving);

    set<ExamId> meeting; // exams that I will meet
    makeIntersection(examsByRooms[nextR], examsByPeriods[nextP], meeting);

    leaving.erase(leaving.find(ex)); // I am not included in the leaving

    // room exclusive

    if(leaving.size() == 1 && instance.examIsRoomExclusive[*leaving.begin()])
        costs.hard.roomConstraint--;
    if(leaving.size() > 0 && instance.examIsRoomExclusive[ex])
        costs.hard.roomConstraint--;

    if(meeting.size() == 1 && instance.examIsRoomExclusive[*meeting.begin()])
        costs.hard.roomConstraint++;
    if(meeting.size() > 0 && instance.examIsRoomExclusive[ex])
        costs.hard.roomConstraint++;

    // overcapacity
    int mySize = instance.exams[ex].students.size();

    int leavingSum = 0; // = sum(instance.exams[j].students.size() for j in leaving)
    int meetingSum = 0; // = sum(instance.exams[j].students.size() for j in meeting)

    for(ExamId j : leaving)
        leavingSum += instance.exams[j].students.size();

    for(ExamId j : meeting)
        meetingSum += instance.exams[j].students.size();

    int capLeaving = instance.rooms[prevR].capacity;
    int capMeeting = instance.rooms[nextR].capacity;

    if(leavingSum <= capLeaving && capLeaving < leavingSum + mySize)
        costs.hard.overCapacity -= leavingSum + mySize - capLeaving;

    if(meetingSum <= capMeeting && capMeeting < meetingSum + mySize)
        costs.hard.overCapacity += meetingSum + mySize - capMeeting;

    if(1) {
        cout << "Leaving " << leaving << " Meeting " << meeting << endl;
        if(leavingSum <= capLeaving && capLeaving < leavingSum + mySize)
            cout << "[overCapacity -" << leavingSum + mySize - capLeaving << endl;

        if(meetingSum <= capMeeting && capMeeting < meetingSum + mySize)
            cout << "[overCapacity +" << meetingSum + mySize - capMeeting << endl;
    }

    // mixedDurations
    int myDuration = instance.exams[ex].duration;
    std::set<Minutes> durationsLeaving; // = set(instance.exams[j].duration for j in leaving)
    std::set<Minutes> durationsMeeting; // = set(instance.exams[j].duration for j in meeting)

    for(ExamId j : leaving)
        durationsLeaving.insert(instance.exams[j].duration);
    for(ExamId j : meeting)
        durationsMeeting.insert(instance.exams[j].duration);

    // counters

    int NLeavingWithoutMe = durationsLeaving.size();
    int NMeetingWithoutMe = durationsMeeting.size();
    int NLeavingWithMe = durationsLeaving.size() + (durationsLeaving.count(myDuration) ? 0 : 1);
    int NMeetingWithMe = durationsMeeting.size() + (durationsMeeting.count(myDuration) ? 0 : 1);

    if(NLeavingWithoutMe > 0)
        costs.soft.mixedDuration += (NLeavingWithoutMe - NLeavingWithMe) * weightings.nonMixedDurations;
    if(NMeetingWithoutMe > 0)
        costs.soft.mixedDuration += (NMeetingWithMe - NMeetingWithoutMe) * weightings.nonMixedDurations;

    /*
    // no counters
    if(durationsLeaving.size() && durationsLeaving.count(myDuration))
        costs.soft.mixedDuration -= weightings.nonMixedDurations;
    if(durationsMeeting.size() && durationsMeeting.count(myDuration))
        costs.soft.mixedDuration += weightings.nonMixedDurations;
    */
}

int seconds(const Time &time) {
    return std::get<2>(time) + std::get<1>(time) * 60 + std::get<0>(time) * 3600;
}

ostream &operator <<(ostream &out, const InstitutionalWeightings &w){
    return out << "TWO IN A ROW       " << w.twoInARow << "$" << endl
               << "TWO IN A DAY       " << w.twoInADay << "$" << endl
               << "NON MIXED DURATION " << w.nonMixedDurations << "$" << endl
               << "PERIOD SPREAD      " << w.periodSpread << "P" << endl
               << "FRONTLOAD          " << w.frontload << endl;
}

ostream &operator <<(ostream &out, const FrontLoadParams &w){
    return out << w.largestExams << "e " << w.time << "p " << w.penalty << "$";
}

}
}
