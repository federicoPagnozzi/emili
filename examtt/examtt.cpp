#include "examtt.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <functional>
#include <tuple>
#include <memory>
#include <set>
#include <unordered_set>
#include <cassert>

using namespace std;

namespace emili
{
namespace ExamTT
{

int ExamTTSolution::numberOfClones = 0, ExamTTSolution::numberOfTotalCompute = 0;

template<typename Type, unsigned N, unsigned Last>
struct tuple_operations {

    static void print(std::ostream& out, const Type& value) {
        out << std::get<N>(value) << ", ";
        tuple_operations<Type, N + 1, Last>::print(out, value);
    }
};

template<typename Type, unsigned N>
struct tuple_operations<Type, N, N> {

    static void print(std::ostream& out, const Type& value) {
        out << std::get<N>(value);
    }
};

template<typename... Types>
std::ostream& operator<<(std::ostream& out, const std::tuple<Types...>& value) {
    out << "(";
    tuple_operations<std::tuple<Types...>, 0, sizeof...(Types) - 1>::print(out, value);
    out << ")";
    return out;
}

template <typename U, typename V>
std::ostream& operator <<(std::ostream& out, const std::pair<U,V>& value) {
    return out << "(" << value.first << ", " << value.second << ")";
}

/*
op = '-'
def make(N):
    letters = [chr(ord('A') + i) for i in range(N)]
    return '''template <{typenames}>\nstd::tuple<{types}> operator{op}(std::tuple<{types}> const& a, std::tuple<{types}> const& b) {{\n    return std::make_tuple({body});\n}}'''.format(
        typenames=', '.join(map('typename {}'.format, letters)),
        types=', '.join(letters),
        op=op,
        body=', '.join(map(('std::get<{0}>(a) ' + op + ' std::get<{0}>(b)').format, map(str,range(N))))
    )
print('\n'.join(map(make, range(10))))
*/

template <typename T, typename U>
struct KeyValue {
    T key;
    U value;
};

template <typename T, typename U>
bool operator <(KeyValue<T,U> const& a, KeyValue<T,U> const& b) {
    return make_pair(a.key, a.value) < make_pair(b.key, b.value);
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

template <typename U>
struct VectorPrinter {
    U const& self;
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

template <typename U>
VectorPrinter<U> print(U const& t) {
    return {t, 0};
}

template <typename U>
VectorPrinter<U> print(U const& t, unsigned int flags) {
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

const Prefix EPrefix{"e", 5};
const Prefix PPrefix{"p", 4};
const Prefix RPrefix{"r", 4};

bool ExamTT::periodInTheEnd(PeriodId p) const {
    return ! correctPeriod(p + institutionalWeightings.frontload.time);
}

int ExamTT::sizeOfExam(ExamId e) const {
    return exams[e].students.size();
}

namespace {
void removeDuplicate(std::vector<int> & vec) {
    std::set<int> s(vec.begin(), vec.end());
    vec.assign(s.begin(), s.end());
}
}

void ExamTT::ExamsRelated::removeDuplicates() {
    removeDuplicate(coincidences);
    removeDuplicate(exclusions);
    removeDuplicate(afters);
    removeDuplicate(befores);
    removeDuplicate(hasStudentsInCommon);
}

void ExamTT::buildStructures()  {
    int E = exams.size();

    examsRelated.resize(E);

    for(Two<ExamId> p : examsCoincidence) {
        coincidenceOfExam(p.first).push_back(p.second);
        coincidenceOfExam(p.second).push_back(p.first);
    }

    for(Two<ExamId> p : examsExclusion) {
        exclusionOfExam(p.first).push_back(p.second);
        exclusionOfExam(p.second).push_back(p.first);
    }

    for(Two<ExamId> p : examsAfter) {
        afterOfExam(p.first).push_back(p.second);
        beforeOfExam(p.second).push_back(p.first);
    }

    for(ExamId i = 0; i < E; i++)
    for(ExamId j = i + 1; j < E; j++) {
        if(numberStudentsInCommon[i][j]) {
            hasStudentsInCommonOfExam(i).push_back(j);
            hasStudentsInCommonOfExam(j).push_back(i);
        }
    }

    for(ExamId e = 0; e < E; e++)
        examsRelated[e].removeDuplicates();

    // frontLoadExams
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

    // daysToPeriod
    daysToPeriod.push_back({0,-1});

    for(int p = 0; p < P(); p++) {
        if(periods[p].dateid != periods[daysToPeriod.back().first].dateid) {
            daysToPeriod.back().second = p;
            daysToPeriod.push_back({p,-1});
        }
    }

    daysToPeriod.back().second = P();
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

int ExamTT::numberOfExamsOfStudent(int student) {
    int s = 0;
    for(Exam& e : exams)
        s += e.students.count(student);
    return s;
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

    // will be featured based, now it's just a constant
    hardWeight = 1000.0;
}

bool ExamTT::periodsSorted() const {
    return is_sorted(periods.begin(), periods.end(), &Period::compareDateTime);
}

bool ExamTT::allPeriodsConsecutives() const {
    for(size_t i = 0; i < periods.size() - 1; i++)
        if(periods[i].dateid == periods[i+1].dateid)
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
    cout << "path" << this << endl;
    InstanceParser parser(instance_path);
    if(parser.file)
        parser.parse(*this);
    else
        throw invalid_argument(string(instance_path) + " does not exist !");
}

double ExamTT::evaluateSolution(Solution & raw_solution) {
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
                   " + CO(" << hard.periodConstraintCoincidence << ")" << " pairs of exams" << endl
                << "  simultaneousExams " << setw(10) << hard.simultaneousExams << " students" << endl
                << "       overCapacity " << setw(10) << hard.overCapacity << " students" << endl;
}

std::ostream& operator <<(std::ostream& out, SoftCostComponents::Printer printer)
{
    auto const& soft = printer.self;
    auto const& instance = printer.instance;

    auto div0 = [](int a, int b) -> std::string {
        return b == 0 ? to_string(a) + "/0" : to_string(a / b);
    };

    return out << endl << "* Soft " << soft.sum() << endl
               << "    periodsPenalty " << setw(10) << soft.periodsPenalty << endl
               << "      roomsPenalty " << setw(10) << soft.roomsPenalty << endl
               << "    twoExamsInARow " << setw(10) << soft.twoExamsInARow << " = "
                    << setw(5) << div0(soft.twoExamsInARow, instance.institutionalWeightings.twoInARow) << "x"
                    << setw(5) << instance.institutionalWeightings.twoInARow << " soft cost / student" << endl
               << "    twoExamsInADay " << setw(10) << soft.twoExamsInADay << " = "
                    << setw(5) << div0(soft.twoExamsInADay, instance.institutionalWeightings.twoInADay) << "x"
                    << setw(5) << instance.institutionalWeightings.twoInADay << " soft cost / student" << endl
               << "      periodSpread " << setw(10) << soft.periodSpread << " = "
                    << setw(5) << soft.periodSpread << "x    1 soft cost / student" << endl
               << "         frontload " << setw(10) << soft.frontload << " = "
                    << setw(5) << div0(soft.frontload, instance.institutionalWeightings.frontload.penalty) << "x"
                    << setw(5) << instance.institutionalWeightings.frontload.penalty << " soft cost / big exam too late" << endl
               << "     mixedDuration " << setw(10) << soft.mixedDuration << " = "
                    << setw(5) << div0(soft.mixedDuration, instance.institutionalWeightings.nonMixedDurations) << "x"
                    << setw(5) << instance.institutionalWeightings.nonMixedDurations << " soft cost / period / room " << endl;
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
    if(! hasStructures)
        const_cast<ExamTTSolution*>(this)->buildStructures(instance);

    int E = instance.exams.size();
    int P = instance.periods.size();
    int R = instance.rooms.size();

    ExamTTSolution const& sol = *this;

    out << "* Timeline" << endl;
    int day = 0;
    for(PeriodId p = 0; p < P; p++) {
        if(p == 0 || instance.periods[p].dateid != instance.periods[p-1].dateid)
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
            // std::set<ExamId> exams = intersection(sol.examsByPeriods[p], sol.examsByRooms[r]);
            std::list<ExamId> const& exams = examsByPeriodRoom[p][r];

            out << setw(5) << " " << RPrefix(r)
                << " : " << print(exams);

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

            for(ExamId i : exams) if(periods[i] != -1) {
                for(ExamId j : instance.afterOfExam(i)) if(periods[j] != -1)
                    if(periods[i] <= periods[j])
                        out << string(12, ' ') << "(" << i << "," << j << ")[after]" << endl;

                for(ExamId j : instance.exclusionOfExam(i)) if(periods[j] != -1)
                    if(periods[i] == periods[j])
                        out << string(12, ' ') << "(" << i << "," << j << ")[exclu]" << endl;

                for(ExamId j : instance.coincidenceOfExam(i)) if(periods[j] != -1)
                    if(periods[i] != periods[j])
                        out << string(12, ' ') << "(" << i << "," << j << ")[coinc]" << endl;
            }

            for(ExamId i : exams) if(periods[i] != -1) {
                for(ExamId j : instance.hasStudentsInCommonOfExam(i)) if(periods[j] != -1) {
                    Two<PeriodId> periodsId = periodsOf(i,j);
                    Two<Period const&> periodObjects = instance.periodsOf(periodsId);

                    Cost cost = instance.numberStudentsInCommon[i][j];
                    if(periodsId.first == periodsId.second) {
                        out << string(12, ' ') << "(" << i << "," << j << ")[simul:" << cost << "]" << endl;
                    } else {
                        NumberOfPeriods diff = abs(periodsId.first - periodsId.second);

                        if(periodObjects.first.dateid == periodObjects.second.dateid) {
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
                for(ExamId i : exams) if(isAssigned(i))
                    if(instance.isInFrontLoad[i])
                        out << string(12, ' ') << "(" << i << ")" << " [frontload]" << endl;

            if(exams.size() > 1)
                for(ExamId i : exams) if(isAssigned(i))
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
             << (p != -1 && instance.periods[p].penalty ? " " + to_string(instance.periods[p].penalty) + ".P$" : string(""))
             << (r != -1 && instance.rooms[r].penalty ? " " + to_string(instance.rooms[r].penalty) + ".R$" : string(""))
             << (costExclusive[i] ? " [rE]" : "")
             << (p != -1 && instance.exams[i].duration > instance.periods[p].duration ? " [eD]" : "")
             << endl;
    }

    out << endl << "* Sorted Assignements" << endl;
    for(ExamId i = 0; i < E; i++)
        out << "" << RPrefix(get<0>(byRoom[i])) << PPrefix(get<1>(byRoom[i])) << EPrefix(get<2>(byRoom[i]))
            << " |" << PPrefix(get<0>(byPeriod[i])) << RPrefix(get<1>(byPeriod[i])) << EPrefix(get<2>(byPeriod[i])) << endl;

    printTimelineTo(instance, out);

    /*
    out << endl << "* ExamsByPeriods " << endl;
    for(PeriodId i = 0; i < P; i++)
        out << PPrefix(i) << setw(7) << "[" + to_string(sol.examsByPeriods[i].size()) + "]"
             << " " << print(sol.examsByPeriods[i], nocomma) << endl;

    out << endl << "* ExamsByRooms" << endl;
    for(RoomId i = 0; i < R; i++)
        out << RPrefix(i) << setw(7) << "[" + to_string(sol.examsByRooms[i].size()) + "]"
             << " " << print(sol.examsByRooms[i], nocomma) << endl;
    */

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

namespace {
template <typename T>
void mapToSortedInt(std::map<T, int> & M) {
    int x = 0;
    for(auto& p : M)
        p.second = x++;
}
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

    for(Exam& e : i.exams)
        i.colorsOfMinute.insert(make_pair(e.duration, i.colorsOfMinute.size()));

    mapToSortedInt(i.colorsOfMinute);

    for(Exam& e : i.exams)
        e.durationColor = i.colorsOfMinute[e.duration];

    i.periods.resize(readHeaderInt("Periods"));

    for(Period& period : i.periods){
        readline(line);
        CommaSpaceSeparated tok(line);
        period.dateRaw = makeDateFromITCFormat(tok.next());
        period.time = makeTimeFromITCFormat(tok.next());
        period.duration = tok.nextInt();
        period.penalty = tok.nextInt();
    }

    std::map<Date, Day> dateToDateId;
    for(Period& period : i.periods)
        dateToDateId.insert(make_pair(period.dateRaw, dateToDateId.size()));

    mapToSortedInt(dateToDateId);

    for(Period& period : i.periods)
        period.dateid = dateToDateId[period.dateRaw];

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
                if(w.frontload.time >= (int)i.periods.size()) {
                    cout << "Warning: Too many periods for the Frontload Constraint" << endl;
                    w.frontload.time = i.periods.size();
                }
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

    constexpr bool STOP_WHEN_DUPLICATE = false;

    // AFTER
    {
        std::vector<Two<ExamId>> res;
        for(size_t a = 0; a < i.examsAfter.size(); a++) {
            Two<ExamId> p = i.examsAfter[a];
            bool hasDuplicate = false;
            for(size_t b = a + 1; b < i.examsAfter.size(); b++) {
                Two<ExamId> q = i.examsAfter[b];
                if(p.first == q.first && p.second == q.second) {
                    hasDuplicate = true;
                    break;
                }
            }
            if(! hasDuplicate)
                res.push_back(p);
        }

        if(STOP_WHEN_DUPLICATE) {
            if(i.examsAfter.size() != res.size()) {
                ostringstream oss;
                oss << "Exams after, exist a duplicate";
                throw invalid_argument(oss.str());
            }
        } else {
            i.examsAfter.swap(res);
        }
    }

    // EXCLUSIVE, EXAM_COINCIDENCE remove duplicate
    {
        for(std::vector<Two<ExamId>>* v : {&i.examsCoincidence, &i.examsExclusion}) {
            std::vector<Two<ExamId>> res;
            for(size_t a = 0; a < v->size(); a++) {
                Two<ExamId> p = (*v)[a];
                bool hasDuplicate = false;
                for(size_t b = a+1; b < v->size(); b++){
                    Two<ExamId> q = (*v)[b];
                    if(p.first == q.first && p.second == q.second || p.second == q.first && p.first == q.second) {
                        hasDuplicate = true;
                        break;
                    }
                }
                if(! hasDuplicate)
                    res.push_back(p);
            }
            if(STOP_WHEN_DUPLICATE) {
                if(v->size() != res.size()) {
                    ostringstream oss;
                    oss << "Exams coincidence or exclusion, exists duplicate";
                    throw invalid_argument(oss.str());
                }
            } else {
                (*v).swap(res);
            }
        }
    }

    // remove aberations : x COINCIDENCE x
    {
        std::vector<Two<ExamId>> res;
        for(auto p : i.examsCoincidence)
            if(p.first != p.second)
                res.push_back(p);
        res.swap(i.examsCoincidence);
    }

    // test for impossibility x EXCLUSION x, x AFTER x
    {
        for(auto p : i.examsExclusion)
            if(p.first == p.second) {
                ostringstream oss;
                oss << "Exam exclusion " << p.first << " EXCLUSION " << p.second;
                throw invalid_argument(oss.str());
            }

        for(auto p : i.examsAfter)
            if(p.first == p.second) {
                ostringstream oss;
                oss << "Exam after " << p.first << " AFTER " << p.second;
                throw invalid_argument(oss.str());
            }
    }

    /*
    for(std::vector<Two<ExamId>>* v : {&i.examsAfter, &i.examsCoincidence, &i.examsExclusion}) {
        cout << "--" << endl;
        for(auto& x : *v)
            cout << x << endl;
    }
    */

    if(i.institutionalWeightings.periodSpread > i.periods.size()) {
        cout << "Warning: " << ("Period Spread longer than the number of periods") << endl;
        i.institutionalWeightings.periodSpread = i.periods.size();
    }

    i.compute();
    i.buildStructures();
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

int ExamTTSolution::sizeOfPartialSolution() const {
    return (int)rooms.size() - (int)unAssignedExamList.size();
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

void ExamTTSolution::initFromAssign(ExamTT const& instance, std::vector<std::pair<int,int>> assign) {
    int E = instance.exams.size();
    int P = instance.periods.size();
    int R = instance.rooms.size();
    periods.resize(E);
    rooms.resize(E);

    if(assign.size() != E)
        throw invalid_argument("assign.size != E");

    for(ExamId i = 0; i < E; i++)
        assignement(i) = assign[i];

    computeCost(instance);
    buildStructures(instance);
}

void ExamTTSolution::initFromZeroPeriodAndZeroRoom(ExamTT const& instance) {
    const int E = instance.E();

    std::vector<std::pair<int,int>> assign(E, std::pair<int,int>(0,0));
    initFromAssign(instance, assign);
}

void ExamTTSolution::initFromPeriodsAndZeroRoom(ExamTT const& instance, const std::vector<int> &initPeriods) {
    const int E = instance.E();

    std::vector<std::pair<int,int>> assign(E);
    for(int i = 0; i < E; i++)
        assign[i] = {initPeriods[i], 0};
    initFromAssign(instance, assign);
}

void ExamTTSolution::initRandom(InstanceRef instance) {
    Random r;
    initRandom(instance, r);
}

void ExamTTSolution::initUnassigned(ExamTTSolution::InstanceRef instance) {
    const int E = instance.E(), P = instance.P(), R = instance.R();
    periods.assign(E, -1);
    rooms.assign(E, -1);

    for(int i = 0; i < E; i++)
        unAssignedExamList.push_back(i);

    computeCost(instance);
    // costs = CostComponents::zero();
    // refreshSolutionValue(instance);
    buildStructures(instance);
}

void ExamTTSolution::initRandom(ExamTT const& instance, Random & r) {
    const int E = instance.E(), P = instance.P(), R = instance.R();
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
    if(! hasStructures)
        buildStructures(instance);

    if(USE_DELTA)
        updateMove(instance, e, nextP, nextR, this->costs);

    applyMove(instance, e, nextP, nextR);

    if(! USE_DELTA)
        computeCost(instance);

    refreshSolutionValue(instance);
}

void ExamTTSolution::refreshSolutionValue(InstanceRef instance) {
    setSolutionValue(costs.total(instance.hardWeight));
}

void ExamTTSolution::applyMove(InstanceRef instance, ExamId e, PeriodId nextP, RoomId nextR) {
    PeriodId prevP = periods[e];
    RoomId prevR = rooms[e];

    if(USE_COLOR_STRUCTURE) {
        durationColorUsed[prevP][prevR][instance.exams[e].durationColor]--;
        durationColorUsed[nextP][nextR][instance.exams[e].durationColor]++;
    }

    examsByPeriodRoom[nextP][nextR].splice(
        examsByPeriodRoom[nextP][nextR].end(),
        examsByPeriodRoom[prevP][prevR],
        examsByPeriodRoomIterators[e]
    );
    // iterators are still valid !

    periods[e] = nextP;
    rooms[e] = nextR;
}

void ExamTTSolution::removeExam(ExamTT const& instance, ExamId e) {
    if(! USE_DELTA)
        throw std::invalid_argument("! USE_DELTA");

    int prevP = periods[e];
    int prevR = rooms[e];

    if(1) {
        updateRemove(instance, e);
        applyRemoveExam(instance, e, prevP, prevR);
    } else {
        applyRemoveExam(instance, e, prevP, prevR);
        updateRemoveOutOfStructures(instance, e, prevP, prevR);
    }

    periods[e] = -1;
    rooms[e] = -1;

    refreshSolutionValue(instance);
}

inline void ExamTTSolution::addExam(ExamTTSolution::InstanceRef instance, ExamId e, Assignement p) {
    addExam(instance, e, p.first, p.second);
}

void ExamTTSolution::applyRemoveExam(InstanceRef instance, ExamId e, PeriodId prevP, RoomId prevR) {
    if(USE_COLOR_STRUCTURE)
        durationColorUsed[prevP][prevR][instance.exams[e].durationColor]--;

    unAssignedExamList.splice(
        unAssignedExamList.end(),
        examsByPeriodRoom[prevP][prevR],
        examsByPeriodRoomIterators[e]
    );
}

void ExamTTSolution::addExam(InstanceRef instance, ExamId e, PeriodId p, RoomId r) {
    if(! USE_DELTA)
        throw std::invalid_argument("! USE_DELTA");

    updateAdd(instance, e, p, r);
    applyAddExam(instance, e, p, r);

    refreshSolutionValue(instance);
}

void ExamTTSolution::applyAddExam(InstanceRef instance, ExamId e, PeriodId p, RoomId r) {
    if(USE_COLOR_STRUCTURE)
        durationColorUsed[p][r][instance.exams[e].durationColor]++;

    examsByPeriodRoom[p][r].splice(
        examsByPeriodRoom[p][r].end(),
        unAssignedExamList,
        examsByPeriodRoomIterators[e]
    );

    periods[e] = p;
    rooms[e] = r;
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
    refreshSolutionValue(instance);
}

void ExamTTSolution::computeCost(ExamTT const& instance, CostComponents& costs) const {
    numberOfTotalCompute++;

    if(unAssignedExamList.size()) {
        computeCostPartial(instance, costs);
        return;
    }

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

    // pairs of exams with students in common, edges in graph

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

            if(periodObjects.first.dateid == periodObjects.second.dateid) {
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

void ExamTTSolution::computeCostPartial(ExamTT const& instance, CostComponents& costs) const {
    numberOfTotalCompute++;
    int E = instance.exams.size();
    int R = instance.rooms.size();
    int P = instance.periods.size();

    InstitutionalWeightings const& institutionalWeightings = instance.institutionalWeightings;
    FrontLoadParams const& frontload = institutionalWeightings.frontload;

    // period constraints : pair of exams

    costs.hard.periodConstraintAfter = 0;
    costs.hard.periodConstraintExclusion = 0;
    costs.hard.periodConstraintCoincidence = 0;

    for(Two<ExamId> exams : instance.examsAfter) if(isAssigned(exams.first) && isAssigned(exams.second))
        if(periods[exams.first] <= periods[exams.second])
            costs.hard.periodConstraintAfter++;

    for(Two<ExamId> exams : instance.examsExclusion) if(isAssigned(exams.first) && isAssigned(exams.second))
        if(periods[exams.first] == periods[exams.second])
            costs.hard.periodConstraintExclusion++;

    for(Two<ExamId> exams : instance.examsCoincidence) if(isAssigned(exams.first) && isAssigned(exams.second))
        if(periods[exams.first] != periods[exams.second])
            costs.hard.periodConstraintCoincidence++;

    // pairs of exams with students in common, edges in graph

    costs.hard.simultaneousExams = 0;

    costs.soft.twoExamsInARow = 0;
    costs.soft.twoExamsInADay = 0;
    costs.soft.periodSpread = 0;

    for(ExamId i = 0; i < E; i++)
    if(isAssigned(i))
    for(ExamId j = i + 1; j < E; j++)
    if(isAssigned(j))
    if(instance.numberStudentsInCommon[i][j]) {
        Two<PeriodId> periodsId = periodsOf(i,j);
        Two<Period const&> periodObjects = instance.periodsOf(periodsId);

        Cost cost = instance.numberStudentsInCommon[i][j];
        if(periodsId.first == periodsId.second) {
            costs.hard.simultaneousExams += cost;
        } else {
            int diff = abs(periodsId.first - periodsId.second);

            if(periodObjects.first.dateid == periodObjects.second.dateid) {
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

    for(ExamId e : instance.frontLoadExams) if(isAssigned(e))
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
    for(ExamId i = 0; i < E; i++) if(isFullyAssigned(i)) {
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

    for(ExamId i : instance.examsRoomsExclusive) if(isFullyAssigned(i)) {
        auto a = assignement(i);

        for(ExamId j = 0; j < E; j++) if(isFullyAssigned(j)) {
            if(a == assignement(j) && i != j) {
                costs.hard.roomConstraint++;
                break;
            }
        }
    }
}

void ExamTT::presentation(ostream & log) const {
    auto& inst = *this;
    auto& instance = *this;

    int E = inst.exams.size(),
        P = inst.periods.size(),
        R = inst.rooms.size();

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

    log << "* Exams : minutes | id | size | connected exams : student in commons : {exam1 exam2} [studentsWith1 studentsWith2]" << endl;
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
            << setw(3) << "[" + to_string(instance.hasStudentsInCommonOfExam(i).size()) << "] " << EPrefix(i) << " : "
            << /* print(M, nocomma | nobraces) << */ print(instance.hasStudentsInCommonOfExam(i), nocomma) << " " << print(M2, nocomma) << endl;
    }

    std::set<ExamId> setFrontloadExams(instance.frontLoadExams.begin(), instance.frontLoadExams.end());
    log << "* Frontload exams (" << instance.frontLoadExams.size() << ") : " << print(setFrontloadExams, nocomma) << " : ";
    for(ExamId i : instance.frontLoadExams)
        log << EPrefix(i) << "(" << inst.exams[i].students.size() << ") ";
    log << endl;

    log << "* Frontload periods (" << inst.institutionalWeightings.frontload.time << ") : "
        << PPrefix(P - inst.institutionalWeightings.frontload.time) << " - " << PPrefix(P-1) << endl;

    log << "* ROOM_EX Exams that must be room exclusive (" << inst.examsRoomsExclusive.size() << ")" << endl;
    for(ExamId e : inst.examsRoomsExclusive)
        log << e << endl;

    log << "* InstitutionalWeightings" << endl
         << inst.institutionalWeightings << endl;
}

void ExamTT::testDelta(ExamTTSolution& sol, std::ostream& log, int N, bool checkEachMove) const {
    auto& inst = *this;

    int E = inst.exams.size(),
        P = inst.periods.size(),
        R = inst.rooms.size();

    Random ran;
    for(int i = 0; i < N; i++) {
        int e = ran.randrange(E), p = ran.randrange(P), r = ran.randrange(R);
        int bp = sol.periods[e], br = sol.rooms[e];
        sol.move(inst, e,p,r);

        if(checkEachMove) {
            if(! sol.costs.exactlyEqual(sol.computeAndGetCost(inst))) {
                log << "Error " << i << " Real " << sol.computeAndGetCost(inst).print(inst)
                     << "vs diffed " << sol.costs.print(inst) << endl
                     << " e p r; bp br = " << e << " " << p << " " << r << ";" << bp << " " << br << endl;

                sol.move(inst, e,bp,br);
                log << "Before" << endl << sol.computeAndGetCost(inst).print(inst); sol.printTimelineTo(inst, log);

                sol.move(inst, e,p,r);
                log << "After" << endl << sol.computeAndGetCost(inst).print(inst); sol.printTimelineTo(inst, log);

                cout << "Error in random move " << i << " see log." << endl;
                exit(1);
            }
        }
    }


    auto real = sol.computeAndGetCost(inst);
    log << "End of random moves" << endl
        << "Real" << real.print(inst)
        << "Diffed" << sol.costs.print(inst) << endl;

    if(! sol.costs.exactlyEqual(real))
        cout << "Error " << endl;
}

string lower(string s) {
    // first time we Want a copy...
    for(char& c : s)
        if('A' <= c && c <= 'Z')
            c += 'a' - 'A';
    return s;
}

void interactiveSolution(ExamTTSolution& sol, ExamTT const& inst) {
    auto& instance = inst;

    bool isHelp = true;
    string str;
    for(;;) {
        if(isHelp)
            cout
             << "Change solution, C in front means compute and print difference in cost" << endl
             << "[p | cp] Period <exam> <period>      : move <exam> to period <period>" << endl
             << "[r | cr] Room <exam> <room>          : move <exam> to room <room>" << endl
             << "[m | cm] Move <exam> <period> <room> : move <exam> to <period> <room>" << endl
             << "[s | cs] Swap <exam1> <exam2>        : swap <exam> and <exam>" << endl
             << "     [d] Display                     : display current solution" << endl
             << "     [w] Write <filename>            : write current solution to <file>" << endl
             << "    [wv] WriteVerbose <filename>     : write current solution verbosly <file>" << endl
             << "    [dc] DCost                       : display current solution cost" << endl
             << "    [rc] ReCalculate                 : recalculate cost of current solution" << endl
             << "     [h] Help                        : display help" << endl
             << "     [q] Quit                        : quit" << endl
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
        bool isQuit         = name == "q" || name == "quit" || name == "exit";
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
                    cout << "Wrong exam id or period id !";
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

        } else if(isQuit){
            break;
        } else {
            cout << "Unknown command '" << name << "'" << endl;
        }
    }
}

void test::delta(ExamTT const& inst, int N, bool checkEachMove) {
    ofstream log("log");
    inst.presentation(log);

    ExamTTSolution sol;
    Random ran;
    sol.initRandom(inst, ran);

    log << endl << "* Solution" << endl;
    sol.printTo(inst, log);

    inst.testDelta(sol, log, N, checkEachMove);

    log.close();

    cout << "Printed to log" << endl;
}

void test::deltaRemoveAdd(ExamTT const& inst, int N, bool checkEachMove, int G, bool checkEachMovePartial) {
    ExamTTSolution sol;
    Random ran;
    sol.initRandom(inst, ran);

    int E = inst.exams.size(),
        P = inst.periods.size(),
        R = inst.rooms.size();

    if(G > E)
        throw std::invalid_argument("G > E");

    vector<int> bps(G), brs(G);
    vector<int> ps(G), rs(G);
    std::vector<ExamId> es(G);

    for(int i = 0; i < N; i++) {

        if(G == 1) {
            int e = ran.randrange(E), p = ran.randrange(P), r = ran.randrange(R);
            int bp = sol.periods[e], br = sol.rooms[e];

            sol.removeExam(inst, e);
            sol.addExam(inst, e, p, r);

            if(checkEachMove) {
                if(! sol.costs.exactlyEqual(sol.computeAndGetCost(inst))) {
                    cerr << "Error " << i << " Real " << sol.computeAndGetCost(inst).print(inst)
                         << "vs diffed " << sol.costs.print(inst) << " = " << (sol.costs - sol.computeAndGetCost(inst)).print(inst) << endl
                         << " e p r; bp br = " << e << " " << p << " " << r << ";" << bp << " " << br << endl;

                    exit(1);
                }
            }
        } else {
            std::set<ExamId> removed;

            for(int j = 0; j < G; j++) {
                int e = ran.randrange(E);
                while(removed.count(e))
                    e = ran.randrange(E);

                bps[j] = sol.periods[e];
                brs[j] = sol.rooms[e];
                es[j] = e;

                sol.removeExam(inst, e);

                removed.insert(e);
            }

            for(int j = 0; j < G; j++) {
                int p = ps[j] = ran.randrange(P), r = rs[j] = ran.randrange(R);
                sol.addExam(inst, es[j], p, r);
            }

            if(checkEachMove) {
                if(! sol.costs.exactlyEqual(sol.computeAndGetCost(inst))) {
                    cerr << "Error " << i << " Real " << sol.computeAndGetCost(inst).print(inst)
                         << "vs diffed " << sol.costs.print(inst) << " = " << (sol.costs - sol.computeAndGetCost(inst)).print(inst) << endl
                         << " es ps rs; bps brs = " << es << " " << ps << " " << rs << "; " << bps << " " << brs << endl;

                    exit(1);
                }
            }
        }
    }

    auto real = sol.computeAndGetCost(inst);
    if(! sol.costs.exactlyEqual(real))
        cerr << "Error " << endl;
}

void test::iterateVsComputeStep(ExamTTSolution* solution, Neighborhood* neigh) {
    ExamTTSolution* sol1 = (ExamTTSolution*) solution->clone();
    std::set<std::pair<std::vector<int>, std::vector<int>>> first, second;

    neigh->iterate(solution, [sol1,&first](){
        first.insert(std::make_pair(sol1->periods, sol1->rooms));
    });

    neigh->reset();

    ExamTTSolution* sol2 = (ExamTTSolution*) solution->clone();
    for(auto it = neigh->begin(sol2); it != neigh->end(); ++it)
        second.insert(std::make_pair(sol2->periods, sol2->rooms));

    if(first != second)
        throw std::invalid_argument("two ways of iterating differ ! " + to_string(first.size()) + " vs " + to_string(second.size()) + " neighbors");
}

void test::constructUnassignedVsRemoveAll(ExamTT& instance) {
    ExamTTSolution sol1;
    sol1.initFromZeroPeriodAndZeroRoom(instance);
    for(int i = 0; i < instance.E(); ++i)
        sol1.removeExam(instance, i);

    ExamTTSolution sol2;
    sol2.initUnassigned(instance);

    std::set<int> A1(sol1.unAssignedExamList.begin(), sol1.unAssignedExamList.end());
    std::set<int> A2(sol2.unAssignedExamList.begin(), sol2.unAssignedExamList.end());

    if(!(true
         && A1.size() == sol1.unAssignedExamList.size()
         && A2.size() == sol2.unAssignedExamList.size()
         && A1.size() == instance.E()
         && A2.size() == instance.E()
         && A1 == A2
         && sol1.periods == sol2.periods
         && sol1.rooms == sol2.rooms
         && std::all_of(sol1.periods.begin(), sol1.periods.end(), [](int i){ return i == -1; })
         && std::all_of(sol1.rooms.begin(), sol1.rooms.end(), [](int i){ return i == -1; })
    )) {
        throw std::invalid_argument("error");
    }
}

namespace test {
void interactive(ExamTT const& inst) {
    ExamTTSolution sol;
    Random ran;
    sol.initRandom(inst, ran);
    interactiveSolution(sol, inst);
}
}

ExamTTSolution::~ExamTTSolution()
{

}

const void* ExamTTSolution::getRawData() const {
    return (const void*) this;
}

void ExamTTSolution::setRawData(const void *data) {
    numberOfClones++;

    ExamTTSolution const * other = (ExamTTSolution const *) data;

    setSolutionValue(const_cast<ExamTTSolution*>(other)->getSolutionValue());

    costs = other->costs; // O(1)
    periods = other->periods; // O(P)
    rooms = other->rooms; // O(R)

    // structures

    if(USE_LAZY_STRUCTURES) {
        hasStructures = false;
        return;
    }

    examsByPeriodRoomIterators.resize(periods.size()); // O(1) E

    // inner copy

    examsByPeriodRoom = other->examsByPeriodRoom; // O(P R + E)
    durationColorUsed = other->durationColorUsed; // O(P R NoD)
    unAssignedExamList = other->unAssignedExamList;

    // O(P R + E) because sum(list.size() for row in examsByPeriodRoom for list in row) == E
    for(auto & a : examsByPeriodRoom)
        for(auto & b : a)
            for(auto it = b.begin(); it != b.end(); ++it)
                examsByPeriodRoomIterators[*it] = it;
}

Solution* ExamTTSolution::clone() {
    // The difference with setRawData is that new solution is empty. Here it doesn't change anything

    ExamTTSolution* other = new ExamTTSolution(getSolutionValue());

    other->setRawData(this);

    return other;
}

void ExamTTSolution::swap(Solution * rawOther) {
    ExamTTSolution* other = (ExamTTSolution*) rawOther;

    // swap cost
    auto c = other->getSolutionValue();
    other->setSolutionValue(getSolutionValue());
    setSolutionValue(c);

    // swap vectors
    std::swap(other->costs, costs);
    other->periods.swap(periods);
    other->rooms.swap(rooms);

    /*
    other->examsByPeriodsIterators.swap(examsByPeriodsIterators);
    other->examsByRoomsIterators.swap(examsByRoomsIterators);
    other->examsByPeriods.swap(examsByPeriods);
    other->examsByRooms.swap(examsByRooms);
    */

    other->examsByPeriodRoom.swap(examsByPeriodRoom);
    other->examsByPeriodRoomIterators.swap(examsByPeriodRoomIterators);
    other->durationColorUsed.swap(durationColorUsed);
}

string ExamTTSolution::getSolutionRepresentation() {
    ostringstream oss;

    oss << "assign = [";

    if(periods.size() && rooms.size()) {
        auto it1 = periods.begin();
        auto it2 = rooms.begin();
        oss << "[" << *it1++ << ", " << *it2++ << "]";
        while(it1 != periods.end() && it2 != rooms.end())
            oss << ", [" << *it1++ << ", " << *it2++ << "]";
    }

    oss << "]";

    oss << "; hard = " << costs.hard.make_tuple() << "; soft = " << costs.soft.make_tuple() << ";";
    oss << " cost_tuple = (" << costs.hard.sum() << ", " << costs.soft.sum() << ");";

    return oss.str();
}

/*
 * Methods
 */

void ExamTTSolution::buildStructures(InstanceRef instance)
{
    int P = instance.periods.size();
    int E = instance.exams.size();
    int R = instance.rooms.size();

    examsByPeriodRoomIterators.resize(E);
    examsByPeriodRoom.assign(P, std::vector<std::list<ExamId>>(R, std::list<ExamId>()));

    durationColorUsed.assign(P, std::vector<MapVec<Color,int>>(R, std::vector<int>(instance.numberOfDurations(), 0)));

    if(unAssignedExamList.size()) {
        if(unAssignedExamList.size() == E) {

        } else {

            for(ExamId e = 0; e < E; e++)
                if(isFullyAssigned(e))
                    examsByPeriodRoomIterators[e] = examsByPeriodRoom[periods[e]][rooms[e]].insert(examsByPeriodRoom[periods[e]][rooms[e]].end(), e);

            for(auto it = unAssignedExamList.begin(); it != unAssignedExamList.end(); ++it)
                examsByPeriodRoomIterators[*it] = it;

            durationColorUsed.assign(P, std::vector<MapVec<Color,int>>(R, std::vector<int>(instance.numberOfDurations(), 0)));
            for(ExamId e = 0; e < E; e++)
                durationColorUsed[periods[e]][rooms[e]][instance.exams[e].durationColor]++;
        }
    } else {

        for(ExamId e = 0; e < E; e++)
            examsByPeriodRoomIterators[e] = examsByPeriodRoom[periods[e]][rooms[e]].insert(examsByPeriodRoom[periods[e]][rooms[e]].end(), e);

        for(ExamId e = 0; e < E; e++)
            durationColorUsed[periods[e]][rooms[e]][instance.exams[e].durationColor]++;
    }

    hasStructures = true;
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

    if(prevP == nextP && prevR == nextR) // if(make_pair(prevP, prevR) == make_pair(nextP, nextR))
        return;

    InstitutionalWeightings const& weightings = instance.institutionalWeightings;

    Instance::ExamsRelated const& related = instance.examsRelated[ex];

    int s = nextP > prevP ? +1 :
            nextP < prevP ? -1 : 0;

    // if the period has changed
    if(s) {

        // period penalty
        costs.soft.periodsPenalty -= instance.periods[prevP].penalty;
        costs.soft.periodsPenalty += instance.periods[nextP].penalty;

        // after
        if(1) {
            // remove
            for(ExamId j : related.afters)
                if(periods[j] >= prevP)
                    costs.hard.periodConstraintAfter--;

            for(ExamId j : related.befores)
                if(periods[j] <= prevP)
                    costs.hard.periodConstraintAfter--;

            // add
            for(ExamId j : related.afters)
                if(periods[j] >= nextP)
                    costs.hard.periodConstraintAfter++;

            for(ExamId j : related.befores)
                if(periods[j] <= nextP)
                    costs.hard.periodConstraintAfter++;

        } else {
            int m = min(prevP, nextP);
            int M = max(prevP, nextP);

            for(ExamId j : related.afters) // very small vector : [min_max_mean([len(i.related[j].afters) for j in range(i.E)]) for i in I]
                if(m <= periods[j] && periods[j] < M)
                    costs.hard.periodConstraintAfter -= s;

            for(ExamId j : related.befores)
                if(m < periods[j] && periods[j] <= M)
                    costs.hard.periodConstraintAfter += s;
        }

        // exclusion/coincidences

        //// remove
        for(ExamId j : related.exclusions) // very small vector : [min_max_mean([len(i.related[j].exclusions) for j in range(i.E)]) for i in I]
            if(periods[j] == prevP)
                costs.hard.periodConstraintExclusion--;

        for(ExamId j : related.coincidences) // very small vector : [min_max_mean([len(i.related[j].coincidences) for j in range(i.E)]) for i in I]
            if(periods[j] != prevP)
                costs.hard.periodConstraintCoincidence--;

        //// add
        for(ExamId j : related.exclusions)
            if(periods[j] == nextP)
                costs.hard.periodConstraintExclusion++;

        for(ExamId j : related.coincidences)
            if(periods[j] != nextP)
                costs.hard.periodConstraintCoincidence++;
    }

    // room penalty
    costs.soft.roomsPenalty -= instance.rooms[prevR].penalty;
    costs.soft.roomsPenalty += instance.rooms[nextR].penalty;

    // about students in common : (twoExamsInARow, twoExamsInADay, periodSpread, simultaneousExams)

    for(ExamId j : related.hasStudentsInCommon) if(periods[j] != -1) {
        Cost cost = instance.numberStudentsInCommon[ex][j];

        int relatedP = periods[j];

        // remove

        if(prevP == relatedP) {
            costs.hard.simultaneousExams -= cost;
        } else {
            int prevDiff = abs(relatedP - prevP);

            if(instance.periods[relatedP].dateid == instance.periods[prevP].dateid) {
                if(prevDiff == 1)
                    costs.soft.twoExamsInARow -= cost * weightings.twoInARow;
                else
                    costs.soft.twoExamsInADay -= cost * weightings.twoInADay;
            }

            if(prevDiff <= weightings.periodSpread)
                costs.soft.periodSpread -= cost;
        }

        // add

        if(nextP == relatedP) {
            costs.hard.simultaneousExams += cost;
        } else {
            int nextDiff = abs(relatedP - nextP);

            if(instance.periods[relatedP].dateid == instance.periods[nextP].dateid) {
                if(nextDiff == 1)
                    costs.soft.twoExamsInARow += cost * weightings.twoInARow;
                else
                    costs.soft.twoExamsInADay += cost * weightings.twoInADay;
            }

            if(nextDiff <= weightings.periodSpread)
                costs.soft.periodSpread += cost;
        }
    }

    // frontload
    if(instance.isInFrontLoad[ex]) {
        if(instance.periodInTheEnd(prevP))
            costs.soft.frontload -= weightings.frontload.penalty;
        if(instance.periodInTheEnd(nextP))
            costs.soft.frontload += weightings.frontload.penalty;
    }

    // excessiveDuration
    {
        if(instance.exams[ex].duration > instance.periods[prevP].duration)
            costs.hard.excessiveDuration--;
        if(instance.exams[ex].duration > instance.periods[nextP].duration)
            costs.hard.excessiveDuration++;
    }

    // comparing assignement before and after

    // room exclusive
    std::list<ExamId> const& leavingAndMe = examsByPeriodRoom[prevP][prevR];
    std::list<ExamId> const& meeting = examsByPeriodRoom[nextP][nextR];

    if(leavingAndMe.size() - 1 == 1 && instance.examIsRoomExclusive[leavingAndMe.front() == ex ? leavingAndMe.back() : leavingAndMe.front()])
        costs.hard.roomConstraint--;
    if(leavingAndMe.size() - 1 > 0 && instance.examIsRoomExclusive[ex])
        costs.hard.roomConstraint--;

    if(meeting.size() == 1 && instance.examIsRoomExclusive[*meeting.begin()])
        costs.hard.roomConstraint++;
    if(meeting.size() > 0 && instance.examIsRoomExclusive[ex])
        costs.hard.roomConstraint++;


    // overcapacity
    {
        int mySize = instance.exams[ex].students.size();

        // remove
        int capLeaving = instance.rooms[prevR].capacity;

        int leavingSum = 0; // = sum(instance.exams[j].students.size() for j in leaving)
        for(ExamId j : leavingAndMe)
            leavingSum += instance.exams[j].students.size();
        leavingSum -= instance.exams[ex].students.size();

        if(leavingSum + mySize <= capLeaving) // was not in overcapacity
            ;
        else if(leavingSum > capLeaving) // was in overcapacity, still in overcapacity
            costs.hard.overCapacity -= mySize;
        else // was in overcapacity, not in overcapacity now
            costs.hard.overCapacity -= leavingSum + mySize - capLeaving;

        // add
        int capMeeting = instance.rooms[nextR].capacity;

        int meetingSum = 0; // = sum(instance.exams[j].students.size() for j in meeting)

        std::list<ExamId> const& meeting = examsByPeriodRoom[nextP][nextR];

        for(ExamId j : meeting)
            meetingSum += instance.exams[j].students.size();

        if(meetingSum > capMeeting) // already over capacity
            costs.hard.overCapacity += mySize;
        else if(meetingSum + mySize > capMeeting)
            costs.hard.overCapacity += meetingSum + mySize - capMeeting;
        else
            ; // no overcapacity created
    }

    // mixedDurations
    {
        int myDuration = instance.exams[ex].duration;

        if(USE_COLOR_STRUCTURE) {
            // no heap allocated, incremental info

            Color myDurationColor = instance.exams[ex].durationColor;
            auto const& leavingDuration = durationColorUsed[prevP][prevR];
            auto const& meetingDuration = durationColorUsed[nextP][nextR];

            const_cast<MapVec<Color,int>&>(leavingDuration)[myDurationColor]--;

            int countLeaving = 0;
            for(int x : leavingDuration)
                if(x > 0)
                    countLeaving++;

            if(leavingDuration[myDurationColor] == 0 && countLeaving > 0)
                costs.soft.mixedDuration -= weightings.nonMixedDurations;

            const_cast<MapVec<Color,int>&>(meetingDuration)[myDurationColor]++;
            int countMeeting = 0;
            for(int x : meetingDuration)
                if(x > 0)
                    countMeeting++;

            if(meetingDuration[myDurationColor] == 1 && countMeeting > 1)
                costs.soft.mixedDuration += weightings.nonMixedDurations;
            const_cast<MapVec<Color,int>&>(meetingDuration)[myDurationColor]--;

            const_cast<MapVec<Color,int>&>(leavingDuration)[myDurationColor]++;

        } else {
            int NLeavingWithoutMe, NMeetingWithoutMe, NLeavingWithMe, NMeetingWithMe;

            if(1) {
                // set implementation
                std::set<Minutes> durationsLeaving;
                std::set<Minutes> durationsMeeting;

                for(ExamId j : leavingAndMe) if(j != ex)
                    durationsLeaving.insert(instance.exams[j].duration);

                for(ExamId j : meeting)
                    durationsMeeting.insert(instance.exams[j].duration);

                NLeavingWithoutMe = durationsLeaving.size();
                NMeetingWithoutMe = durationsMeeting.size();
                NLeavingWithMe = durationsLeaving.size() + (durationsLeaving.count(myDuration) ? 0 : 1);
                NMeetingWithMe = durationsMeeting.size() + (durationsMeeting.count(myDuration) ? 0 : 1);
            } else {
                // N² implementation, no heap allocated
                int durationsLeavingDuplicates = 0;
                bool myDurationInLeaving = false;
                for(auto i = leavingAndMe.begin(); i != leavingAndMe.end(); ++i) {
                    if(*i == ex)
                        continue;
                    if(instance.exams[*i].duration == myDuration)
                        myDurationInLeaving = true;

                    auto j = i;
                    ++j;
                    for( ; j != leavingAndMe.end(); ++j) {
                        if(*j == ex)
                            continue;

                        if(instance.exams[*i].duration == instance.exams[*j].duration) {
                            ++durationsLeavingDuplicates;
                            break;
                        }
                    }
                }

                int durationsMeetingDuplicates = 0;
                bool myDurationInMeeting = false;
                for(auto i = meeting.begin(); i != meeting.end(); ++i) {
                    if(instance.exams[*i].duration == myDuration)
                        myDurationInMeeting = true;
                    auto j = i;
                    ++j;
                    for( ; j != meeting.end(); ++j) {
                        if(instance.exams[*i].duration == instance.exams[*j].duration) {
                            ++durationsMeetingDuplicates;
                            break;
                        }
                    }
                }

                int durationsLeavingSize = leavingAndMe.size() - 1 - durationsLeavingDuplicates;
                int durationsMeetingSize = meeting.size() - durationsMeetingDuplicates;

                NLeavingWithoutMe = durationsLeavingSize;
                NMeetingWithoutMe = durationsMeetingSize;
                NLeavingWithMe = durationsLeavingSize + (myDurationInLeaving ? 0 : 1);
                NMeetingWithMe = durationsMeetingSize + (myDurationInMeeting ? 0 : 1);
            }

            if(NLeavingWithoutMe > 0)
                costs.soft.mixedDuration -= (NLeavingWithMe - NLeavingWithoutMe) * weightings.nonMixedDurations;
            if(NMeetingWithoutMe > 0)
                costs.soft.mixedDuration += (NMeetingWithMe - NMeetingWithoutMe) * weightings.nonMixedDurations;
        }
    }
}

void ExamTTSolution::updateRemoveOutOfStructures(InstanceRef instance, ExamId ex, PeriodId prevP, RoomId prevR) {
    InstitutionalWeightings const& weightings = instance.institutionalWeightings;

    Instance::ExamsRelated const& related = instance.examsRelated[ex];

    // period penalty
    costs.soft.periodsPenalty -= instance.periods[prevP].penalty;

    // remove
    for(ExamId j : related.afters) if(periods[j] != -1)
        if(periods[j] >= prevP)
            costs.hard.periodConstraintAfter--;

    for(ExamId j : related.befores) if(periods[j] != -1)
        if(periods[j] <= prevP)
            costs.hard.periodConstraintAfter--;

    // exclusion/coincidences

    //// remove
    for(ExamId j : related.exclusions) if(periods[j] != -1) // very small vector : [min_max_mean([len(i.related[j].exclusions) for j in range(i.E)]) for i in I]
        if(periods[j] == prevP)
            costs.hard.periodConstraintExclusion--;

    for(ExamId j : related.coincidences) if(periods[j] != -1) // very small vector : [min_max_mean([len(i.related[j].coincidences) for j in range(i.E)]) for i in I]
        if(periods[j] != prevP)
            costs.hard.periodConstraintCoincidence--;

    // room penalty
    costs.soft.roomsPenalty -= instance.rooms[prevR].penalty;

    // about students in common : (twoExamsInARow, twoExamsInADay, periodSpread, simultaneousExams)

    for(ExamId j : related.hasStudentsInCommon) {
        Cost cost = instance.numberStudentsInCommon[ex][j];

        int relatedP = periods[j];

        // remove

        if(prevP == relatedP) {
            costs.hard.simultaneousExams -= cost;
        } else {
            int prevDiff = abs(relatedP - prevP);

            if(instance.periods[relatedP].dateid == instance.periods[prevP].dateid) {
                if(prevDiff == 1)
                    costs.soft.twoExamsInARow -= cost * weightings.twoInARow;
                else
                    costs.soft.twoExamsInADay -= cost * weightings.twoInADay;
            }

            if(prevDiff <= weightings.periodSpread)
                costs.soft.periodSpread -= cost;
        }
    }

    // frontload
    if(instance.isInFrontLoad[ex])
        if(instance.periodInTheEnd(prevP))
            costs.soft.frontload -= weightings.frontload.penalty;

    // excessiveDuration
    if(instance.exams[ex].duration > instance.periods[prevP].duration)
        costs.hard.excessiveDuration--;

    // comparing assignement before and after

    // room exclusive
    std::list<ExamId> const& leaving = examsByPeriodRoom[prevP][prevR];

    if(leaving.size() == 1 && instance.examIsRoomExclusive[leaving.front()])
        costs.hard.roomConstraint--;
    if(leaving.size() > 0 && instance.examIsRoomExclusive[ex])
        costs.hard.roomConstraint--;

    // overcapacity
    {
        int mySize = instance.exams[ex].students.size();

        // remove
        int capLeaving = instance.rooms[prevR].capacity;

        int leavingSum = 0; // = sum(instance.exams[j].students.size() for j in leaving)
        for(ExamId j : leaving)
            leavingSum += instance.exams[j].students.size();

        if(leavingSum + mySize <= capLeaving) // was not in overcapacity
            ;
        else if(leavingSum > capLeaving) // was in overcapacity, still in overcapacity
            costs.hard.overCapacity -= mySize;
        else // was in overcapacity, not in overcapacity now
            costs.hard.overCapacity -= leavingSum + mySize - capLeaving;
    }

    // mixedDurations
    {
        int myDuration = instance.exams[ex].duration;

        if(USE_COLOR_STRUCTURE) {
            // no heap allocated, incremental info

            Color myDurationColor = instance.exams[ex].durationColor;
            auto const& leavingDuration = durationColorUsed[prevP][prevR];

            int countLeaving = 0;
            for(int x : leavingDuration)
                if(x > 0)
                    countLeaving++;

            if(leavingDuration[myDurationColor] == 0 && countLeaving > 0)
                costs.soft.mixedDuration -= weightings.nonMixedDurations;

        } else {
            int NLeavingWithoutMe, NLeavingWithMe;

            if(1) {
                // set implementation
                std::set<Minutes> durationsLeaving;

                for(ExamId j : leaving)
                    durationsLeaving.insert(instance.exams[j].duration);

                NLeavingWithoutMe = durationsLeaving.size();
                NLeavingWithMe = durationsLeaving.size() + (durationsLeaving.count(myDuration) ? 0 : 1);
            } else {
                // N² implementation, no heap allocated
                throw std::invalid_argument("not implemented");
            }

            if(NLeavingWithoutMe > 0)
                costs.soft.mixedDuration -= (NLeavingWithMe - NLeavingWithoutMe) * weightings.nonMixedDurations;
        }
    }
}

void ExamTTSolution::updateRemove(InstanceRef instance, ExamId ex) {
    PeriodId prevP = periods[ex];
    RoomId prevR = rooms[ex];

    InstitutionalWeightings const& weightings = instance.institutionalWeightings;

    Instance::ExamsRelated const& related = instance.examsRelated[ex];

    // period penalty
    costs.soft.periodsPenalty -= instance.periods[prevP].penalty;

    // remove
    for(ExamId j : related.afters) if(periods[j] != -1)
        if(periods[j] >= prevP)
            costs.hard.periodConstraintAfter--;

    for(ExamId j : related.befores) if(periods[j] != -1)
        if(periods[j] <= prevP)
            costs.hard.periodConstraintAfter--;

    // exclusion/coincidences

    //// remove
    for(ExamId j : related.exclusions) if(periods[j] != -1) // very small vector : [min_max_mean([len(i.related[j].exclusions) for j in range(i.E)]) for i in I]
        if(periods[j] == prevP)
            costs.hard.periodConstraintExclusion--;

    for(ExamId j : related.coincidences) if(periods[j] != -1) // very small vector : [min_max_mean([len(i.related[j].coincidences) for j in range(i.E)]) for i in I]
        if(periods[j] != prevP)
            costs.hard.periodConstraintCoincidence--;

    // room penalty
    costs.soft.roomsPenalty -= instance.rooms[prevR].penalty;

    // about students in common : (twoExamsInARow, twoExamsInADay, periodSpread, simultaneousExams)

    for(ExamId j : related.hasStudentsInCommon) if(periods[j] != -1) {
        Cost cost = instance.numberStudentsInCommon[ex][j];

        int relatedP = periods[j];

        // remove

        if(prevP == relatedP) {
            costs.hard.simultaneousExams -= cost;
        } else {
            int prevDiff = abs(relatedP - prevP);

            if(instance.periods[relatedP].dateid == instance.periods[prevP].dateid) {
                if(prevDiff == 1)
                    costs.soft.twoExamsInARow -= cost * weightings.twoInARow;
                else
                    costs.soft.twoExamsInADay -= cost * weightings.twoInADay;
            }

            if(prevDiff <= weightings.periodSpread)
                costs.soft.periodSpread -= cost;
        }
    }

    // frontload
    if(instance.isInFrontLoad[ex])
        if(instance.periodInTheEnd(prevP))
            costs.soft.frontload -= weightings.frontload.penalty;

    // excessiveDuration
    if(instance.exams[ex].duration > instance.periods[prevP].duration)
        costs.hard.excessiveDuration--;

    // comparing assignement before and after

    // room exclusive
    std::list<ExamId> const& leavingAndMe = examsByPeriodRoom[prevP][prevR];

    if(leavingAndMe.size() - 1 == 1 && instance.examIsRoomExclusive[leavingAndMe.front() == ex ? leavingAndMe.back() : leavingAndMe.front()])
        costs.hard.roomConstraint--;
    if(leavingAndMe.size() - 1 > 0 && instance.examIsRoomExclusive[ex])
        costs.hard.roomConstraint--;

    // overcapacity
    {
        int mySize = instance.exams[ex].students.size();

        // remove
        int capLeaving = instance.rooms[prevR].capacity;

        int leavingSum = 0; // = sum(instance.exams[j].students.size() for j in leaving)
        for(ExamId j : leavingAndMe)
            leavingSum += instance.exams[j].students.size();
        leavingSum -= instance.exams[ex].students.size();

        if(leavingSum + mySize <= capLeaving) // was not in overcapacity
            ;
        else if(leavingSum > capLeaving) // was in overcapacity, still in overcapacity
            costs.hard.overCapacity -= mySize;
        else // was in overcapacity, not in overcapacity now
            costs.hard.overCapacity -= leavingSum + mySize - capLeaving;
    }

    // mixedDurations
    {
        int myDuration = instance.exams[ex].duration;

        if(USE_COLOR_STRUCTURE) {
            // no heap allocated, incremental info

            Color myDurationColor = instance.exams[ex].durationColor;
            auto const& leavingDuration = durationColorUsed[prevP][prevR];

            const_cast<MapVec<Color,int>&>(leavingDuration)[myDurationColor]--;

            int countLeaving = 0;
            for(int x : leavingDuration)
                if(x > 0)
                    countLeaving++;

            if(leavingDuration[myDurationColor] == 0 && countLeaving > 0)
                costs.soft.mixedDuration -= weightings.nonMixedDurations;

            const_cast<MapVec<Color,int>&>(leavingDuration)[myDurationColor]++;

        } else {
            int NLeavingWithoutMe, NLeavingWithMe;

            if(1) {
                // set implementation
                std::set<Minutes> durationsLeaving;

                for(ExamId j : leavingAndMe) if(j != ex)
                    durationsLeaving.insert(instance.exams[j].duration);

                NLeavingWithoutMe = durationsLeaving.size();
                NLeavingWithMe = durationsLeaving.size() + (durationsLeaving.count(myDuration) ? 0 : 1);
            } else {
                // N² implementation, no heap allocated
                throw std::invalid_argument("not implemented");
            }

            if(NLeavingWithoutMe > 0)
                costs.soft.mixedDuration -= (NLeavingWithMe - NLeavingWithoutMe) * weightings.nonMixedDurations;
        }
    }
}

void ExamTTSolution::updateAdd(InstanceRef instance, ExamId ex, PeriodId nextP, RoomId nextR) {
    InstitutionalWeightings const& weightings = instance.institutionalWeightings;

    Instance::ExamsRelated const& related = instance.examsRelated[ex];

    // period penalty
    costs.soft.periodsPenalty += instance.periods[nextP].penalty;

    // add
    for(ExamId j : related.afters) if(periods[j] != -1)
        if(periods[j] >= nextP)
            costs.hard.periodConstraintAfter++;

    for(ExamId j : related.befores) if(periods[j] != -1)
        if(periods[j] <= nextP)
            costs.hard.periodConstraintAfter++;

    // exclusion/coincidences

    //// add
    for(ExamId j : related.exclusions) if(periods[j] != -1)
        if(periods[j] == nextP)
            costs.hard.periodConstraintExclusion++;

    for(ExamId j : related.coincidences) if(periods[j] != -1)
        if(periods[j] != nextP)
            costs.hard.periodConstraintCoincidence++;

    // room penalty
    costs.soft.roomsPenalty += instance.rooms[nextR].penalty;

    // about students in common : (twoExamsInARow, twoExamsInADay, periodSpread, simultaneousExams)

    for(ExamId j : related.hasStudentsInCommon) if(periods[j] != -1) {
        Cost cost = instance.numberStudentsInCommon[ex][j];

        int relatedP = periods[j];

        // add

        if(nextP == relatedP) {
            costs.hard.simultaneousExams += cost;
        } else {
            int nextDiff = abs(relatedP - nextP);

            if(instance.periods[relatedP].dateid == instance.periods[nextP].dateid) {
                if(nextDiff == 1)
                    costs.soft.twoExamsInARow += cost * weightings.twoInARow;
                else
                    costs.soft.twoExamsInADay += cost * weightings.twoInADay;
            }

            if(nextDiff <= weightings.periodSpread)
                costs.soft.periodSpread += cost;
        }
    }

    // frontload
    if(instance.isInFrontLoad[ex])
        if(instance.periodInTheEnd(nextP))
            costs.soft.frontload += weightings.frontload.penalty;

    // excessiveDuration
    if(instance.exams[ex].duration > instance.periods[nextP].duration)
        costs.hard.excessiveDuration++;

    // comparing assignement before and after

    // room exclusive
    std::list<ExamId> const& meeting = examsByPeriodRoom[nextP][nextR];

    if(meeting.size() == 1 && instance.examIsRoomExclusive[*meeting.begin()])
        costs.hard.roomConstraint++;
    if(meeting.size() > 0 && instance.examIsRoomExclusive[ex])
        costs.hard.roomConstraint++;

    // overcapacity
    {
        int mySize = instance.exams[ex].students.size();

        // add
        int capMeeting = instance.rooms[nextR].capacity;

        int meetingSum = 0; // = sum(instance.exams[j].students.size() for j in meeting)

        std::list<ExamId> const& meeting = examsByPeriodRoom[nextP][nextR];

        for(ExamId j : meeting)
            meetingSum += instance.exams[j].students.size();

        if(meetingSum > capMeeting) // already over capacity
            costs.hard.overCapacity += mySize;
        else if(meetingSum + mySize > capMeeting)
            costs.hard.overCapacity += meetingSum + mySize - capMeeting;
        else
            ; // no overcapacity created
    }

    // mixedDurations
    {
        int myDuration = instance.exams[ex].duration;

        if(USE_COLOR_STRUCTURE) {
            // no heap allocated, incremental info

            Color myDurationColor = instance.exams[ex].durationColor;
            auto const& meetingDuration = durationColorUsed[nextP][nextR];

            const_cast<MapVec<Color,int>&>(meetingDuration)[myDurationColor]++;
            int countMeeting = 0;
            for(int x : meetingDuration)
                if(x > 0)
                    countMeeting++;

            if(meetingDuration[myDurationColor] == 1 && countMeeting > 1)
                costs.soft.mixedDuration += weightings.nonMixedDurations;
            const_cast<MapVec<Color,int>&>(meetingDuration)[myDurationColor]--;
        } else {
            int NMeetingWithoutMe, NMeetingWithMe;

            if(1) {
                // set implementation
                std::set<Minutes> durationsMeeting;

                for(ExamId j : meeting)
                    durationsMeeting.insert(instance.exams[j].duration);

                NMeetingWithoutMe = durationsMeeting.size();
                NMeetingWithMe = durationsMeeting.size() + (durationsMeeting.count(myDuration) ? 0 : 1);
            } else {
                // N² implementation, no heap allocated
                throw std::invalid_argument("not implemented");
            }

            if(NMeetingWithoutMe > 0)
                costs.soft.mixedDuration += (NMeetingWithMe - NMeetingWithoutMe) * weightings.nonMixedDurations;
        }
    }

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

void MoveNeighborhood::iterate(Solution *rawStep, std::function<void ()> yield) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;
    int E = instance.E(), P = instance.P(), R = instance.R();

    for(exam = 0; exam < E; exam++) {
        bperiod = sol->periods[exam];
        broom = sol->periods[exam];
        for(period = 0; period < P; period++) {
            for(room = 0; room < R; room++) {
                sol->move(instance, exam, period, room);
                yield();
                sol->move(instance, exam, bperiod, broom);
            }
        }
    }
}

Solution *MoveNeighborhood::computeStep(Solution *rawStep) {
    emili::iteration_increment();
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;

    int E = instance.E();

    if(exam == E)
        return nullptr;

    bperiod = sol->periods[exam];
    broom = sol->rooms[exam];

    sol->move(instance, exam, period, room);

    return sol;
}

void MoveNeighborhood::reverseLastMove(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;

    int R = instance.rooms.size();
    int P = instance.periods.size();

    sol->move(instance, exam, bperiod, broom);

    if(++room == R) {
        room = 0;
        if(++period == P) {
            period = 0;
            ++exam;
        }
    }
}

MoveNeighborhood::MoveNeighborhood(const ExamTT &instance_)
    : instance(instance_)
{
    reset();
}

Neighborhood::NeighborhoodIterator MoveNeighborhood::begin(Solution *base) {
    return emili::Neighborhood::begin(base);
}

Solution *MoveNeighborhood::step(Solution *currentSolution) {
    return this->computeStep(currentSolution);
}

void MoveNeighborhood::reset() {
    exam = 0;
    period = 0;
    room = 0;
}

int MoveNeighborhood::size() {
    int E = instance.exams.size();
    int R = instance.rooms.size();
    int P = instance.periods.size();

    return E * P * R;
}

Solution *MoveNeighborhood::random(Solution *currentSolution) {
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution->clone();

    int E = instance.exams.size(),
        P = instance.periods.size(),
        R = instance.rooms.size();

    Random r;
    sol->move(instance, r.randrange(E), r.randrange(P), r.randrange(R));
    return sol;
}

void MoveNeighborhood::randomStep(Solution *currentSolution)
{
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution;

    int E = instance.exams.size(),
        P = instance.periods.size(),
        R = instance.rooms.size();

    Random r;
    re = r.randrange(E);
    rp = sol->periods[re];
    rr = sol->rooms[re];
    sol->move(instance, re, r.randrange(P), r.randrange(R));
}

void MoveNeighborhood::reverseLastRandomStep(Solution *currentSolution)
{
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution;
    sol->move(instance, re, rp, rr);
}

Solution *RandomInitialSolution::generateSolution()
{
    ExamTT const& instance = (ExamTT const&) this->instance;

    ExamTTSolution* sol = new ExamTTSolution;

    sol->initRandom(instance, random);

    cout << "Initial solution " << sol->getSolutionRepresentation() << endl;

    return sol;
}

Solution *RandomInitialSolution::generateEmptySolution()
{
    return new ExamTTSolution;
}

Solution *ConstructorInitialSolution::generateSolution()
{
    ExamTT const& instance = (ExamTT const&) this->instance;

    auto sol = constructor->constructFull();

    cout << "Initial solution " << sol->getSolutionRepresentation() << endl;

    return sol;
}

Solution *ConstructorInitialSolution::generateEmptySolution()
{
    return new ExamTTSolution;
}

Solution *ZeroInitialSolution::generateSolution()
{
    ExamTT const& instance = (ExamTT const&) this->instance;

    ExamTTSolution* sol = new ExamTTSolution;
    sol->periods.assign(instance.exams.size(), 0);
    sol->rooms.assign(instance.exams.size(), 0);

    for(size_t i = 0; i < firstAssign.size(); ++i)
        sol->assignement(i) = firstAssign[i];

    sol->buildStructures(instance);
    sol->computeCost(instance);

    cout << "Initial solution " << sol->getSolutionRepresentation() << endl;

    return sol;
}

Solution *ZeroInitialSolution::generateEmptySolution()
{
    return new ExamTTSolution;
}

void SwapNeighborhood::iterate(Solution *rawStep, std::function<void ()> yield) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;
    int E = instance.students.size();

    for(e1 = 0; e1 < E; e1++) {
        for(e2 = e1 + 1; e2 < E; e2++) {
            sol->swap(instance, e1, e2);
            yield();
            sol->swap(instance, e1, e2);
        }
    }
}

Solution *SwapNeighborhood::computeStep(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;
    int E = instance.students.size();

    if(e1 == E)
        return nullptr;

    sol->swap(instance, e1, e2);

    return sol;
}

void SwapNeighborhood::reverseLastMove(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;
    int E = instance.students.size();

    sol->swap(instance, e1, e2);

    if(++e2 == E)
        e2 = ++e1 + 1;
}

SwapNeighborhood::SwapNeighborhood(const ExamTT &instance_)
    : instance(instance_)
{
    reset();
}

Solution *SwapNeighborhood::step(Solution *currentSolution)
{
    return computeStep(currentSolution);
}

void SwapNeighborhood::reset() {
    e1 = 0;
    e2 = 0;
}

int SwapNeighborhood::size() {
    int E = instance.students.size();
    return E * (E-1) / 2;
}

Solution *SwapNeighborhood::random(Solution *currentSolution) {
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution->clone();

    int E = instance.exams.size();

    Random r;
    sol->swap(instance, r.randrange(E), r.randrange(E));
    return sol;
}

void SwapNeighborhood::randomStep(Solution *currentSolution)
{
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution;
    int E = instance.exams.size();

    Random r;
    re1 = r.randrange(E);
    re2 = r.randrange(E);
    sol->swap(instance, re1, re2);
}

void SwapNeighborhood::reverseLastRandomStep(Solution* currentSolution)
{
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution;
    sol->swap(instance, re1, re2);
}

KempeChainNeighborhood::KempeChainNeighborhood(const ExamTT &instance_)
    : instance(instance_)
{
    reset();
}

void KempeChainNeighborhood::iterate(Solution* rawStep, std::function<void ()> yield) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;
    int E = instance.exams.size(), R = instance.rooms.size(), P = instance.periods.size();

    for(exam = 0; exam < E; exam++) {
        t0 = sol->periods[exam];
        for(t1 = 0; t1 < P; t1++) {
            if(t1 != sol->periods[exam]) {
                chain.clear();
                createChain(sol, exam);
                swapPeriodsInChain(sol);
                yield();
                swapPeriodsInChain(sol);
            }
        }
    }
}

void KempeChainNeighborhood::reset() {
    exam = 0;
    t1 = -1; // so that t1+1 = 0
    chain.clear();
}

int KempeChainNeighborhood::size() {
    int E = instance.exams.size();
    int P = instance.periods.size();

    return E * (P-1); // at most
}

Solution* KempeChainNeighborhood::computeStep(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;
    int E = instance.exams.size(), R = instance.rooms.size(), P = instance.periods.size();

    ++t1;
    if(t1 == sol->periods[exam])
        ++t1;

    if(t1 == P) {
        ++exam;
        t1 = sol->periods[exam] == 0 ? 1 : 0;
    }

    if(exam == instance.E())
        return nullptr;

    t0 = sol->periods[exam];

    chain.clear();
    createChain(sol, exam);
    swapPeriodsInChain(sol);

    return sol;
}

void KempeChainNeighborhood::swapPeriodsInChain(ExamTTSolution* sol) {
    for(ExamId i : chain)
        sol->movePeriod(instance, i, t0 == sol->periods[i] ? t1 : t0);
}

void KempeChainNeighborhoodFastIter::swapPeriodsInChainList(ExamTTSolution* sol) {
    for(ExamId i : chainList)
        sol->movePeriod(instance, i, t0 == sol->periods[i] ? t1 : t0);
}

void KempeChainNeighborhood::createChain(ExamTTSolution* sol, ExamId x) {
    if(chain.count(x))
        return;

    if(sol->periods[x] == t0 || sol->periods[x] == t1) {
        chain.insert(x);
        for(ExamId j : instance.hasStudentsInCommonOfExam(x))
            createChain(sol, j);
    }
}

void KempeChainNeighborhood::reverseLastMove(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;
    swapPeriodsInChain(sol);
}

Solution *KempeChainNeighborhood::step(Solution *currentSolution) {
    return computeStep(currentSolution);
}

Solution* KempeChainNeighborhood::random(Solution *currentSolution) {
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution->clone();
    int E = instance.E(), P = instance.P();

    Random r;
    exam = r.randrange(E);
    t0 = sol->periods[exam];
    t1 = r.randrange(P-1);
    t1 += t1 >= t0;

    chain.clear();
    createChain(sol, exam);
    swapPeriodsInChain(sol);
    return sol;
}

void KempeChainNeighborhood::randomStep(Solution *currentSolution) {
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution;
    int E = instance.E(), P = instance.P();

    Random r;
    exam = r.randrange(E);
    t0 = sol->periods[exam];
    t1 = r.randrange(P-1);
    t1 += t1 >= t0;

    chain.clear();
    createChain(sol, exam);
    swapPeriodsInChain(sol);
}

void KempeChainNeighborhood::reverseLastRandomStep(Solution *currentSolution) {
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution;

    swapPeriodsInChain(sol);
}

void stats::kempe_print_iteration(ExamTT& instance, std::vector<int> initPeriods, bool useFastIter, bool useIterate) {
    emili::Neighborhood* neigh = useFastIter ? (emili::Neighborhood*)new KempeChainNeighborhoodFastIter(instance) : (emili::Neighborhood*)new KempeChainNeighborhood(instance);
    ExamTTSolution* sol = new ExamTTSolution;

    sol->initFromPeriodsAndZeroRoom(instance, initPeriods);

    cout << sol->periods << endl;
    if(useIterate) {
        neigh->iterate(sol, [sol](){
            cout << sol->periods << endl;
        });
    } else {
        for(Solution* a : neigh->stdIterate(sol)) {
            ExamTTSolution* c = dynamic_cast<ExamTTSolution*>(a);
            cout << c->periods << endl;
        }
    }

    delete neigh;
    delete sol;
}

void stats::kempe_compare_size_fast_iter(ExamTT &instance) {
    ExamTTSolution sol;
    sol.initRandom(instance);
    KempeChainNeighborhood nx(instance);
    KempeChainNeighborhoodFastIter ny(instance);
    int a = 0, b = 0, c = 0, d = 0;
    int N = 50;
    for(int i = 0; i < N; i++) {
        a = 0, c = 0;
        nx.reset();
        nx.iterate(&sol, [&a]{
            a++;
        });
        nx.reset();
        for(auto s : nx.stdIterate(&sol))
            c++;
    }

    cout <<
        "  Normal (iterate, computeStep): " << a << "," << c << endl
    ;

    for(int i = 0; i < N; i++) {
        b = 0, d = 0;
        ny.reset();
        ny.iterate(&sol, [&b]{
            b++;
        });
        ny.reset();
        for(auto s : ny.stdIterate(&sol))
            d++;
    }
    cout <<
        "FastIter (iterate, computeStep): " << b << "," << d << endl
    ;
}

void test::kempe_iteration_vs_random(ExamTT& instance, std::vector<int> initPeriods, int N, bool useFastIter, bool useIterate) {
    ExamTTSolution* sol = new ExamTTSolution;
    sol->initFromPeriodsAndZeroRoom(instance, initPeriods);
    kempe_iteration_vs_random(instance, *sol, N, useFastIter, useIterate);
}

void test::kempe_iteration_vs_random(ExamTT& instance, ExamTTSolution& s, int N, bool useFastIter, bool useIterate) {

    ExamTTSolution* sol = &s;

    emili::Neighborhood* neigh = useFastIter ? (emili::Neighborhood*)new KempeChainNeighborhoodFastIter(instance) : (emili::Neighborhood*)new KempeChainNeighborhood(instance);

    std::set<std::vector<int>> iteratedNeighborhood;

    if(useIterate) {
        neigh->iterate(sol, [sol, &iteratedNeighborhood]{
            iteratedNeighborhood.insert(sol->periods);
        });
    } else {
        for(Solution* a : neigh->stdIterate(sol)) {
            ExamTTSolution* c = dynamic_cast<ExamTTSolution*>(a);
            iteratedNeighborhood.insert(c->periods);
        }
    }

    delete neigh;
    neigh = useFastIter ? (emili::Neighborhood*) new KempeChainNeighborhoodFastIter(instance) : (emili::Neighborhood*) new KempeChainNeighborhood(instance);

    auto costBefore = sol->costs;

    for(int i = 0; i < N; i++) {
        neigh->randomStep(sol);

        if(! iteratedNeighborhood.count(sol->periods)) {
            ostringstream oss;
            oss << "random neighborhood not in iterated neighborhood: ";
            if(sol->periods.size() < 20)
                oss << "&periods=" << sol->periods;
            if(iteratedNeighborhood.size() < 100)
                oss << "&neighb=" << iteratedNeighborhood;
            throw std::invalid_argument(oss.str());
        }

        neigh->reverseLastRandomStep(sol);

        if(costBefore != sol->costs) {
            throw std::invalid_argument("reverse last random step not reversing the step");
        }
    }
}

MixedMoveSwapNeighborhood::MixedMoveSwapNeighborhood(ExamTT const& instance_, double swapRate_)
    : move(instance_), swap(instance_), swapRate(swapRate_)
{

}

Solution* MixedMoveSwapNeighborhood::random(Solution *currentSolution) {
    return emili::generateRealRandomNumber() <= swapRate ? swap.random(currentSolution) : move.random(currentSolution);
}

void MixedMoveSwapNeighborhood::randomStep(Solution *currentSolution)
{
    chosenSwap = emili::generateRealRandomNumber() <= swapRate;
    if(chosenSwap)
        swap.randomStep(currentSolution);
    else
        move.randomStep(currentSolution);
}

void MixedMoveSwapNeighborhood::reverseLastRandomStep(Solution *currentSolution)
{
    if(chosenSwap)
        swap.reverseLastRandomStep(currentSolution);
    else
        move.reverseLastRandomStep(currentSolution);
}

MixedRandomNeighborhood::MixedRandomNeighborhood(std::vector<Neighborhood *> ns, std::vector<int> weights) : neighborhoods(ns) {
    if(weights.size() != ns.size())
        throw std::invalid_argument("weights.size != neighborhoods.size");
    if(weights.size() == 0)
        throw std::invalid_argument("weights.size = 0");
    cumul.resize(weights.size());

    cumul[0] = weights[0];
    for(int i = 1; i < weights.size(); i++)
        cumul[i] = cumul[i-1] + weights[i];

    _size = 0;
    for(auto n : neighborhoods)
        _size += n->size();
}

Solution* MixedRandomNeighborhood::random(Solution *currentSolution) {
    auto x = emili::generateRandRange(cumul.back());
    i = 0;
    while(x >= cumul[i])
        i++;
    return neighborhoods[i]->random(currentSolution);
}

void MixedRandomNeighborhood::randomStep(Solution *currentSolution)
{
    auto x = emili::generateRandRange(cumul.back());
    i = 0;
    while(x >= cumul[i])
        i++;
    neighborhoods[i]->randomStep(currentSolution);
}

void MixedRandomNeighborhood::reverseLastRandomStep(Solution *currentSolution)
{
    neighborhoods[i]->reverseLastRandomStep(currentSolution);
}

MixedRandomNeighborhoodProba::MixedRandomNeighborhoodProba(std::vector<Neighborhood *> ns, std::vector<float> weights) : neighborhoods(ns) {
    if(weights.size() != ns.size() - 1)
        throw std::invalid_argument("weights.size != neighborhoods.size - 1");
    if(ns.size() == 0)
        throw std::invalid_argument("n.size = 0");
    cumul.resize(ns.size());

    if(ns.size() > 1) {
        cumul[0] = weights[0];
        for(int i = 1; i < weights.size(); i++)
            cumul[i] = cumul[i-1] + weights[i];
    }
    cumul.back() = 1.1; // in theory it's 1 but in case generateRandomNumber returns 1...

    _size = 0;
    for(auto n : neighborhoods)
        _size += n->size();
}

Solution* MixedRandomNeighborhoodProba::random(Solution *currentSolution) {
    auto x = emili::generateRealRandomNumber();
    i = 0;
    while(x >= cumul[i])
        i++;
    return neighborhoods[i]->random(currentSolution);
}

void MixedRandomNeighborhoodProba::randomStep(Solution *currentSolution) {
    auto x = emili::generateRealRandomNumber();
    i = 0;
    while(x >= cumul[i])
        i++;
    neighborhoods[i]->randomStep(currentSolution);
}

void MixedRandomNeighborhoodProba::reverseLastRandomStep(Solution *currentSolution) {
    neighborhoods[i]->reverseLastRandomStep(currentSolution);
}

BruteForce::BruteForce(Instance &i)
    : BruteForce(i, {})
{

}

BruteForce::BruteForce(Instance &i, std::vector<std::pair<PeriodId, RoomId> > firstAssign_)
    : LocalSearch(
          * new ZeroInitialSolution(i, firstAssign_),
          * new WhileTrueTermination(),
          * new MoveNeighborhood(instance)
    )
    , instance(i)
    , firstAssign(firstAssign_)
{

}

Solution *BruteForce::search(Solution *initial)
{
    constexpr bool DEBUG = true;
    constexpr bool SEE_ALL_SOL = false;

    int E = instance.exams.size();
    int P = instance.periods.size();
    int R = instance.rooms.size();

    ExamTTSolution* sol = (ExamTTSolution*) initial;

    bestSoFar = sol->clone();

    std::function<void (int)> fNoDebug = [this, sol, E, P, R, &fNoDebug](int i){
        if(i >= E) {
            if(sol->costs < ((ExamTTSolution*)bestSoFar)->costs)
                *bestSoFar = *sol;
        } else {
            for(int p = 0; p < P; p++) {
                for(int r = 0; r < R; r++) {
                    sol->move(instance, i, p, r);
                    fNoDebug(i + 1);
                }
            }
        }
    };

    std::function<void (int)> fDebug = [this, sol, E, P, R, &fDebug, &fNoDebug](int i){
        if(i >= E) {
            if(SEE_ALL_SOL)
                std::cout << sol->getSolutionValue() << "\t" << sol->getSolutionRepresentation() << endl;
            if(sol->costs < ((ExamTTSolution*)bestSoFar)->costs)
                *bestSoFar = *sol;
        } else {
            auto nextF = SEE_ALL_SOL ? fDebug : i == firstAssign.size() + sizeOfDebug-1 ? fNoDebug : fDebug;
            for(int p = 0; p < P; p++) {
                for(int r = 0; r < R; r++) {
                    if(! SEE_ALL_SOL) {
                        if(firstAssign.size() <= i && i < firstAssign.size() + sizeOfDebug)
                            std::cout << std::string(i, ' ') << i << " " << p << '/' << P << " " << r << '/' << R << "; best_cost=" << bestSoFar->getSolutionValue() << "; cur=" << sol->getSolutionRepresentation() << endl;
                    }
                    sol->move(instance, i, p, r);

                    nextF(i + 1);
                }
            }
        }
    };

    if(DEBUG)
        fDebug(firstAssign.size());
    else
        fNoDebug(firstAssign.size());

    return bestSoFar;
}

FixedRandomDestructor::FixedRandomDestructor(ExamTT const& instance_, const int G)
    : instance(instance_)
    , inserted(G)
{
    if(!(G <= instance.E()))
        throw std::invalid_argument("assert G <= instance.E()");
}

Solution *FixedRandomDestructor::destruct(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;

    Random ran;

    int E = instance.E();

    for(size_t i = 0; i < inserted.size(); i++) {
        int e = ran.randrange(E);
        while(! sol->isAssigned(e))
            e = ran.randrange(E);

        inserted[i] = e;
        sol->removeExam(instance, e);
    }

    return sol;
}

Solution* IteratedGreedyNeighborhood::computeStep(Solution *step){
    throw std::invalid_argument("not implemented");
}

void IteratedGreedyNeighborhood::reverseLastMove(Solution *step){
    throw std::invalid_argument("not implemented");
}

void IteratedGreedyNeighborhood::reset() {
    inserted.resize(G);
}

Solution *IteratedGreedyNeighborhood::random(Solution *currentSolution) {
    auto s = currentSolution->clone();
    randomStep(s);
    return s;
}

void IteratedGreedyNeighborhood::randomStep(Solution *currentSolution) {
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution;

    Random ran;

    int E = instance.E();

    for(int i = 0; i < G; i++) {
        int e = ran.randrange(E);
        while(! sol->isAssigned(e))
            e = ran.randrange(E);

        inserted[i] = e;
        sol->removeExam(instance, e);
    }

    // greedy construct, random order
    ran.shuffle(inserted);

    for(ExamId e : inserted)
        sol->addExam(instance, e, BestInsertHeuristic(instance).searchPosition(sol, e)); // greedy insert
}

void IteratedGreedyNeighborhood::reverseLastRandomStep(Solution *currentSolution) {
    throw std::invalid_argument("not implemented");
}

namespace {
template <typename T>
inline T pop_from_set(std::set<T> & S) {
    auto it = S.begin();
    T r = *it;
    S.erase(it);
    return r;
}
}

/**
 * Implement iteration as :
    for i in range(P)
        for j in range(i+1,P)
            chains = components(C[i] | C[j])
            for chain in chains:
                apply(chain, i, j)
                reverse(chain, i, j)
                ...

 * Implementation :
    examsByPeriod = [set() for p in range(P)]
    for e in range(E)
        examsByPeriod[periods[e]].insert(e)

    for i in range(P)
        for j in range(i+1,P)
            A,B = examsByPeriod[i], examsByPeriod[j] // copy
            while A and B
                e = A.pop() if A else B.pop()
                chain = createChainDestruct(e,A,B) # elements go from A and B to chain
                apply(chain, i, j)
                reverse(chain, i, j)
                ...
 */
void KempeChainNeighborhoodFastIter::iterate(Solution *base, std::function<void ()> yield) {
    ExamTTSolution* sol = (ExamTTSolution*) base;

    int E = instance.E(), P = instance.P();

    examsByPeriod.assign(P, std::set<ExamId>());
    for(ExamId e = 0; e < E; e++)
        examsByPeriod[sol->periods[e]].insert(e);

    for(t0 = 0; t0 < P; t0++) {
        for(t1 = t0+1; t1 < P; t1++) {
            A = examsByPeriod[t0]; // copy
            A.insert(examsByPeriod[t1].begin(), examsByPeriod[t1].end());
            while(A.size()) {
                chainList.clear();
                createChainDestruct(*A.begin());
                swapPeriodsInChainList(sol);
                yield(); // may throw an exception
                swapPeriodsInChainList(sol);
            }
        }
    }
}

void KempeChainNeighborhoodFastIter::reset() {
    t1 = -1;
    t0 = 0;
    chainList.clear();
    A.clear();
}

Solution *KempeChainNeighborhoodFastIter::computeStep(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;

    int E = instance.E(), P = instance.P(), R = instance.R();

    if(examsByPeriod.empty()) {
        examsByPeriod.resize(P); // assign each color to a set of exam of that color
        for(ExamId e = 0; e < E; e++)
            examsByPeriod[sol->periods[e]].insert(e);
    }

    if(A.empty()) {
        for(;;) {
            if(t1 == -1) {
                t1 = 1;
                t0 = 0;
            }
            else if(++t1 == P) {
                if(++t0 >= P-1)
                    return nullptr;
                t1 = t0 + 1;
            }

            A = examsByPeriod[t0]; // copy
            A.insert(examsByPeriod[t1].begin(), examsByPeriod[t1].end());

            if(A.size())
                break;
        }
    }

    chainList.clear();
    createChainDestruct(*A.begin());

    swapPeriodsInChainList(sol);

    return sol;
}

void KempeChainNeighborhoodFastIter::createChainDestruct(ExamId e) {
    auto itA = A.find(e);
    if(itA != A.end()) {
        A.erase(itA);
        chainList.push_back(e);
        for(ExamId j : instance.hasStudentsInCommonOfExam(e))
            createChainDestruct(j);
    }
}

void KempeChainNeighborhoodFastIter::reverseLastMove(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;

    swapPeriodsInChainList(sol);
}

RandomOrderInserter::RandomOrderInserter(const Instance &instance_, InsertHeuristic *insertHeuristic_)
    : BaseConstructor(instance_)
    , insertHeuristic(insertHeuristic_)
{

}

Solution *RandomOrderInserter::construct(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;

    // instead of doing a copy, can re use std::vector from destructor
    inserted.assign(sol->unAssignedExamList.begin(), sol->unAssignedExamList.end());

    Random().shuffle(inserted);

    for(ExamId e : inserted)
        sol->addExam(instance, e, insertHeuristic->searchPosition(sol, e));

    return sol;
}

Solution *DegreeInserter::construct(Solution *raw) {
    ExamTTSolution* sol = (ExamTTSolution*) raw;

    data.resize(sol->unAssignedExamList.size());
    int i = 0;
    for(ExamId e : sol->unAssignedExamList)
        data[i++] = {instance.hasStudentsInCommonOfExam(e).size(), e};
    std::sort(data.begin(), data.end(), std::greater<Two<int>>());

    for(auto p : data) {
        ExamId e = p.second;
        sol->addExam(instance, e, insertHeuristic->searchPosition(sol, e));
    }

    return sol;
}

Assignement BestInsertHeuristic::searchPosition(ExamTTSolution * sol, ExamId e) {
    int P = instance.P(), R = instance.R();

    auto before = sol->getSolutionValue();

    int bestP = 0, bestR = 0;
    double bestDelta = 1e9;

    for(int p = 0; p < P; p++)
    for(int r = 0; r < R; r++) {
        sol->addExam(instance, e, p, r);

        auto delta = sol->getSolutionValue() - before;
        if(delta < bestDelta) {
            bestDelta = delta;
            bestP = p;
            bestR = r;
        }

        sol->removeExam(instance, e);
    }
    return {bestP, bestR};
}



NRandomDaysDestructor::NRandomDaysDestructor(const Instance &instance_, const int N_)
    : instance(instance_)
    , N(N_)
{
    if(!(N <= instance.Days()))
        throw std::invalid_argument("assert N <= instance.Days()");
    days.resize(instance.Days());
    for(size_t i = 0; i < days.size(); i++)
        days[i] = i;
}

Solution* NRandomDaysDestructor::destruct(Solution *solution) {
    ExamTTSolution* sol = (ExamTTSolution*) solution;

    Random().shuffleFront(days, N);

    for(int i = 0; i < N; i++) {
        Two<int> day = instance.daysToPeriod[days[i]];
        for(ExamId e = 0; e < instance.E(); e++)
            if(sol->isAssigned(e))
                if(day.first <= sol->periods[e] && sol->periods[e] < day.second)
                    sol->removeExam(instance, e);
    }

    return solution;
}

NGroupedDaysDestructor::NGroupedDaysDestructor(const Instance &instance_, const int N_)
    : instance(instance_)
    , N(N_)
{
    if(!(N <= instance.Days()))
        throw std::invalid_argument("assert N <= instance.Days()");
}

Solution* NGroupedDaysDestructor::destruct(Solution *solution) {
    ExamTTSolution* sol = (ExamTTSolution*) solution;

    int i = Random().randrange(instance.Days() - N);
    int j = i+N;
    int di = instance.daysToPeriod[i].first;
    int dj = instance.daysToPeriod[j].second;

    for(ExamId e = 0; e < instance.E(); e++)
        if(sol->isAssigned(e))
            if(di <= sol->periods[e] && sol->periods[e] < dj)
                sol->removeExam(instance, e);

    return solution;
}

NGroupedPeriodsDestructor::NGroupedPeriodsDestructor(const Instance &instance_, const int N_)
    : instance(instance_)
    , N(N_)
{
    if(!(N <= instance.P()))
        throw std::invalid_argument("assert N <= instance.P()");
}

Solution* NGroupedPeriodsDestructor::destruct(Solution *solution) {
    ExamTTSolution* sol = (ExamTTSolution*) solution;

    int i = Random().randrange(instance.P() - N);
    int j = i+N;

    for(ExamId e = 0; e < instance.E(); e++)
        if(sol->isAssigned(e))
            if(i <= sol->periods[e] && sol->periods[e] < j)
                sol->removeExam(instance, e);

    return solution;
}

Solution *DSaturInserter::construct(Solution *partial) {
    ExamTTSolution* sol = (ExamTTSolution*) partial;

    // dsat(v) = number of different "colors" in neighborhood of v = saturation
    throw std::invalid_argument("not imp");

    return sol;
}


}
}
