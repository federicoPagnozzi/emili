#include "examtt.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <functional>
#include <tuple>
#include <memory>

using namespace std;

namespace emili
{
namespace ExamTT
{

int ExamTTSolution::numberOfClones = 0;

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

const Prefix EPrefix{"e", 5};
const Prefix PPrefix{"p", 4};
const Prefix RPrefix{"r", 4};

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

    // will be featured based, now it's just a constant
    hardWeight = 1000.0;
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

    // EXCLUSIVE, EXAM_COINCIDENCE
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

    constexpr bool USE_DELTA = true;

    if(USE_DELTA) {
        updateMove(instance, e, nextP, nextR, this->costs);
        setSolutionValue(costs.total(instance.hardWeight));
    }

    examsByPeriods[prevP].erase(examsByPeriodsIterators[e]);
    examsByRooms[prevR].erase(examsByRoomsIterators[e]);

    examsByPeriodsIterators[e] = examsByPeriods[nextP].insert(e).first;
    examsByRoomsIterators[e] = examsByRooms[nextR].insert(e).first;

    periods[e] = nextP;
    rooms[e] = nextR;

    if(! USE_DELTA) {
        computeCost(instance);
    }
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
    setSolutionValue(costs.total(instance.hardWeight));
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

void ExamTT::presentation(ostream & log) {
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

    log << "* ROOM_EX Exams that must be room exclusive (" << inst.examsRoomsExclusive.size() << ")" << endl;
    for(ExamId e : inst.examsRoomsExclusive)
        log << e << endl;

    log << "* InstitutionalWeightings" << endl
         << inst.institutionalWeightings << endl;
}

void ExamTT::testDelta(ExamTTSolution& sol, std::ostream& log) {
    auto& inst = *this;
    constexpr bool CHECK_EACH_MOVE = true;

    int E = inst.exams.size(),
        P = inst.periods.size(),
        R = inst.rooms.size();

    Random ran;
    ran.seed(89);
    for(int i = 0; i < 1000; i++) {
        int e = ran.randrange(E), p = ran.randrange(P), r = ran.randrange(R);
        int bp = sol.periods[e], br = sol.rooms[e];
        sol.move(inst, e,p,r);

        if(CHECK_EACH_MOVE) {
            if(! sol.costs.exactlyEqual(sol.computeAndGetCost(inst))) {
                log << "Error " << "Real " << sol.computeAndGetCost(inst).print(inst)
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

void interactiveSolution(ExamTTSolution& sol, ExamTT& inst) {
    auto& instance = inst;

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
             << "Quit                        : quit" << endl
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

        } else if(isQuit){
            break;
        } else {
            cout << "Unknown command '" << name << "'" << endl;
        }
    }
}

void test() {
    // string filename = "../simple1.exam"; // "../itc2007-exam-instances/exam_comp_set3.exam"; // args[0];
    // string filename = "../itc2007-exam-instances/exam_comp_set3.exam";
    string filename = "../examtt-instances/Instances/art000.exam";

    ExamTT inst;
    ExamTT& instance = inst;
    InstanceParser parser(filename);

    if(! parser.file) {
        cerr << "NO FILE" << endl;
        return;
    }

    parser.parse(inst);

    ofstream log("log");
    inst.presentation(log);

    ExamTTSolution sol;
    Random ran(89);
    sol.initRandom(inst, ran);

    log << endl << "* Solution" << endl;
    sol.printTo(inst, log);

    inst.testDelta(sol, log);

    log.close();

    interactiveSolution(sol, inst);
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

    costs = other->costs;
    periods = other->periods;
    rooms = other->rooms;

    examsByPeriodsIterators.resize(periods.size());
    examsByRoomsIterators.resize(rooms.size());

    // Using inner copy

    examsByPeriods = other->examsByPeriods;
    examsByRooms = other->examsByRooms;

    for(set<ExamId>& S : examsByPeriods)
        for(auto it = S.begin() ; it != S.end(); ++it)
            examsByPeriodsIterators[*it] = it;

    for(set<ExamId>& S : examsByRooms)
        for(auto it = S.begin() ; it != S.end(); ++it)
            examsByRoomsIterators[*it] = it;
}

Solution* ExamTTSolution::clone() {
    numberOfClones++;

    ExamTTSolution* other = new ExamTTSolution(getSolutionValue());

    other->costs = costs;
    other->periods = periods;
    other->rooms = rooms;

    other->examsByPeriodsIterators.resize(periods.size());
    other->examsByRoomsIterators.resize(rooms.size());

    other->examsByPeriods = examsByPeriods;
    other->examsByRooms = examsByRooms;

    for(set<ExamId>& S : other->examsByPeriods)
        for(auto it = S.begin(); it != S.end(); ++it)
            other->examsByPeriodsIterators[*it] = it;

    for(set<ExamId>& S : other->examsByRooms)
        for(auto it = S.begin(); it != S.end(); ++it)
            other->examsByRoomsIterators[*it] = it;

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
    other->examsByPeriodsIterators.swap(examsByPeriodsIterators);
    other->examsByRoomsIterators.swap(examsByRoomsIterators);
    other->examsByPeriods.swap(examsByPeriods);
    other->examsByRooms.swap(examsByRooms);
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

    examsByPeriodsIterators.resize(E);
    examsByRoomsIterators.resize(E);
    examsByPeriods.assign(P, set<ExamId>());
    examsByRooms.assign(R, set<ExamId>());

    // O(n log n)

    for(ExamId e = 0; e < E; e++) {
        examsByPeriodsIterators[e] = examsByPeriods[periods[e]].insert(e).first;
        examsByRoomsIterators[e] = examsByRooms[rooms[e]].insert(e).first;
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

    InstitutionalWeightings const& weightings = instance.institutionalWeightings;

    // logging
    if(0) {
        cout << endl << "* Move " << EPrefix(ex) << " from "
             << PPrefix(prevP) << "," << RPrefix(prevR) << " to "
             << PPrefix(nextP) << "," << RPrefix(nextR) << endl;

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
    }

    Instance::ExamsRelated const& related = instance.examsRelated[ex];

    int s = nextP > prevP ? +1 :
            nextP < prevP ? -1 : 0;

    // if the period has changed
    if(s) {

        // period penalty
        costs.soft.periodsPenalty -= instance.periods[prevP].penalty;
        costs.soft.periodsPenalty += instance.periods[nextP].penalty;

        // after
        if(0) {
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

            for(ExamId j : related.afters)
                if(m <= periods[j] && periods[j] < M)
                    costs.hard.periodConstraintAfter -= s;

            for(ExamId j : related.befores)
                if(m < periods[j] && periods[j] <= M)
                    costs.hard.periodConstraintAfter += s;
        }

        // exclusion/coincidences

        //// remove
        for(ExamId j : related.exclusions)
            if(periods[j] == prevP)
                costs.hard.periodConstraintExclusion--;

        for(ExamId j : related.coincidences)
            if(periods[j] != prevP)
                costs.hard.periodConstraintCoincidence--;

        //// add
        for(ExamId j : related.exclusions)
            if(periods[j] == nextP)
                costs.hard.periodConstraintExclusion++;

        for(ExamId j : related.coincidences)
            if(periods[j] != nextP)
                costs.hard.periodConstraintCoincidence++;

        // for coincidence/exclusion, here is a different method : by counting intersections
        if(0){
            costs.hard.periodConstraintCoincidence += countIntersection(examsByPeriods[prevP], related.coincidences);
            costs.hard.periodConstraintCoincidence -= countIntersection(examsByPeriods[nextP], related.coincidences);

            costs.hard.periodConstraintExclusion -= countIntersection(examsByPeriods[prevP], related.exclusions);
            costs.hard.periodConstraintExclusion += countIntersection(examsByPeriods[nextP], related.exclusions);
        }

        // logging
        if(0) {
            cout << "Passing through" << endl;

            for(int i = prevP + s; (nextP - i) * s >= 0; i += s) {
                cout << setw(4) << i << " ";

                cout << (s == 1 ? " -" : " +") << countIntersection(examsByPeriods[i], related.afters)
                     << (s == 1 ? " +" : " -") << countIntersection(examsByPeriods[i], related.befores) << " ";

                cout << "AF" << intersection(examsByPeriods[i], related.afters);
                cout << "BE" << intersection(examsByPeriods[i], related.befores);
                cout << "CO" << intersection(examsByPeriods[i], related.coincidences);
                cout << "EX" << intersection(examsByPeriods[i], related.exclusions);
                cout << " -- " << examsByPeriods[i] << " ";
                cout << endl;

                // old way
                // costs.hard.periodConstraintAfter -= s * countIntersection(examsByPeriods[i], instance.afterOfExam[e]);
                // costs.hard.periodConstraintAfter += s * countIntersection(examsByPeriods[i], instance.beforeEqualOfExam[e]);
            }
        }
    }

    // room penalty
    costs.soft.roomsPenalty -= instance.rooms[prevR].penalty;
    costs.soft.roomsPenalty += instance.rooms[nextR].penalty;

    // about students in common : (twoExamsInARow, twoExamsInADay, periodSpread, simultaneousExams)

    for(ExamId j : related.hasStudentsInCommon) {
        Cost cost = instance.numberStudentsInCommon[ex][j];

        int relatedP = periods[j];

        // remove

        Two<Period const&> before = {instance.periods[relatedP], instance.periods[prevP]};
        if(prevP == relatedP) {
            costs.hard.simultaneousExams -= cost;
        } else {
            int prevDiff = abs(relatedP - prevP);

            if(before.first.date == before.second.date) {
                if(prevDiff == 1)
                    costs.soft.twoExamsInARow -= cost * weightings.twoInARow;
                else
                    costs.soft.twoExamsInADay -= cost * weightings.twoInADay;
            }

            if(prevDiff <= weightings.periodSpread)
                costs.soft.periodSpread -= cost;
        }

        // add

        Two<Period const&> after = {instance.periods[relatedP], instance.periods[nextP]};
        if(nextP == relatedP) {
            costs.hard.simultaneousExams += cost;
        } else {
            int nextDiff = abs(relatedP - nextP);

            if(after.first.date == after.second.date) {
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

    set<ExamId> leaving = intersection(examsByRooms[prevR], examsByPeriods[prevP]); // exams that I am leaving
    leaving.erase(leaving.find(ex)); // I am not included in the leaving

    set<ExamId> meeting = intersection(examsByRooms[nextR], examsByPeriods[nextP]); // exams that I will meet

    // room exclusive
    {
        if(leaving.size() == 1 && instance.examIsRoomExclusive[*leaving.begin()])
            costs.hard.roomConstraint--;
        if(leaving.size() > 0 && instance.examIsRoomExclusive[ex])
            costs.hard.roomConstraint--;

        if(meeting.size() == 1 && instance.examIsRoomExclusive[*meeting.begin()])
            costs.hard.roomConstraint++;
        if(meeting.size() > 0 && instance.examIsRoomExclusive[ex])
            costs.hard.roomConstraint++;
    }

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

        // add
        int capMeeting = instance.rooms[nextR].capacity;

        int meetingSum = 0; // = sum(instance.exams[j].students.size() for j in meeting)
        for(ExamId j : meeting)
            meetingSum += instance.exams[j].students.size();

        if(meetingSum > capMeeting) // already over capacity
            costs.hard.overCapacity += mySize;
        else if(meetingSum + mySize > capMeeting)
            costs.hard.overCapacity += meetingSum + mySize - capMeeting;
        else
            ; // no overcapacity created

        if(0) {
            cout << "Leaving " << leaving << " Meeting " << meeting << endl;
        }
    }

    // mixedDurations
    int myDuration = instance.exams[ex].duration;

    std::set<Minutes> durationsLeaving;
    std::set<Minutes> durationsMeeting;

    for(ExamId j : leaving)
        durationsLeaving.insert(instance.exams[j].duration);
    for(ExamId j : meeting)
        durationsMeeting.insert(instance.exams[j].duration);

    // using counters

    int NLeavingWithoutMe = durationsLeaving.size();
    int NMeetingWithoutMe = durationsMeeting.size();
    int NLeavingWithMe = durationsLeaving.size() + (durationsLeaving.count(myDuration) ? 0 : 1);
    int NMeetingWithMe = durationsMeeting.size() + (durationsMeeting.count(myDuration) ? 0 : 1);

    if(NLeavingWithoutMe > 0)
        costs.soft.mixedDuration += (NLeavingWithoutMe - NLeavingWithMe) * weightings.nonMixedDurations;
    if(NMeetingWithoutMe > 0)
        costs.soft.mixedDuration -= (NMeetingWithoutMe - NMeetingWithMe) * weightings.nonMixedDurations;

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

Solution *MoveNeighborhood::computeStep(Solution *rawStep) {
    emili::iteration_increment();
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;

    int E = instance.exams.size();

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
    rp = r.randrange(P);
    rr = r.randrange(R);
    sol->move(instance, re, rp, rr);
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
    : instance(instance)
{
    reset();
}

void KempeChainNeighborhood::reset() {
    exam = 0;
    t1 = -1; // so that +1 = 0
    chain.clear();
}

int KempeChainNeighborhood::size() {
    int E = instance.exams.size();
    int P = instance.periods.size();

    return E * (P-1); // at most
}

Solution* KempeChainNeighborhood::computeStep(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;

    int E = instance.exams.size();
    int R = instance.rooms.size();
    int P = instance.periods.size();

    ++t1;
    if(t1 == sol->periods[exam])
        ++t1;
    if(t1 == P)
        ++exam;

    if(exam == instance.E())
        return nullptr;

    t0 = sol->periods[exam];

    chain.clear();
    createChain(sol, exam);

    // currently have a lot of duplicates and we construct a chain per computeStep
    // should do like this :
    // for i in range(P)
    //   for j in range(i+1,P)
    //     chains = components(C[i] | C[j])
    //     for chain in chains:
    //       apply(chain, i, j)
    //       reverse(chain, i, j)

    for(ExamId i : chain)
        sol->movePeriod(instance, i, t0 == sol->periods[i] ? t1 : t0);
}

void KempeChainNeighborhood::createChain(ExamTTSolution* sol, ExamId x) {
    if(chain.count(x))
        return;

    if(sol->periods[x] == t0 || sol->periods[x] == t1)
        for(ExamId j : instance.hasStudentsInCommonOfExam(x))
            createChain(sol, j);
}

void KempeChainNeighborhood::reverseLastMove(Solution *rawStep) {
    ExamTTSolution* sol = (ExamTTSolution*) rawStep;

    for(ExamId i : chain)
        sol->movePeriod(instance, i, t0 == sol->periods[i] ? t1 : t0);
}

Solution *KempeChainNeighborhood::step(Solution *currentSolution) {
    return computeStep(currentSolution);
}

Solution* KempeChainNeighborhood::random(Solution *currentSolution) {
    ExamTTSolution* sol = (ExamTTSolution*) currentSolution->clone();

    int E = instance.exams.size();
    int P = instance.periods.size();

    // when is random called ?

    Random r;
    exam = r.randrange(E);
    t0 = sol->periods[exam];
    t1 = r.randrange(P-1);
    t1 += t1 >= t0;

    chain.clear();
    createChain(sol, exam);
    return sol;
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

Solution *BSUSA::search(Solution *initial) {
    auto sol = initial->clone();
    searchInPlace(sol);
    return sol;
}

void BSUSA::searchInPlace(Solution *initial) {
    auto best = initial->clone();

    // auto temp = ;
    neighbh->reset();

    auto accept = [](Solution* sol, double delta){
        return delta < 0;
    };

    int accepted = 0;
    int iterations = 0;

    do {
        // best = exploration->nextSolution();
        auto costBefore = initial->getSolutionValue();
        neighbh->randomStep(initial);
        auto delta = initial->getSolutionValue() - costBefore;

        if(accept(initial, delta)) {
            // initial is the same
            accepted++;
        } else {
            // we'll try another random
            neighbh->reverseLastRandomStep(initial);
        }

        // temp = update(temp)
        // accept.setTemp(temp)
    } while (true);

    best->swap(initial);
    delete best;
}

}
}
