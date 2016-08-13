#ifndef ExamTT_H
#include "../emilibase.h"
#define ExamTT_H

#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <list>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>
#include <functional>

// SA
#include "../SA/sa.h"

namespace emili
{
namespace ExamTT
{


// Output iterator that counts the number of times ++ is called
// can be casted to a int to get the value or get attribute .count
struct CountIterator {
    int count = 0; // contain the number of times ++ has been applied
    int dummyWrite; // contain last value written

    int&           operator*()          { return dummyWrite; }
    int const&     operator*() const    { return dummyWrite; }
    CountIterator& operator++()         { ++count; return *this; }
                   operator int() const { return count; }
};

/**
 * for(auto p : VectorExclusivePair(4))
 *   print(p.first, p.second) # (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
 */
template <typename T>
struct VectorExclusivePair {
    struct Iterator {
        T i, j;
        int N;

        std::pair<T,T> operator*(){
            return {i,j};
        }

        Iterator& operator++(){
            if(++j >= N)
                j = ++i + 1;
            return *this;
        }

        bool operator ==(Iterator const& a) const {
            return a.i == i && a.j == j;
        }

        bool operator !=(Iterator const& a) const {
            return !(a.i == i && a.j == j);
        }
    };

    int N;

    VectorExclusivePair(int N_) : N(N_) {}

    Iterator begin() {
        return {0,1,N};
    }

    Iterator end() {
        return {N,N+1,N};
    }
};

struct Random {

    Random() {

    }

    Random(long s) {
        seed(s);
    }

    void seed(long s) {
       emili::initializeRandom(s);
    }

    float random() {
        return emili::generateRealRandomNumber();
    }

    int randrange(int N) {
        return emili::generateRandRange(N);
    }

    int randrange(int a, int b) {
        return emili::generateRandRange(a,b);
    }

    int randint(int a, int b) {
        return emili::generateRandInt(a,b);
    }

    template <typename T>
    void shuffle(std::vector<T>& vec) {
        shuffleFront(vec, vec.size());
    }

    /**
     * <L> last elem of vec will be random
     */
    template <typename T>
    void shuffleEnd(std::vector<T>& vec, int L) {
        const int N = vec.size();
        for(int i = 0; i < L; i++)
            std::swap(vec[N - 1 - i], vec[randrange(N - i)]);
    }

    /**
     * <L> first elem of vec will be random
     */
    template <typename T>
    void shuffleFront(std::vector<T>& vec, int L) {
        const int N = vec.size();
        for(int i = 0; i < L; i++)
            std::swap(vec[i], vec[randrange(i, N)]);
    }
};

template <typename T>
int countIntersection(std::set<T> const& one, std::set<T> const& two) {
    return set_intersection(one.begin(), one.end(), two.begin(), two.end(), CountIterator()); // .count;

    int n = 0;
    auto a = one.begin(); auto A = one.end();
    auto b = two.begin(); auto B = two.end();
    if(a != A && b != B) {
        for(;;){
            if(*a < *b) {
                if(++a == A)
                    break;
            }
            else if(*b < *a) {
                if(++b == B)
                    break;
            }
            else {
                ++n, ++a, ++b;
                if(a == A || b == B)
                    break;
            }
        }
    }
    return n;
}

template <typename T, typename U>
void makeIntersection(std::set<T> const& one, std::set<T> const& two, U& res) {
    set_intersection(one.begin(), one.end(), two.begin(), two.end(), inserter(res, res.end()));
}

template <typename T, typename U>
void makeIntersectionBack(std::set<T> const& one, std::set<T> const& two, U& res) {
    set_intersection(one.begin(), one.end(), two.begin(), two.end(), std::back_inserter(res));
}

template <typename T>
void makeIntersectionWithReserve(std::set<T> const& one, std::set<T> const& two, std::vector<T>& res) {
    res.reserve(std::min(one.size(), two.size()));
    set_intersection(one.begin(), one.end(), two.begin(), two.end(), back_inserter(res));
}

template <typename T>
std::set<T> intersection(std::set<T> const& one, std::set<T> const& two) {
    std::set<T> res;
    set_intersection(one.begin(), one.end(), two.begin(), two.end(), inserter(res, res.end()));
    return res;
}

template <typename T>
std::vector<T> intersectionBack(std::set<T> const& one, std::set<T> const& two) {
    std::vector<T> res;
    makeIntersectionBack(one, two, res);
    return res;
}

template <typename T>
std::vector<T> intersectionWithReserve(std::set<T> const& one, std::set<T> const& two) {
    std::vector<T> res;
    makeIntersectionWithReserve(one, two, res);
    return res;
}

typedef std::tuple<int,int,int> Date; // y m d
typedef std::tuple<int,int,int> Time; // h m s

int seconds(Time const& time);

typedef int StudentId;
typedef int ExamId;
typedef int PeriodId;
typedef int RoomId;

typedef int Minutes;
typedef int Color;
typedef int DurationColor;
typedef int Day;
typedef int Cost; // A soft or hard cost
typedef int HCost; // Hard
typedef int SCost; // Soft

typedef int NumberOfStudents;
typedef int NumberOfPeriods;
typedef int NumberOfRooms;
typedef int NumberOfExams;

typedef std::pair<PeriodId, RoomId> Assignement;
typedef std::pair<PeriodId&, RoomId&> AssignementRef;
typedef std::pair<PeriodId const&, RoomId const&> AssignementConstRef;

class Exam {
public:
    Minutes duration;
    DurationColor durationColor;
    std::set<StudentId> students;

    static bool compareDuration(Exam const& e1, Exam const& e2) {
        return e1.duration < e2.duration;
    }

    static bool compareSizeStudents(Exam const& e1, Exam const& e2) {
        return e1.students.size() < e2.students.size();
    }
};

class Room {
public:
    NumberOfStudents capacity;
    SCost penalty;
};

class Period {
public:
    Date dateRaw;
    Day dateid; // day
    Time time;
    Minutes duration;
    SCost penalty;

    static bool compareDateTime(Period const& p1, Period const& p2) {
        return std::make_tuple(p1.dateid, p2.time) < std::make_tuple(p2.dateid, p2.time);
    }

    static bool compareDuration(Period const& p1, Period const& p2) {
        return p1.duration < p2.duration;
    }
};

struct MetaInfo {
    int Exams;
    int Periods;
    int Rooms;
    int Students;
    float ConflictDensity;
    int conflictsComponents;
    std::vector<float> conflictsComponentsDensities;
};

struct FrontLoadParams {
    NumberOfExams largestExams; // number of exams
    NumberOfPeriods time;
    SCost penalty;
};

std::ostream& operator <<(std::ostream& out, FrontLoadParams const& w);

struct InstitutionalWeightings {
    SCost twoInARow = 0,
          twoInADay = 0,
          nonMixedDurations = 0;
    NumberOfPeriods periodSpread = 0;
    FrontLoadParams frontload = {0,0,0};
};

std::ostream& operator <<(std::ostream& out, InstitutionalWeightings const& w);

template <typename T>
using Two = std::pair<T,T>;

template <typename T, typename U>
using MapVec = std::vector<U>;

class ExamTTSolution;

/* emili classes */
class ExamTT: public emili::Problem
{
public:
    ///////////////
    // Variables //
    ///////////////

    MapVec<ExamId, Exam> exams;
    MapVec<PeriodId, Period> periods;
    MapVec<RoomId, Room> rooms;

    // Hard constraints

    // Periods related
    std::vector<Two<ExamId>> examsCoincidence, // pair of exams that must be in same period
                             examsExclusion, // pair of exams that must not be in same period
                             examsAfter; // pair of exams where first.end <= second.begin

    // Rooms related
    std::set<ExamId> examsRoomsExclusive; // exams that must not share a room
    MapVec<ExamId, bool> examIsRoomExclusive;

    // Weights of institution
    InstitutionalWeightings institutionalWeightings;

    // used to compute scalar cost from (Hc,SC) cost
    double hardWeight = 0;

    // factor used at the end at the algorithm
    // irace will be based on this value
    // better integration of HC/SC problems should output the pair at the end of the program
    // then it's depend on the user (here irace) to use it's hardWeight
    // the hardWeight may be dependent on the instance (feature based)
    // Irace will need that a single hardWeight is used per instance,
    // do not compare different costs with different hardWeight on the same instance !
    double hardWeightFinal = -1; // < 0 means = hardWeight

    // Quick numbers of information
    MetaInfo meta;
    std::set<StudentId> students;
    MapVec<ExamId, MapVec<ExamId, int>> numberStudentsInCommon;

    // structures
    struct ExamsRelated {
        std::vector<ExamId> coincidences,
                         exclusions,
                         afters,
                         befores,
                         hasStudentsInCommon;

        void removeDuplicates();
    };

    MapVec<ExamId,ExamsRelated> examsRelated;

    std::vector<std::pair<PeriodId,PeriodId>> daysToPeriod;

    std::vector<ExamId>& coincidenceOfExam(ExamId i)         { return examsRelated[i].coincidences; }
    std::vector<ExamId>& exclusionOfExam(ExamId i)           { return examsRelated[i].exclusions; }
    std::vector<ExamId>& afterOfExam(ExamId i)               { return examsRelated[i].afters; }
    std::vector<ExamId>& beforeOfExam(ExamId i)              { return examsRelated[i].befores; }
    std::vector<ExamId>& hasStudentsInCommonOfExam(ExamId i) { return examsRelated[i].hasStudentsInCommon; }

    // const versions
    std::vector<ExamId> const& coincidenceOfExam(ExamId i)         const { return examsRelated[i].coincidences; }
    std::vector<ExamId> const& exclusionOfExam(ExamId i)           const { return examsRelated[i].exclusions; }
    std::vector<ExamId> const& afterOfExam(ExamId i)               const { return examsRelated[i].afters; }
    std::vector<ExamId> const& beforeOfExam(ExamId i)              const { return examsRelated[i].befores; }
    std::vector<ExamId> const& hasStudentsInCommonOfExam(ExamId i) const { return examsRelated[i].hasStudentsInCommon; }

    std::map<Minutes, Color> colorsOfMinute;

    std::vector<ExamId> frontLoadExams; // all exams that are in the largest
    MapVec<ExamId,bool> isInFrontLoad; // isInFrontLoad[i] iff i in frontLoadExams

    bool periodInTheEnd(PeriodId p) const;
    int sizeOfExam(ExamId e) const;

    void buildStructures();

    void finaliseSolution(Solution* solution) override;

    /////////////
    // Methods //
    /////////////
public:
    // compute meta and other informations
    void compute();

    int problemSize() override {
        return E();
    }

    // true if all periods are sorted by (date/hour) and function to sort it
    bool periodsSorted() const;
    // if we want to sort the periods, we should update the int pointers

    // true iff in a day, all periods are consecutive, periods
    // assert sorted
    bool allPeriodsConsecutives() const;

    // true iff all exams, periods, and rooms refer to an existing Exam, Period, Room
    bool correctIndexes() const;

    bool correctExam(ExamId) const;
    bool correctPeriod(PeriodId) const;
    bool correctRoom(RoomId) const;

    int numberOfExamsOfStudent(int student);

    // true iff all penalties are 0 or positive
    bool nonNegativePenalties();

    // access methods
    Two<Room&> roomsOf(Two<RoomId>);
    Two<Exam&> examsOf(Two<ExamId>);
    Two<Period&> periodsOf(Two<PeriodId>);

    Two<Room&> roomsOf(RoomId a, RoomId b);
    Two<Exam &> examsOf(ExamId a, ExamId b);
    Two<Period&> periodsOf(PeriodId a, PeriodId b);

    // const versions for c++ purist
    Two<Room const&> roomsOf(Two<RoomId>) const;
    Two<Exam const&> examsOf(Two<ExamId>) const;
    Two<Period const&> periodsOf(Two<PeriodId>) const;

    Two<Room const&> roomsOf(RoomId a, RoomId b) const;
    Two<Exam const&> examsOf(ExamId a, ExamId b) const;
    Two<Period const&> periodsOf(PeriodId a, PeriodId b) const;

    VectorExclusivePair<ExamId> pairsOfExams();

    int E() const { return exams.size(); }
    int P() const { return periods.size(); }
    int R() const { return rooms.size(); }

    int numberOfDurations() const { return colorsOfMinute.size(); }
    int numberOfDays() const { return this->daysToPeriod.size(); }

    int Dur() const { return numberOfDurations(); }
    int Days() const { return numberOfDays(); }

    void presentation(std::ostream&) const;
    void testDelta(ExamTTSolution&, std::ostream &log, int N = 1000, bool checkEachMove = true) const;

public:
    ExamTT() {

    }
    ExamTT(ExamTT const&) = delete;
    ExamTT(char* instance_path);
    virtual double evaluateSolution(Solution & solution);
};

typedef ExamTT Instance;

struct HardCostComponents {
    HCost
        simultaneousExams, // students
        overCapacity, // students (seats)
        excessiveDuration, // exams
        periodConstraintAfter, // pairs
        periodConstraintCoincidence, // pairs
        periodConstraintExclusion, // pairs
        roomConstraint; // exams

    static HardCostComponents zero() {
        return {0,0,0,0,0,0,0};
    }

    HCost periodConstraint() const {
        return periodConstraintAfter + periodConstraintCoincidence + periodConstraintExclusion;
    }

    HCost sum() const {
        return simultaneousExams + overCapacity + excessiveDuration + periodConstraint() + roomConstraint;
    }

    typedef std::tuple<HCost, HCost, HCost, HCost, HCost> SmallTuple;
    typedef std::tuple<HCost, HCost, HCost, HCost, HCost, HCost, HCost> Tuple;

    SmallTuple small_tuple() const {
        return std::make_tuple(simultaneousExams, overCapacity, excessiveDuration, periodConstraint(), roomConstraint);
    }

    Tuple make_tuple() const {
        return std::make_tuple(simultaneousExams, overCapacity, excessiveDuration, periodConstraintAfter, periodConstraintCoincidence, periodConstraintExclusion, roomConstraint);
    }

    bool exactlyEqual(HardCostComponents const& other) const {
        return make_tuple() == other.make_tuple();
    }

    struct Printer {
        HardCostComponents const& self;
        Instance const& instance;
    };

    Printer print(Instance const& instance) const {
        return Printer{*this, instance};
    }

    HardCostComponents operator -(HardCostComponents const& two) const {
        return {
            simultaneousExams           - two.simultaneousExams,
            overCapacity                - two.overCapacity,
            excessiveDuration           - two.excessiveDuration,
            periodConstraintAfter       - two.periodConstraintAfter,
            periodConstraintCoincidence - two.periodConstraintCoincidence,
            periodConstraintExclusion   - two.periodConstraintExclusion,
            roomConstraint              - two.roomConstraint,
        };
    }
};

struct SoftCostComponents {
    SCost twoExamsInARow,
          twoExamsInADay,
          periodSpread,
          periodsPenalty,
          roomsPenalty,
          frontload,
          mixedDuration;

    static SoftCostComponents zero() {
        return {0,0,0,0,0,0,0};
    }

    SCost sum() const {
        return twoExamsInARow + twoExamsInADay + periodSpread + periodsPenalty + roomsPenalty + frontload + mixedDuration;
    }

    typedef std::tuple<SCost, SCost, SCost, SCost> SmallTuple;
    typedef std::tuple<SCost, SCost, SCost, SCost, SCost, SCost, SCost> Tuple;

    SmallTuple small_tuple() const {
        return std::make_tuple(twoExamsInARow + twoExamsInADay + periodSpread, periodsPenalty + roomsPenalty, frontload, mixedDuration);
    }

    Tuple make_tuple() const {
        return std::make_tuple(twoExamsInARow, twoExamsInADay, periodSpread, periodsPenalty, roomsPenalty, frontload, mixedDuration);
    }

    bool exactlyEqual(SoftCostComponents const& other) const {
        return make_tuple() == other.make_tuple();
    }

    struct Printer {
        SoftCostComponents const& self;
        Instance const& instance;
    };

    Printer print(Instance const& instance) const {
        return Printer{*this, instance};
    }

    SoftCostComponents operator -(SoftCostComponents const& two) const {
        return {
            twoExamsInARow - two.twoExamsInARow,
            twoExamsInADay - two.twoExamsInADay,
            periodSpread   - two.periodSpread,
            periodsPenalty - two.periodsPenalty,
            roomsPenalty   - two.roomsPenalty,
            frontload      - two.frontload,
            mixedDuration  - two.mixedDuration,
        };
    }
};

std::ostream& operator <<(std::ostream& out, HardCostComponents::Printer printer);
std::ostream& operator <<(std::ostream& out, SoftCostComponents::Printer printer);

struct CostComponents {
    HardCostComponents hard;
    SoftCostComponents soft;

    static CostComponents zero() {
        return {HardCostComponents::zero(), SoftCostComponents::zero()};
    }

    double total(double hardWeight) const {
        return soft.sum() + hardWeight * hard.sum();
    }

    bool exactlyEqual(CostComponents const& other) const {
        return hard.exactlyEqual(other.hard) && soft.exactlyEqual(other.soft);
    }

    bool operator ==(CostComponents const& other) const {
        return exactlyEqual(other);
    }

    bool operator !=(CostComponents const& other) const {
        return ! exactlyEqual(other);
    }

    typedef std::tuple<HardCostComponents::Tuple, SoftCostComponents::Tuple> Tuple;
    typedef std::tuple<HCost, SCost> ComparableTuple;

    // (1, 0, 0, 2), (10, 0, ..., 5)
    Tuple make_tuple() const {
        return std::make_tuple(hard.make_tuple(), soft.make_tuple());
    }

    // (3, 160)
    ComparableTuple make_comparable_tuple() const {
        return std::make_tuple(hard.sum(), soft.sum());
    }

    bool operator <(CostComponents const& other) const {
        return make_comparable_tuple() < other.make_comparable_tuple();
    }

    struct Printer {
        CostComponents const& self;
        Instance const& instance;
    };

    Printer print(const Instance &instance) const {
        return Printer{*this, instance};
    }

    CostComponents operator -(CostComponents const& two) const {
        return {hard - two.hard, soft - two.soft};
    }
};

std::ostream& operator <<(std::ostream& out, CostComponents::Printer printer);

/*
 * Emili Solution
 */
class ExamTTSolution: public emili::Solution
{
public:
    ~ExamTTSolution();
    virtual const void* getRawData() const;
    virtual void setRawData(const void* data);

    typedef ExamTT const& InstanceRef;

    // assignement
    MapVec<ExamId, PeriodId> periods;
    MapVec<ExamId, RoomId> rooms;

    // meta
    CostComponents costs;

    // structures

    MapVec<PeriodId, MapVec<RoomId, std::list<ExamId>>> examsByPeriodRoom;
    std::list<ExamId> unAssignedExamList;

    MapVec<ExamId, std::list<ExamId>::iterator> examsIterator;
    MapVec<PeriodId, MapVec<RoomId, MapVec<Color,int>>> durationColorUsed;

    static constexpr bool USE_LAZY_STRUCTURES = false; // will not copy structures on clone, but will build structures when modified
    static constexpr bool USE_COLOR_STRUCTURE = true;
    static constexpr bool USE_DELTA = true;

    bool hasStructures = false;

    void buildStructures(InstanceRef);

    static int numberOfClones, numberOfTotalCompute;

    /**
     * return the difference in cost for a move/swap
     */
    CostComponents differenceCostMove(InstanceRef, ExamId ex, PeriodId nextP, RoomId nextR) const;
    CostComponents differenceCostSwap(InstanceRef, ExamId e1, ExamId e2) const;

    /**
     * Compute the difference in cost for <move> and update costs (+= -=).
     * Structures will be updated AFTER
     * update(...)
     * apply(...)
     */
    void updateMove(InstanceRef, ExamId ex, PeriodId nextP, RoomId nextR, CostComponents& costs) const;

    /**
     * @brief Compute the difference in cost for <remove/add>, does not update structures.
     * Structures will be updated AFTER
     * update(...)
     * apply(...)
     */
    void updateRemove(InstanceRef, ExamId);
    void updateAdd(InstanceRef, ExamId, PeriodId, RoomId);

    /**
     * @brief same as updateRemove but the exam is already removed from structures
     * apply(...)
     * update(...)
     */
    void updateRemoveOutOfStructures(InstanceRef, ExamId, PeriodId, RoomId);

    struct Printer { ExamTTSolution const& self; InstanceRef instance; };
    struct Writer { ExamTTSolution const& self; };

    Printer printer(InstanceRef i) const { return Printer{*this, i}; }
    Writer writer() const { return Writer{*this}; }

    // pretty print for user
    void printTo(InstanceRef, std::ostream &out) const;
    void printTimelineTo(InstanceRef, std::ostream& out) const;

    // itc format
    void writeTo(std::ostream &cout) const;

    // methods
    void initFromZeroPeriodAndZeroRoom(ExamTT const& instance);
    void initFromPeriodsAndZeroRoom(InstanceRef, std::vector<int> const& periods);
    void initFromAssign(InstanceRef, std::vector<std::pair<int,int>> assign);
    void initRandom(InstanceRef, Random&);

    void initRandom(InstanceRef);
    void initUnassigned(InstanceRef);

    void move(InstanceRef, ExamId, PeriodId, RoomId);
    void swap(InstanceRef, ExamId, ExamId);
    void movePeriod(InstanceRef, ExamId, PeriodId);
    void moveRoom(InstanceRef, ExamId, RoomId);

    void removeExam(InstanceRef, ExamId);
    void addExam(InstanceRef, ExamId, PeriodId, RoomId);
    void addExam(InstanceRef, ExamId, Assignement);;

    void refreshSolutionValue(InstanceRef);

    void applyMove(InstanceRef, ExamId, PeriodId, RoomId);
    void applyRemoveExam(InstanceRef, ExamId, PeriodId prevP, RoomId prevR);
    void applyAddExam(InstanceRef, ExamId, PeriodId, RoomId);

    // compute
    /**
     * @brief full compute of cost, if unassigned list is not empty, will change to computeCostPartial
     * @param costs
     */
    void computeCost(InstanceRef, CostComponents& costs) const;
    void computeCost(InstanceRef);
    CostComponents computeAndGetCost(InstanceRef) const;

    void computeCostPartial(InstanceRef, CostComponents& costs) const;

    // read
    AssignementRef assignement(ExamId exam);

    Two<AssignementRef> assignementOf(Two<ExamId> exam);
    Two<PeriodId&> periodsOf(Two<ExamId> exams);
    Two<RoomId&> roomsOf(Two<ExamId> exams);

    Two<AssignementRef> assignementOf(ExamId exam1, ExamId exam2);
    Two<PeriodId&> periodsOf(ExamId exam1, ExamId exam2);
    Two<RoomId&> roomsOf(ExamId exam1, ExamId exam2);

    // const version, for c++ purist
    AssignementConstRef assignement(ExamId exam) const;

    Two<AssignementConstRef> assignementOf(Two<ExamId> exam) const;
    Two<PeriodId const&> periodsOf(Two<ExamId> exams) const;
    Two<RoomId const&> roomsOf(Two<ExamId> exams) const;

    Two<AssignementConstRef> assignementOf(ExamId exam1, ExamId exam2) const;
    Two<PeriodId const&> periodsOf(ExamId exam1, ExamId exam2) const;
    Two<RoomId const&> roomsOf(ExamId exam1, ExamId exam2) const;

    int sizeOfPartialSolution() const;
    inline bool isAssigned(ExamId e) const { return ~periods[e]; /* periods[e] != -1; */ }
    inline bool isFullyAssigned(ExamId e) const { return ~periods[e] & ~rooms[e]; /* periods[e] != -1 && rooms[e] != -1; */ }

public:
    ExamTTSolution(double solution_value = 0):emili::Solution(solution_value) { }

    /* O(n) */
    Solution* clone() override;

    /* O(1) */
    void swap(Solution*) override;

    std::string getSolutionRepresentation() override;
};

struct InstanceParser {
    std::ifstream file;
    InstanceParser(std::string filename) : file(filename.c_str()) {}
    void open(std::string filename) { file.open(filename.c_str()); }

    // Throw exception if something bad happened
    void parse(ExamTT& inst);
    void parse(const ExamTT &instance, ExamTTSolution& sol);

    // Read next line as [{name}:{i}] and return i
    int readHeaderInt(std::string name);
    // Read next line as [{name}] and return void
    void readHeaderVoid(std::string name);
    // Read next line as [{name}] and return name
    std::string readHeaderName();

    // read line and throw exception if no line
    void readline(std::string& line);
    std::string readline();
};

std::ostream& operator <<(std::ostream&, ExamTTSolution::Printer);
std::ostream& operator <<(std::ostream&, ExamTTSolution::Writer);

struct RandomInitialSolution : emili::InitialSolution {
    Random random;

    RandomInitialSolution(ExamTT& inst) : emili::InitialSolution(inst) {}

    Solution* generateSolution() override;
    Solution* generateEmptySolution() override;
};

struct ZeroInitialSolution : emili::InitialSolution {
    ZeroInitialSolution(ExamTT& inst) : ZeroInitialSolution(inst, {}) {}

    ZeroInitialSolution(ExamTT& inst, std::vector<Assignement> firstAssign_)
        : emili::InitialSolution(inst), firstAssign(firstAssign_) {}

    std::vector<Assignement> firstAssign;

    Solution* generateSolution() override;
    Solution* generateEmptySolution() override;
};

struct MoveNeighborhood : emili::Neighborhood {
protected:
    ExamTT const& instance;

    ExamId exam, re;
    PeriodId bperiod, period, rp;
    RoomId broom, room, rr;

    void iterate(Solution *base, std::function<void ()> yield) override;
    Solution* computeStep(Solution *step) override;
    void reverseLastMove(Solution *step) override;
public:
    MoveNeighborhood(ExamTT const& instance);

    NeighborhoodIterator begin(emili::Solution* base) override;
    Solution* step(Solution *currentSolution) override;
    void reset() override;
    int size() override;

    Solution* random(Solution *currentSolution) override;
    void randomStep(Solution* currentSolution) override;
    void reverseLastRandomStep(Solution *currentSolution) override;
};

struct SwapNeighborhood : emili::Neighborhood {
protected:
    ExamTT const& instance;

    ExamId e1;
    ExamId e2;

    ExamId re1, re2;

    void iterate(Solution *base, std::function<void ()> yield) override;
    Solution* computeStep(Solution *rawStep) override;
    void reverseLastMove(Solution *rawStep) override;

public:
    SwapNeighborhood(ExamTT const& instance);

    Solution* step(Solution *currentSolution) override;
    void reset() override;
    int size() override;

    Solution* random(Solution *currentSolution) override;
    void randomStep(Solution* currentSolution) override;
    void reverseLastRandomStep(Solution *currentSolution) override;
};

/**
 * Implements iteration as :
    for exam in range(E)
        t0 = periods[exam]
        for t1 in range(P)
            if t1 != t0
                chain = createChain(exam, t0, t1)
                swapColors(chain, t0, t1)
                yield()
                swapColors(chain, t0, t1)

    This creates multiple duplicates, see FastIter

 * Implements random as :
   e = randrange(E)
   t0 = periods[e]
   t1 = randrange_different(P, t0)
   chain = createChain(e, t0, t1)
   swapInChain(chain, t0, t1)
 */
struct KempeChainNeighborhood : emili::Neighborhood {
protected:
    ExamTT const& instance;

    ExamId exam;
    PeriodId t0, t1;
    std::set<ExamId> chain;

    /**
     * assert period[x] in (t0,t1)
     * insert x to chain
     * repeat with neighbors of x with period[j] in (t0,t1)
     */
    void createChain(ExamTTSolution* sol, ExamId x);

    Solution* computeStep(Solution *rawStep) override;
    void reverseLastMove(Solution *rawStep) override;

    void swapPeriodsInChain(ExamTTSolution *sol);

public:
    KempeChainNeighborhood(ExamTT const& instance);

    void iterate(Solution*, std::function<void()>) override;
    Solution* step(Solution *currentSolution) override;
    void reset() override;
    int size() override;

    Solution* random(Solution *currentSolution) override;
    void randomStep(Solution* currentSolution) override;
    void reverseLastRandomStep(Solution *currentSolution) override;
};

/**
 * Implement iteration as :
    for i in range(P)
        for j in range(i+1,P)
            chains = components(C[i] | C[j]) # C[i] = all exams at period i
            for chain in chains:
                swapColors(chain, i, j)
                yield()
                swapColors(chain, i, j)
 */
struct KempeChainNeighborhoodFastIter : KempeChainNeighborhood {
private:
    std::vector<std::set<ExamId>> examsByPeriod;
    std::set<ExamId> A;
    std::list<ExamId> chainList;

    void swapPeriodsInChainList(ExamTTSolution *sol);
public:
    KempeChainNeighborhoodFastIter(ExamTT const& instance) : KempeChainNeighborhood(instance) {}

    void reset() override;
    void iterate(Solution* base, std::function<void()> yield) override;
    Solution* computeStep(Solution *rawStep) override;
    void reverseLastMove(Solution *rawStep) override;
    void createChainDestruct(ExamId e);
};

struct RandomNeighborhood : emili::Neighborhood {
protected:
    emili::Solution* computeStep(Solution *step) override {
        // undefined for random neigh
        return nullptr;
    }

    void reverseLastMove(Solution *step) override {
        // undefined for random neigh
    }

public:

    Solution* step(Solution *currentSolution) override {
        // undefined for random neigh
        return nullptr;
    }

    void reset() override {
        // the reset is done at the first random() call
    }
};

struct MixedMoveSwapNeighborhood : RandomNeighborhood {
public:
    MoveNeighborhood move;
    SwapNeighborhood swap;

    double swapRate;

    bool chosenSwap;

    Solution* computeStep(Solution *rawStep) override {
        throw std::invalid_argument("undefined for random neighborhood");
    }

    void reverseLastMove(Solution *rawStep) override {
        throw std::invalid_argument("undefined for random neighborhood");
    }
public:
    MixedMoveSwapNeighborhood(ExamTT const& instance, double swapRate);

    Solution* step(Solution *currentSolution) override {
        throw std::invalid_argument("undefined for random neighborhood");
    }

    void reset() override {
        // undefined for random neighborhood
    }

    int size() override {
        return move.size() + swap.size();
    }

    Solution* random(Solution *currentSolution) override;
    void randomStep(Solution* currentSolution) override;
    void reverseLastRandomStep(Solution *currentSolution) override;
};

struct MixedRandomNeighborhood : RandomNeighborhood {
protected:
    std::vector<emili::Neighborhood*> neighborhoods;
    std::vector<int> cumul;
    int i = 0;
    int _size = 0;

public:

    MixedRandomNeighborhood(std::vector<emili::Neighborhood*> ns, std::vector<int> weights);

    ~MixedRandomNeighborhood() {
        for(auto p : neighborhoods)
            delete p;
    }

    void clearVector() {
        neighborhoods.clear();
    }

    int size() override {
        return _size;
    }

    Solution* random(Solution* currentSolution) override;
    void randomStep(Solution* currentSolution) override;
    void reverseLastRandomStep(Solution *currentSolution) override;
};

struct MixedRandomNeighborhoodProba : RandomNeighborhood {
protected:
    std::vector<emili::Neighborhood*> neighborhoods;
    std::vector<float> cumul;
    int i = 0;
    int _size = 0;

public:

    // proba.size = neigh.size - 1
    MixedRandomNeighborhoodProba(std::vector<emili::Neighborhood*> ns, std::vector<float> proba);

    ~MixedRandomNeighborhoodProba() {
        for(auto p : neighborhoods)
            delete p;
    }

    void clearVector() {
        neighborhoods.clear();
    }

    int size() override {
        return _size;
    }

    Solution* random(Solution* currentSolution) override;
    void randomStep(Solution* currentSolution) override;
    void reverseLastRandomStep(Solution *currentSolution) override;
};

struct FixedRandomDestructor : Destructor {
private:
    Instance const& instance;
    std::vector<ExamId> inserted;
public:
    FixedRandomDestructor(Instance const&, const int G_);
public:
    Solution* destruct(Solution *solution) override;
};

// TODO: other ideas : N grouped days / grouped periods
struct NRandomDaysDestructor : Destructor {
private:
    Instance const& instance;
    const int N;
    std::vector<int> days;
public:
    NRandomDaysDestructor(Instance const& instance_, const int N_);
public:
    Solution* destruct(Solution *solution) override;
};

struct NGroupedDaysDestructor : Destructor {
private:
    Instance const& instance;
    const int N;
public:
    NGroupedDaysDestructor(Instance const& instance_, const int N_);
public:
    Solution* destruct(Solution *solution) override;
};

struct NGroupedPeriodsDestructor : Destructor {
private:
    Instance const& instance;
    const int N;
public:
    NGroupedPeriodsDestructor(Instance const& instance_, const int N_);
public:
    Solution* destruct(Solution *solution) override;
};

struct NBiggestWeightedDestructor : Destructor {
private:
    Instance const& instance;
    const int N;
public:
    NBiggestWeightedDestructor(Instance const& instance_, const int N_);
public:
    Solution* destruct(Solution *solution) override;
};

class InsertHeuristic {
protected:
    ExamTT const& instance;
public:
    InsertHeuristic(ExamTT const& instance_) : instance(instance_) {}
public:
    virtual Assignement searchPosition(ExamTTSolution*, ExamId)=0;
};

struct BaseConstructor : Constructor {
    Instance const& instance;

    BaseConstructor(Instance const& inst) : instance(inst) {
        const_cast<Behaviour&>(behaviour) = Behaviour::VOID;
    }

    Solution* constructFull() override {
        auto sol = new ExamTTSolution;
        sol->initUnassigned(instance);
        auto c = construct(sol);
        if(c != sol)
            delete sol;
        return c;
    }
};

struct RandomOrderInserter : BaseConstructor {
private:
    InsertHeuristic* insertHeuristic;
    std::vector<ExamId> inserted;
public:
    RandomOrderInserter(Instance const& instance_, InsertHeuristic* insertHeuristic_);
public:
    Solution* construct(Solution *partial) override;
};

struct DegreeInserter : BaseConstructor {
private:
    InsertHeuristic* insertHeuristic;
    std::vector<Two<int>> data;
public:
    DegreeInserter(Instance const& instance_, InsertHeuristic* insertHeuristic_) : BaseConstructor(instance_), insertHeuristic(insertHeuristic_) {}
public:
    Solution* construct(Solution *partial) override;
};

struct DSaturInserter : BaseConstructor {
private:
    InsertHeuristic* insertHeuristic;
public:
    DSaturInserter(Instance const& instance_, InsertHeuristic* insertHeuristic_) : BaseConstructor(instance_), insertHeuristic(insertHeuristic_) {}
public:
    Solution* construct(Solution *partial) override;
};

struct BestInsertHeuristic : InsertHeuristic {
public:
    BestInsertHeuristic(ExamTT const& instance) : InsertHeuristic(instance) {}
public:
    Assignement searchPosition(ExamTTSolution *, ExamId) override;
};

struct ConstructorInitialSolution : emili::InitialSolution {
    Constructor* constructor;
    ConstructorInitialSolution(ExamTT& inst, Constructor* c)
        : emili::InitialSolution(inst), constructor(c) {}

    Solution* generateSolution() override;
    Solution* generateEmptySolution() override;
};

/**
 * @brief Destruct randomly, constructs by greedy method
 */
struct IteratedGreedyNeighborhood : public emili::Neighborhood {
protected:
    ExamTT const& instance;

    Solution* computeStep(Solution *step) override;
    void reverseLastMove(Solution *step) override;

    const int G;
    std::vector<ExamId> inserted;
public:
    IteratedGreedyNeighborhood(ExamTT const& instance_, const int G_) : instance(instance_), G(G_) {
        reset();
    }

    ~IteratedGreedyNeighborhood() {}

    void reset() override;


    int size() override {
        return 0;
    }

    Solution* random(Solution* currentSolution) override;
    void randomStep(Solution* currentSolution) override;
    void reverseLastRandomStep(Solution *currentSolution) override;
};

struct BruteForce : emili::LocalSearch {
    Instance& instance;
    std::vector<std::pair<PeriodId,RoomId>> firstAssign;
    int sizeOfDebug = 2;

    BruteForce(Instance& i);
    BruteForce(Instance& i, std::vector<std::pair<PeriodId,RoomId>> firstAssign);
    void setSizeDebug(int n) { sizeOfDebug = n; }
    Solution* search(Solution* initial) override;
};

namespace stats {
void kempe_print_iteration(ExamTT& instance, std::vector<int> initPeriods, bool useFastIter=false, bool useIterate=false);
void kempe_compare_size_fast_iter(ExamTT& instance);
}

namespace test {
void delta(const ExamTT &inst, int N, bool checkEachMove);
void deltaRemoveAdd(const ExamTT &inst, int N, bool checkEachMove, int G, bool checkEachMovePartial);
void interactive(const ExamTT &inst);

void kempe_iteration_vs_random(ExamTT& instance, std::vector<int> initPeriods, int N=1000, bool useFastIter=false, bool useIterate=false);
void kempe_iteration_vs_random(ExamTT& instance, ExamTTSolution& sol, int N=1000, bool useFastIter=false, bool useIterate=false);

void iterateVsComputeStep(ExamTTSolution* solution, Neighborhood* neigh);
void constructUnassignedVsRemoveAll(ExamTT& instance);
}

}
}
#endif // ExamTT_H