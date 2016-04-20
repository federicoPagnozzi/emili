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

// SA
#include "../SA/sa.h"

namespace emili
{
namespace ExamTT
{


// Output iterator that counts the number of times ++ is called
// can be casted to a int to get the value or get attribute .count
struct CountIterator {
    int count = 0; // contain the numer of times ++ has been applied
    int dummyWrite; // contain last value written

    int&           operator*()          { return dummyWrite; }
    int const&     operator*() const    { return dummyWrite; }
    CountIterator& operator++()         { ++count; return *this; }
                   operator int() const { return count; }
};

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
        return {0,0,N};
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

    // will be featured based, now it's just a constant
    double hardWeight = 0;

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

    /////////////
    // Methods //
    /////////////
public:
    // compute meta and other informations
    void compute();

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

    MapVec<ExamId, std::list<ExamId>::iterator> examsByPeriodRoomIterators;
    MapVec<PeriodId, MapVec<RoomId, MapVec<Color,int>>> durationColorUsed;

    static constexpr bool USE_LAZY_STRUCTURES = true; // will not copy structures on clone, but will build structures when modified
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
    void initFromPeriodsAndZeroRoom(InstanceRef, std::vector<int> const& periods);
    void initFromAssign(InstanceRef, std::vector<std::pair<int,int>> assign);
    void initRandom(InstanceRef, Random&);

    void move(InstanceRef, ExamId, PeriodId, RoomId);
    void swap(InstanceRef, ExamId, ExamId);
    void movePeriod(InstanceRef, ExamId, PeriodId);
    void moveRoom(InstanceRef, ExamId, RoomId);

    void removeExam(InstanceRef, ExamId);
    void addExam(InstanceRef, ExamId, PeriodId, RoomId);

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

    int sizeOfPartialSolution() const;;

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

struct MoveNeighborhood : emili::Neighborhood {
protected:
    ExamTT const& instance;

    ExamId exam, re;
    PeriodId bperiod, period, rp;
    RoomId broom, room, rr;

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

struct KempeChainNeighborhood : emili::Neighborhood {
protected:
    ExamTT const& instance;

    ExamId exam;
    int t0, t1;
    std::set<ExamId> chain;

    void createChain(ExamTTSolution* sol, ExamId x);

    Solution* computeStep(Solution *rawStep) override;
    void reverseLastMove(Solution *rawStep) override;

public:
    KempeChainNeighborhood(ExamTT const& instance);

    Solution* step(Solution *currentSolution) override;
    void reset() override;
    int size() override;

    Solution* random(Solution *currentSolution) override;
    void randomStep(Solution* currentSolution) override;
    void reverseLastRandomStep(Solution *currentSolution) override;
};

struct MixedMoveSwapNeighborhood : emili::Neighborhood {
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
        return 0; // // undefined for random neighborhood
    }

    Solution* random(Solution *currentSolution) override;
    void randomStep(Solution* currentSolution) override;
    void reverseLastRandomStep(Solution *currentSolution) override;
};

struct MixedRandomNeighborhood : emili::Neighborhood {
protected:
    std::vector<emili::Neighborhood*> neighborhoods;
    std::vector<int> cumul;
    int i = 0;

public:

    MixedRandomNeighborhood(std::vector<emili::Neighborhood*> n, std::vector<int> weights);

    ~MixedRandomNeighborhood() {
        for(auto p : neighborhoods)
            delete p;
    }

    void clearVector() {
        neighborhoods.clear();
    }

    void reset() override {
        // undefined for random neighborhood
    }

    int size() override {
        return 0; // undefined for random neighborhood
    }

    Solution* random(Solution* currentSolution) override;
    void randomStep(Solution* currentSolution) override;
    void reverseLastRandomStep(Solution *currentSolution) override;
};

struct MixedRandomNeighborhoodProba : emili::Neighborhood {
protected:
    std::vector<emili::Neighborhood*> neighborhoods;
    std::vector<float> cumul;
    int i = 0;

public:

    // proba.size = neigh.size - 1
    MixedRandomNeighborhoodProba(std::vector<emili::Neighborhood*> n, std::vector<float> proba);

    ~MixedRandomNeighborhoodProba() {
        for(auto p : neighborhoods)
            delete p;
    }

    void clearVector() {
        neighborhoods.clear();
    }

    void reset() override {
        // undefined for random neighborhood
    }

    int size() override {
        return 0; // undefined for random neighborhood
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

struct BSUSA : emili::LocalSearch {
    BSUSA(
        Instance& i,
        emili::InitialSolution* init,
        SAInitTemp* initialTemperature,
        SATempLength* tempLegth,
        SACooling* cooling
    );
    Solution* search(Solution* initial) override;
    void searchInPlace(Solution* initial) override;
};

namespace test {
void delta(const ExamTT &inst, int N, bool checkEachMove);
void deltaRemoveAdd(const ExamTT &inst, int N, bool checkEachMove, int G);
void interactive(const ExamTT &inst);
void kempe(ExamTT& instance, std::vector<int> initPeriods);

void kempe_iteration_vs_random(ExamTT& instance, std::vector<int> initPeriods, int N=1000);
void kempe_iteration_vs_random(ExamTT& instance, ExamTTSolution& sol, int N=1000);
}

}
}
#endif // ExamTT_H
