#ifndef ExamTT_H
#include "../emilibase.h"
#define ExamTT_H

#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <algorithm>

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

    void seed(long s) {
       emili::initializeRandom(s);
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
                if(a == A || ++b == B)
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

template <typename T>
std::set<T> intersection(std::set<T> const& one, std::set<T> const& two) {
    std::set<T> res;
    set_intersection(one.begin(), one.end(), two.begin(), two.end(), inserter(res, res.end()));
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
typedef int Cost; // A soft or hard cost
typedef int HCost; // Hard
typedef int SCost; // Soft

typedef int NumberOfStudents;
typedef int NumberOfPeriods;
typedef int NumberOfRooms;
typedef int NumberOfExams;


class Exam {
public:
    Minutes duration;
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
    Date date;
    Time time;
    Minutes duration;
    SCost penalty;

    static bool compareDateTime(Period const& p1, Period const& p2) {
        return std::make_tuple(p1.date, p2.time) < std::make_tuple(p2.date, p2.time);
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
        std::set<ExamId> coincidences,
                         exclusions,
                         afters,
                         befores,
                         hasStudentsInCommon;
    };

    MapVec<ExamId,ExamsRelated> examsRelated;

    std::set<ExamId>& coincidenceOfExam(ExamId i)         { return examsRelated[i].coincidences; }
    std::set<ExamId>& exclusionOfExam(ExamId i)           { return examsRelated[i].exclusions; }
    std::set<ExamId>& afterOfExam(ExamId i)               { return examsRelated[i].afters; }
    std::set<ExamId>& beforeOfExam(ExamId i)              { return examsRelated[i].befores; }
    std::set<ExamId>& hasStudentsInCommonOfExam(ExamId i) { return examsRelated[i].hasStudentsInCommon; }

    // const versions
    std::set<ExamId> const& coincidenceOfExam(ExamId i)         const { return examsRelated[i].coincidences; }
    std::set<ExamId> const& exclusionOfExam(ExamId i)           const { return examsRelated[i].exclusions; }
    std::set<ExamId> const& afterOfExam(ExamId i)               const { return examsRelated[i].afters; }
    std::set<ExamId> const& beforeOfExam(ExamId i)              const { return examsRelated[i].befores; }
    std::set<ExamId> const& hasStudentsInCommonOfExam(ExamId i) const { return examsRelated[i].hasStudentsInCommon; }

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

public:
    ExamTT() {}
    ExamTT(char* instance_path);
    virtual double evaluateSolution(Solution & solution);
};

typedef ExamTT Instance;

struct HardCostComponents {
    HCost roomConstraint, // exams
        excessiveDuration, // exams
        periodConstraintAfter, // pairs
        periodConstraintCoincidence, // pairs
        periodConstraintExclusion, // pairs
        simultaneousExams, // students
        overCapacity; // students (seats)

    static HardCostComponents zero() {
        return {0,0,0,0,0,0,0};
    }

    HCost periodConstraint() const {
        return periodConstraintAfter + periodConstraintCoincidence + periodConstraintExclusion;
    }

    HCost sum() const {
        return roomConstraint + excessiveDuration + periodConstraint() + simultaneousExams + overCapacity;
    }

    std::tuple<HCost, HCost, HCost, HCost, HCost> small_tuple() const {
        return std::make_tuple(roomConstraint, excessiveDuration, periodConstraint(), simultaneousExams, overCapacity);
    }

    std::tuple<HCost, HCost, HCost, HCost, HCost, HCost, HCost> make_tuple() const {
        return std::make_tuple(roomConstraint, excessiveDuration, periodConstraintAfter, periodConstraintCoincidence, periodConstraintExclusion, simultaneousExams, overCapacity);
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
            roomConstraint              - two.roomConstraint,
            excessiveDuration           - two.excessiveDuration,
            periodConstraintAfter       - two.periodConstraintAfter,
            periodConstraintCoincidence - two.periodConstraintCoincidence,
            periodConstraintExclusion   - two.periodConstraintExclusion,
            simultaneousExams           - two.simultaneousExams,
            overCapacity                - two.overCapacity
        };
    }
};

struct SoftCostComponents {
    SCost periodsPenalty,
          roomsPenalty,
          twoExamsInARow,
          twoExamsInADay,
          periodSpread,
          frontload,
          mixedDuration;

    static SoftCostComponents zero() {
        return {0,0,0,0,0,0,0};
    }

    SCost sum() const {
        return periodsPenalty + roomsPenalty + twoExamsInARow + twoExamsInADay + periodSpread + frontload + mixedDuration;
    }

    std::tuple<SCost, SCost, SCost, SCost> small_tuple() const {
        return std::make_tuple(periodsPenalty + roomsPenalty, twoExamsInARow + twoExamsInADay + periodSpread, frontload, mixedDuration);
    }

    std::tuple<SCost, SCost, SCost, SCost, SCost, SCost, SCost> make_tuple() const {
        return std::make_tuple(periodsPenalty, roomsPenalty, twoExamsInARow, twoExamsInADay, periodSpread, frontload, mixedDuration);
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
            periodsPenalty - two.periodsPenalty,
            roomsPenalty   - two.roomsPenalty,
            twoExamsInARow - two.twoExamsInARow,
            twoExamsInADay - two.twoExamsInADay,
            periodSpread   - two.periodSpread,
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
    virtual const void* getRawData() const;
    virtual void setRawData(const void* data);

    typedef ExamTT const& InstanceRef;

    // assignement
    MapVec<ExamId, PeriodId> periods;
    MapVec<ExamId, RoomId> rooms;

    // meta
    CostComponents costs;

    // structures
    MapVec<PeriodId, std::set<ExamId>> examsByPeriods;
    MapVec<ExamId, std::set<ExamId>::iterator> examsByPeriodsIterators;
    MapVec<RoomId, std::set<ExamId>> examsByRooms;
    MapVec<ExamId, std::set<ExamId>::iterator> examsByRoomsIterators;

    void buildStructures(InstanceRef);

    /**
     * return the difference in cost for a move/swap
     */

    CostComponents differenceCostMove(InstanceRef, ExamId ex, PeriodId nextP, RoomId nextR) const;
    CostComponents differenceCostSwap(InstanceRef, ExamId e1, ExamId e2) const;

    /**
     * Compute the difference in cost for <move> and update costs (+= -=).
     */
    void updateMove(InstanceRef, ExamId ex, PeriodId nextP, RoomId nextR, CostComponents& costs) const;

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
    void initRandom(InstanceRef, Random&);
    void move(InstanceRef, ExamId e, PeriodId p, RoomId r);
    void swap(InstanceRef, ExamId i, ExamId j);

    void movePeriod(InstanceRef, ExamId e, PeriodId p);
    void moveRoom(InstanceRef, ExamId e, RoomId r);

    // compute
    void computeCost(InstanceRef, CostComponents& costs) const;
    void computeCost(InstanceRef);
    CostComponents computeAndGetCost(InstanceRef) const;

    // read
    std::pair<PeriodId &, RoomId &> assignement(ExamId exam);

    Two<std::pair<PeriodId &, RoomId &>> assignementOf(Two<ExamId> exam);
    Two<PeriodId&> periodsOf(Two<ExamId> exams);
    Two<RoomId&> roomsOf(Two<ExamId> exams);

    Two<std::pair<PeriodId&, RoomId&>> assignementOf(ExamId exam1, ExamId exam2);
    Two<PeriodId&> periodsOf(ExamId exam1, ExamId exam2);
    Two<RoomId&> roomsOf(ExamId exam1, ExamId exam2);

    // const version, for c++ purist
    std::pair<PeriodId const&, RoomId const&> assignement(ExamId exam) const;

    Two<std::pair<PeriodId const&, RoomId const&>> assignementOf(Two<ExamId> exam) const;
    Two<PeriodId const&> periodsOf(Two<ExamId> exams) const;
    Two<RoomId const&> roomsOf(Two<ExamId> exams) const;

    Two<std::pair<PeriodId const&, RoomId const&>> assignementOf(ExamId exam1, ExamId exam2) const;
    Two<PeriodId const&> periodsOf(ExamId exam1, ExamId exam2) const;
    Two<RoomId const&> roomsOf(ExamId exam1, ExamId exam2) const;

public:
    ExamTTSolution(double solution_value = 0):emili::Solution(solution_value) { }
    virtual Solution* clone();

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

    ExamId exam;
    PeriodId bperiod, period;
    RoomId broom, room;

    Solution* computeStep(Solution *step) override;
    void reverseLastMove(Solution *step) override;
public:
    MoveNeighborhood(ExamTT const& instance);

    NeighborhoodIterator begin(emili::Solution* base) override;
    Solution* step(Solution *currentSolution) override;
    void reset() override;
    int size() override;

    Solution* random(Solution *currentSolution) override;
};

struct SwapNeighborhood : emili::Neighborhood {
protected:
    ExamTT const& instance;

    ExamId e1;
    ExamId e2;

    Solution* computeStep(Solution *rawStep) override;
    void reverseLastMove(Solution *rawStep) override;

public:
    SwapNeighborhood(ExamTT const& instance);

    Solution* step(Solution *currentSolution) override;
    void reset() override;
    int size() override;

    Solution* random(Solution *currentSolution) override;
};

struct RandomInitialSolution : emili::InitialSolution {
    Random random;

    RandomInitialSolution(ExamTT& inst) : emili::InitialSolution(inst) {}

    Solution* generateSolution() override;
    Solution* generateEmptySolution() override;
};

void test();

}
}
#endif // ExamTT_H
