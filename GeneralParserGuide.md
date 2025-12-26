# EMILI General Parser Guide

A developer guide for understanding the parser system and adding new components to the EMILI framework.

## Overview

EMILI uses a **component-based architecture** where algorithms are assembled from modular pieces. The parser transforms command-line tokens into fully instantiated algorithm objects.

**Command-line syntax:**
```bash
./emili INSTANCE_FILE PROBLEM ALGORITHM_TOKENS [-it|-ro time] [rnds seed] [ps]
```

**Example:**
```bash
./emili Ta086.txt PFSP_MS first neh locmin insert -it 10 rnds 42
```

This creates a First Improvement local search with NEH initialization, local minima termination, and insert neighborhood.

---

## Architecture

```
main.cpp
    │
    ├── GeneralParserE ─────────────────────────────────────┐
    │       │                                               │
    │       ├── TokenManager                                │
    │       │   (handles command-line tokens)               │
    │       │                                               │
    │       └── Builder[] ──────────────────────┐           │
    │           │                               │           │
    │           ├── EmBaseBuilder               │           │
    │           │   (problem-independent:       │           │
    │           │    ILS, Tabu, VND, etc.)      │           │
    │           │                               │           │
    │           └── PfspBuilder                 │           │
    │               (PFSP-specific:             │           │
    │                NEH, Insert, Exchange)     │           │
    │                                           │           │
    └── LocalSearch ◄───────────────────────────┴───────────┘
        (assembled algorithm)
```

**Key classes:**
- `GeneralParserE` - Orchestrates parsing, manages builders
- `TokenManager` - Handles token navigation and parsing
- `Builder` - Abstract factory for creating components
- `Component` - Wrapper holding type info and object pointer

---

## Core Concepts

### TokenManager

Manages command-line arguments with navigation methods:

| Method | Description |
|--------|-------------|
| `nextToken()` | Consume and return next token |
| `checkToken(str)` | Return true and consume if current matches `str` |
| `peek()` / `operator*()` | Return current token without consuming |
| `getInteger()` / `getDecimal()` | Parse numeric values |
| `hasMoreTokens()` | Check if tokens remain |

### Builder

Abstract factory that creates components. Each problem has its own Builder:

| Method | Returns | Purpose |
|--------|---------|---------|
| `isCompatibleWith(problem)` | bool | Check if builder handles this problem |
| `canOpenInstance(problem)` | bool | Check if builder can load instance |
| `openInstance()` | Problem* | Load problem instance from file |
| `buildAlgo()` | LocalSearch* | Build algorithm from tokens |
| `buildNeighborhood()` | Neighborhood* | Build neighborhood |
| `buildInitialSolution()` | InitialSolution* | Build initial solution generator |
| `buildTermination()` | Termination* | Build termination criterion |
| `buildPerturbation()` | Perturbation* | Build perturbation operator |
| `buildAcceptance()` | Acceptance* | Build acceptance criterion |

### Component Types

| Type | Hex Value | Description |
|------|-----------|-------------|
| `COMPONENT_ALGORITHM` | 0xA0 | Local search algorithms |
| `COMPONENT_INITIAL_SOLUTION_GENERATOR` | 0xB1 | Solution construction |
| `COMPONENT_TERMINATION_CRITERION` | 0xB2 | Stopping conditions |
| `COMPONENT_NEIGHBORHOOD` | 0xB3 | Solution neighborhoods |
| `COMPONENT_PERTURBATION` | 0xC1 | Perturbation operators |
| `COMPONENT_ACCEPTANCE` | 0xC2 | Acceptance criteria |
| `COMPONENT_TABU_TENURE` | 0xB4 | Tabu search memory |

---

## Parsing Flow

```
Command: ./emili Ta086.txt PFSP_MS first neh locmin insert

Tokens:  [Ta086.txt] [PFSP_MS] [first] [neh] [locmin] [insert]
              │          │        │       │      │        │
              │          │        ▼       │      │        │
              │          │   buildAlgo()  │      │        │
              │          │   "first" ──►  │      │        │
              │          │        │       │      │        │
              │          │        ├───────┴──────┴────────┤
              │          │        │                       │
              │          │        ▼                       ▼
              │          │   retrieveComponent()     (recursive)
              │          │        │
              │          │        ├─► buildInitialSolution() → "neh" → NEH object
              │          │        ├─► buildTermination() → "locmin" → LocalMinima object
              │          │        └─► buildNeighborhood() → "insert" → Insert object
              │          │
              ▼          ▼
         Load file   Select compatible builders
                     (EmBaseBuilder + PfspBuilder)

Result: FirstImprovementSearch(NEH, LocalMinima, InsertNeighborhood)
```

**Process:**
1. `GeneralParserE::parseParams()` reads problem token
2. Finds all compatible builders via `isCompatibleWith()`
3. Loads instance via first builder that `canOpenInstance()`
4. Calls `buildComponent(COMPONENT_ALGORITHM)` to start recursive build
5. Each builder's `buildAlgo()` checks tokens and retrieves sub-components
6. Returns assembled `LocalSearch` object

---

## Using the Template Directory

The `template/` directory provides starter files for new problems:

```
template/
├── problem_template.h      # Problem, Solution, Neighborhood, Perturbation classes
├── problem_template.cpp    # Implementation stubs
├── problem_builder.h       # Builder class template
└── problem_builder.cpp     # Builder implementation stubs
```

**Classes in `problem_template.h`:**
- `ProblemX` - Problem definition (inherits `emili::Problem`)
- `SolutionProblemX` - Solution representation (inherits `emili::Solution`)
- `InitialSolutionProblemX` - Solution generator (inherits `emili::InitialSolution`)
- `NeighborhoodProblemX` - Neighborhood structure (inherits `emili::Neighborhood`)
- `PerturbationProblemX` - Perturbation operator (inherits `emili::Perturbation`)

Copy and rename these files as a starting point for new problems.

---

## How to Add a New Problem

### Step 1: Create Directory Structure

```bash
mkdir newproblem
```

### Step 2: Implement Problem Classes

Create `newproblem.h` with classes inheriting from base types:

```cpp
#include "../emilibase.h"

namespace emili {
namespace newproblem {

class NewProblem : public emili::Problem {
public:
    virtual double calcObjectiveFunctionValue(Solution& solution);
    virtual double evaluateSolution(Solution& solution);
    virtual int problemSize();
};

class NewSolution : public emili::Solution {
protected:
    virtual const void* getRawData() const;
    virtual void setRawData(const void* data);
public:
    virtual std::string getSolutionRepresentation();
    virtual Solution* clone();
};

// ... InitialSolution, Neighborhood classes

} // namespace newproblem
} // namespace emili
```

### Step 3: Create Builder

Create `newproblemBuilder.h`:

```cpp
#include "../generalParser.h"

namespace prs {

class NewProblemBuilder : public Builder {
public:
    NewProblemBuilder(GeneralParserE& gp, TokenManager& tm)
        : Builder(gp, tm) {}

    virtual bool isCompatibleWith(char* problem_definition) {
        return strcmp(problem_definition, "NEWPROBLEM") == 0;
    }

    virtual bool canOpenInstance(char* problem_definition) {
        return isCompatibleWith(problem_definition);
    }

    virtual emili::Problem* openInstance();
    virtual emili::Neighborhood* buildNeighborhood();
    virtual emili::InitialSolution* buildInitialSolution();
};

} // namespace prs
```

### Step 4: Update CMakeLists.txt

Add directory to source list:

```cmake
aux_source_directory(./newproblem SRC_LIST)
```

### Step 5: Register Builder in main.cpp

```cpp
#include "newproblem/newproblemBuilder.h"

int main(int argc, char* argv[]) {
    prs::GeneralParserE ps(argv, argc);
    prs::EmBaseBuilder emb(ps, ps.getTokenManager());
    prs::NewProblemBuilder npb(ps, ps.getTokenManager());  // Add this

    ps.addBuilder(&emb);
    ps.addBuilder(&npb);  // Add this

    LocalSearch* ls = ps.parseParams();
    // ...
}
```

---

## How to Add a New Algorithm

### Step 1: Define Token

In your builder's `.cpp` file:

```cpp
#define MY_ALGO "myalgo"
```

### Step 2: Implement Algorithm Class

Inherit from appropriate base class:

```cpp
class MyAlgorithm : public emili::LocalSearch {
public:
    MyAlgorithm(InitialSolution& init, Termination& term, Neighborhood& neigh)
        : LocalSearch(init, term, neigh) {}

    virtual Solution* search(Solution* initial);
};
```

### Step 3: Add Token Recognition

In `buildAlgo()`:

```cpp
emili::LocalSearch* NewProblemBuilder::buildAlgo() {
    emili::LocalSearch* ls = nullptr;

    if (tm.checkToken(MY_ALGO)) {
        printTab("My Algorithm");
        emili::InitialSolution* init = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR);
        emili::Termination* term = retrieveComponent(COMPONENT_TERMINATION_CRITERION);
        emili::Neighborhood* neigh = retrieveComponent(COMPONENT_NEIGHBORHOOD);
        ls = new emili::newproblem::MyAlgorithm(*init, *term, *neigh);
    }

    return ls;
}
```

**Usage:** `./emili instance.txt NEWPROBLEM myalgo random locmin myneigh`

---

## How to Add a New Neighborhood

### Step 1: Define Token

```cpp
#define NEIGHBORHOOD_MYNEW "mynew"
```

### Step 2: Implement Neighborhood Class

Required methods:

```cpp
class MyNewNeighborhood : public emili::Neighborhood {
protected:
    int current_i, current_j;  // Move indices

    virtual Solution* computeStep(Solution* step) {
        // Apply move (i, j) to solution, return modified solution
    }

    virtual void reverseLastMove(Solution* step) {
        // Undo the last move applied
    }

public:
    virtual void reset() {
        current_i = 0;
        current_j = 1;
    }

    virtual Solution* random(Solution* currentSolution) {
        // Return a random neighbor (new Solution object)
    }

    virtual int size() {
        // Return neighborhood size (e.g., n*(n-1)/2 for exchange)
    }

    virtual NeighborhoodIterator begin(Solution* base) {
        reset();
        return NeighborhoodIterator(this, base);
    }
};
```

### Step 3: Add Token Recognition

In `buildNeighborhood()`:

```cpp
emili::Neighborhood* NewProblemBuilder::buildNeighborhood() {
    emili::Neighborhood* neigh = nullptr;

    if (tm.checkToken(NEIGHBORHOOD_MYNEW)) {
        printTab("My New Neighborhood");
        neigh = new emili::newproblem::MyNewNeighborhood(*instance);
    }

    return neigh;
}
```

---

## Other Components Quick Reference

### Initial Solution

```cpp
class MyInitialSolution : public emili::InitialSolution {
public:
    MyInitialSolution(Problem& p) : InitialSolution(p) {}
    virtual Solution* generateSolution();       // Build initial solution
    virtual Solution* generateEmptySolution();  // Allocate empty solution
};
```

Token recognition in `buildInitialSolution()`:
```cpp
if (tm.checkToken("myinit")) {
    init = new MyInitialSolution(*instance);
}
```

### Termination Criterion

```cpp
class MyTermination : public emili::Termination {
public:
    virtual bool terminate(Solution* best, Solution* current);
    virtual void reset();
};
```

Common built-in terminations: `locmin`, `time`, `steps`

### Perturbation

```cpp
class MyPerturbation : public emili::Perturbation {
public:
    virtual Solution* perturb(Solution* solution);
};
```

Token recognition in `buildPerturbation()`.

### Acceptance Criterion

```cpp
class MyAcceptance : public emili::Acceptance {
public:
    virtual Solution* accept(Solution* candidate, Solution* current);
};
```

Common built-in: `improve` (only better), `always`, `pmetro` (probabilistic Metropolis)

---

## Key Files Reference

| File | Purpose |
|------|---------|
| `emilibase.h` | Base class definitions (Problem, Solution, Neighborhood, LocalSearch, etc.) |
| `generalParser.h` | TokenManager, Builder, GeneralParserE, Component classes |
| `generalParser.cpp` | EmBaseBuilder implementation (ILS, Tabu, VND, etc.) |
| `pfsp/pfspBuilder.cpp` | PFSP builder - comprehensive example of problem-specific builder |
| `pfsp/permutationflowshop.h` | PFSP classes - example problem implementation |
| `template/` | Starter templates for new problems |
| `main.cpp` | Builder registration example |
| `CMakeLists.txt` | Build configuration - add new source directories here |

---

## Tips

1. **Start with the template**: Copy `template/` files and rename for your problem
2. **Study PFSP**: The `pfsp/` directory is a complete, well-tested reference
3. **Use `printTab()`**: Helps debug parsing by printing matched tokens
4. **Component retrieval**: Use `retrieveComponent()` to get sub-components - it cycles through all compatible builders
5. **Token ordering**: Tokens are consumed in order; design your grammar accordingly
