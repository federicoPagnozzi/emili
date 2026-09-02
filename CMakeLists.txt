cmake_minimum_required(VERSION 3.8)
project(emili)

# Commit id header (replaces the old in-tree sed on generalParser.cpp).
# Refreshed at cmake time; goldens normalize the commit line so staleness
# between reconfigures is harmless.
execute_process(COMMAND git rev-parse HEAD
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_COMMIT_HASH
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET)
if(NOT GIT_COMMIT_HASH)
    set(GIT_COMMIT_HASH "unknown")
endif()
configure_file(${PROJECT_SOURCE_DIR}/git_version.h.in
               ${PROJECT_BINARY_DIR}/git_version.h @ONLY)
include_directories(${PROJECT_BINARY_DIR})
aux_source_directory(. SRC_LIST)
aux_source_directory(./pfsp SRC_LIST)
#aux_source_directory(./template SRC_LIST)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11")
option(USE_SSE "Use SSE compile path if available" OFF)
option(USE_NEW_MAIN "Use new loading paradigm" ON)
option(WITH_STATS "Print solution informations" OFF)
option(IRACE_OPTIMISED_FLAGS "compiles EMILI using the optimization setting tuned by irace" OFF)
option(OPT_OPTIMISED_FLAGS "compiles EMILI using the optimization setting tuned by open tuner" OFF)
option(O3_FLAGS "compiles EMILI using the default -O3 flag" ON)
option(DEBUG_FLAGS "compiles EMILI using the debug flags" OFF)
option(USE_EM_LIB "experimental library based extension mechanism" OFF)
option(PORTABLE "Links the shared library in a static way" OFF)

#IRACE flags
if(IRACE_OPTIMISED_FLAGS)
	message(STATUS "Using irace tuned gcc settings")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -fno-auto-inc-dec -fno-branch-count-reg -fcombine-stack-adjustments -fcompare-elim -fcprop-registers -fdce -fno-defer-pop -fdelayed-branch -fno-dse -fno-forward-propagate -fguess-branch-probability -fno-if-conversion2 -fif-conversion -finline-functions-called-once -fipa-pure-const -fno-ipa-profile -fipa-reference -fmerge-constants -fmove-loop-invariants -fno-reorder-blocks -fshrink-wrap -fno-shrink-wrap-separate -fno-split-wide-types -fssa-backprop -fno-ssa-phiopt -fno-tree-ccp -ftree-ch -ftree-coalesce-vars -fno-tree-copy-prop -fno-tree-dce -fno-tree-dominator-opts -fno-tree-dse -ftree-forwprop -ftree-fre -ftree-phiprop -fno-tree-sink -ftree-slsr -ftree-sra -fno-tree-pta -ftree-ter -fomit-frame-pointer -fthread-jumps -falign-functions=0 -falign-jumps=0 -falign-loops=2 -falign-labels=1 -fcaller-saves -fcrossjumping -fno-cse-follow-jumps -fno-cse-skip-blocks -fno-delete-null-pointer-checks -fno-devirtualize -fdevirtualize-speculatively -fexpensive-optimizations -fgcse -fno-gcse-lm -fhoist-adjacent-loads -finline-small-functions -fno-indirect-inlining -fipa-cp -fno-ipa-bit-cp -fno-ipa-vrp -fno-ipa-sra -fno-ipa-icf -fno-lra-remat -foptimize-sibling-calls -fno-optimize-strlen -fpartial-inlining -fpeephole -fpeephole2 -freorder-blocks-algorithm=stc -fno-reorder-blocks-and-partition -fno-reorder-functions -fno-sched-interblock -fsched-spec -fschedule-insns -fschedule-insns2 -fno-store-merging -fstrict-aliasing -fstrict-overflow -fno-tree-builtin-call-dce -fno-tree-switch-conversion -ftree-tail-merge -fcode-hoisting -fno-tree-pre -fno-tree-vrp -fipa-ra -finline-functions -fno-unswitch-loops -fno-predictive-commoning -fno-gcse-after-reload -fno-tree-loop-vectorize -ftree-loop-distribute-patterns -fno-split-paths -fvect-cost-model=dynamic -ftree-partial-pre -fno-peel-loops -fno-ipa-cp-clone -ffp-contract=fast -finline -fearly-inlining -finline-limit=20 -fkeep-inline-functions -fno-keep-static-functions -fno-keep-static-consts -fno-merge-all-constants -fmodulo-sched -ffunction-cse -fzero-initialized-in-bss -fno-gcse-sm -fgcse-las -faggressive-loop-optimizations -fno-unconstrained-commons -fdeclone-ctor-dtor -fno-devirtualize-at-ltrans -free -flifetime-dse=2 -flive-range-shrinkage -fira-algorithm=priority -fira-region=one -fira-hoist-pressure -fno-ira-loop-pressure -fno-ira-share-save-slots -fno-ira-share-spill-slots -fsched-pressure -fsched-spec-load -fno-sched-spec-load-dangerous -fsched-stalled-insns=1 -fsched-stalled-insns-dep=4 -fno-sched2-use-superblocks -fsched-group-heuristic -fsched-critical-path-heuristic -fsched-spec-insn-heuristic -fno-sched-rank-heuristic -fno-sched-last-insn-heuristic -fno-sched-dep-count-heuristic -freschedule-modulo-scheduled-loops -fselective-scheduling -fno-selective-scheduling2 -fno-sel-sched-pipelining -fno-semantic-interposition -fconserve-stack -ftree-reassoc -fipa-pta -ftree-loop-optimize -fno-tree-loop-if-convert -ftree-loop-distribution -ftree-loop-im -fno-tree-loop-ivcanon -fivopts -ftree-vectorize -fsplit-ivs-in-unroller -fvariable-expansion-in-unroller -fno-prefetch-loop-arrays -fprintf-return-value -fno-limit-function-alignment -fno-web -fno-whole-program -fno-lto -fno-fat-lto-objects -funroll-loops -funroll-all-loops")
elseif(OPT_OPTIMISED_FLAGS)
	message(STATUS "Using open tuner gcc settings")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -fno-dce -fsemantic-interposition -fno-use-linker-plugin -fno-delete-null-pointer-checks -finline-limit=10 -ftree-builtin-call-dce -fno-sched-stalled-insns -fipa-profile -fno-tree-phiprop -fno-optimize-sibling-calls -fno-sel-sched-pipelining -fno-live-range-shrinkage -flto-compression-level=4 -fprintf-return-value -fno-gcse-after-reload -freorder-blocks-algorithm=stc -fno-reorder-functions -fipa-bit-cp -fno-fat-lto-objects -fsched-stalled-insns-dep=4 -ffp-contract=off -fdevirtualize-speculatively -ftree-slsr -fno-tree-pta -ftree-ccp -fira-share-save-slots -fno-tree-fre -fschedule-insns -fsched-pressure -fif-conversion2 -fmodulo-sched-allow-regmoves -fno-ipa-pta -fvariable-expansion-in-unroller -fno-optimize-strlen -fno-sched-spec-load-dangerous -flifetime-dse=0 -fno-delayed-branch -fkeep-static-functions -ftree-partial-pre -fssa-phiopt -faggressive-loop-optimizations -falign-labels=1 -fkeep-inline-functions -fno-prefetch-loop-arrays -fno-limit-function-alignment -fno-tree-sink -ftree-loop-optimize -fselective-scheduling2 -ffunction-cse -fno-split-paths -fno-shrink-wrap -flra-remat -fgcse-sm -fira-region=all -falign-loops=4 -fno-tree-reassoc -fira-hoist-pressure -findirect-inlining -fno-peel-loops -fomit-frame-pointer -flto-partition=none -fno-tree-tail-merge -fno-predictive-commoning -fcode-hoisting -fdevirtualize-at-ltrans -free -fno-store-merging -fno-sched-dep-count-heuristic -fno-tree-slp-vectorize -fcprop-registers -ftree-loop-distribution -fno-split-ivs-in-unroller -fno-sched-spec-load -fivopts -fipa-icf -fdefer-pop -fno-isolate-erroneous-paths-dereference -fstrict-overflow -fno-dse -ftree-switch-conversion -ftree-forwprop -fno-zero-initialized-in-bss -fcrossjumping -fdeclone-ctor-dtor -fno-inline-small-functions -fno-gcse-lm -ftree-dominator-opts -fsched-last-insn-heuristic -fsched-spec-insn-heuristic -frerun-cse-after-loop -fno-peephole -fgcse -fmove-loop-invariants -ftree-loop-im -ftree-ter -fno-branch-count-reg -fearly-inlining -ftree-ch -fno-tree-loop-if-convert -fmerge-constants -fno-sel-sched-pipelining-outer-loops -fira-algorithm=CB -fno-conserve-stack -fno-sched-interblock -ftree-coalesce-vars -ftree-vrp -fno-tree-loop-ivcanon -ftree-loop-distribute-patterns -funswitch-loops -fif-conversion -fno-compare-elim -fno-split-wide-types -fdevirtualize -fweb -fsched-spec -fno-modulo-sched -fauto-inc-dec -fno-ipa-pure-const -freschedule-modulo-scheduled-loops -fstrict-aliasing -fcaller-saves -fno-cse-skip-blocks -fcse-follow-jumps -fno-keep-static-consts -falign-jumps=64 -fcombine-stack-adjustments -fno-unconstrained-commons -fschedule-insns2 -fno-ipa-reference -flto-odr-type-merging -finline -fno-tree-copy-prop -fno-tree-dce -falign-functions=4 -fno-unroll-loops -fguess-branch-probability -fhoist-adjacent-loads -fno-lto -freorder-blocks -fno-reorder-blocks-and-partition -fno-expensive-optimizations -fipa-sra -fno-inline-functions -fno-sched2-use-superblocks -fno-tree-bit-ccp -fno-shrink-wrap-separate -fno-partial-inlining -fno-tree-vectorize -fno-ipa-ra -fno-ssa-backprop -fno-thread-jumps -funroll-all-loops -fno-inline-functions-called-once -fira-share-spill-slots -fno-forward-propagate -fipa-cp-clone -fno-tree-pre -fselective-scheduling -fno-tree-loop-vectorize -fno-merge-all-constants -fipa-cp -ftree-dse -fno-sched-rank-heuristic -fira-loop-pressure -fsched-group-heuristic -ftree-sra -fno-sched-critical-path-heuristic -fvect-cost-model=cheap -fno-peephole2 -fipa-vrp -fno-whole-program -fno-gcse-las")
elseif(O3_FLAGS)
	message(STATUS "Compiling with -O3")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")
elseif(DEBUG_FLAGS)
        message(STATUS "Compiling in debug mode with ASan/UBSan")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra -Werror")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,undefined")
endif(IRACE_OPTIMISED_FLAGS)

# P1.3 warning policy: -Wall on release, -Wall -Wextra -Werror on the sanitizer
# build. Three pervasive stylistic families in frozen PFSP/legacy code are
# suppressed as documented backlog (fix candidates for P2+):
#   -Wreorder-ctor      (~160 init-list orderings, mostly permutationflowshop.h)
#   -Wvla-cxx-extension (28 VLAs in Taillard-accelerated hot paths)
#   -Woverloaded-virtual(~20 name-hidings across the search hierarchy)
# -Wunused-parameter fires on the many intentionally-unused base-class hook
# parameters in emilibase.h, and -Wsign-compare on ~30 int-vs-size_t loop
# indices in frozen PFSP code (core sites are fixed); suppressed with the
# same backlog status.
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    add_compile_options(-Wno-reorder-ctor -Wno-vla-cxx-extension -Wno-overloaded-virtual -Wno-unused-parameter -Wno-sign-compare)
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    # -Wmaybe-uninitialized is a GCC-only flow heuristic that fires on ~20
    # 'min'/'ind'/'best' locals in frozen PFSP hot loops (permutationflowshop.cpp)
    # whose initialization the compiler cannot prove; clang has no equivalent
    # and ASan/UBSan runs are clean. Suppressed with the same backlog status
    # as the families above (P2 Task 7, Linux CI).
    add_compile_options(-Wno-reorder -Wno-vla -Wno-overloaded-virtual -Wno-unused-parameter -Wno-sign-compare -Wno-maybe-uninitialized)
endif()

#SSE
if(USE_SSE)
	message(STATUS "Using sse compile path") 
	add_definitions(-DENABLE_SSE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
else()
        message(STATUS "sse compile path disabled")
endif(USE_SSE)

#Solution stats
if(WITH_STATS)
    message(STATUS "Solution stats enabled") 
    add_definitions(-DWITH_STATS)
else()
    message(STATUS "Solution stats disabled")
endif(WITH_STATS)

if(USE_NEW_MAIN)
    message(STATUS "Component based algorithm loading enabled")
    add_definitions(-DMAIN_NEW)
else()
    message(STATUS "Component based algorithm loading disabled")
endif(USE_NEW_MAIN)


if(USE_EM_LIB)
    message(STATUS "Experimental library based extension mechanism enabled")
    add_definitions(-DEM_LIB)
else()
    message(STATUS "Experimental library based extension mechanism disabled")
endif(USE_EM_LIB)

if(PORTABLE)
	message(STATUS "Using static linking") 
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
else()
        message(STATUS "Using dynamic linking")
endif(PORTABLE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(${PROJECT_NAME} ${SRC_LIST})

# Testing support
enable_testing()

# Add custom target for running tests (using 'check' instead of 'test' to avoid conflicts)
add_custom_target(check
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testing/run_tests.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running EMILI test suite"
    VERBATIM
)

# Make check depend on the emili executable
add_dependencies(check ${PROJECT_NAME})

# Also add a CTest integration for standard 'make test' compatibility
add_test(NAME emili_basic_test
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testing/run_tests.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
