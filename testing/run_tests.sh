#!/bin/bash

# EMILI Testing Framework - Test Runner Script
# Phase 2: Configuration Coverage

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Test configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/../build"
EMILI_BIN="${BUILD_DIR}/emili"
TEST_INSTANCE="${SCRIPT_DIR}/DD_Ta055.txt"
CONFIG_FILE="${SCRIPT_DIR}/configurations.txt"
TIME_LIMIT="0.01"  # 0.01 seconds as specified
LOG_DIR="${SCRIPT_DIR}/test_logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_TIMEOUT=0
TESTS_CRASHED=0

# Test modes
TEST_MODE="${1:-all}"  # all, quick, single
TEST_INDEX="${2:-1}"   # For single mode

# Function to print colored messages
print_pass() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

print_fail() {
    echo -e "${RED}[FAIL]${NC} $1"
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_section() {
    echo -e "${CYAN}=== $1 ===${NC}"
}

# Function to check prerequisites
check_prerequisites() {
    print_info "Checking prerequisites..."
    
    # Check if emili binary exists
    if [ ! -f "${EMILI_BIN}" ]; then
        print_fail "EMILI binary not found at ${EMILI_BIN}"
        print_info "Please build EMILI first:"
        print_info "  mkdir -p build && cd build"
        print_info "  cmake ../ && make"
        exit 1
    fi
    
    # Check if test instance exists
    if [ ! -f "${TEST_INSTANCE}" ]; then
        print_fail "Test instance not found at ${TEST_INSTANCE}"
        exit 1
    fi
    
    # Check if configuration file exists
    if [ ! -f "${CONFIG_FILE}" ]; then
        print_fail "Configuration file not found at ${CONFIG_FILE}"
        exit 1
    fi
    
    # Create log directory if it doesn't exist
    mkdir -p "${LOG_DIR}"
    
    print_pass "All prerequisites met"
}

# Function to detect timeout command
detect_timeout_command() {
    if command -v timeout >/dev/null 2>&1; then
        echo "timeout"
    elif command -v gtimeout >/dev/null 2>&1; then
        echo "gtimeout"
    else
        echo ""
    fi
}

# Function to run a single test configuration
run_single_test() {
    local config="$1"
    local test_name="$2"
    local config_index="$3"
    
    # Skip empty lines and comments
    if [[ -z "$config" ]] || [[ "$config" =~ ^# ]]; then
        return 0
    fi
    
    ((TESTS_RUN++))
    
    # Extract problem type from configuration
    local problem_type=$(echo "$config" | awk '{print $1}')
    
    print_info "Test ${config_index}: ${problem_type} configuration"
    if [ ${#config} -gt 80 ]; then
        echo "  Config: ${config:0:77}..."
    else
        echo "  Config: ${config}"
    fi
    
    # Create temporary files for output
    local stdout_file="${LOG_DIR}/test_${config_index}_stdout_${TIMESTAMP}.log"
    local stderr_file="${LOG_DIR}/test_${config_index}_stderr_${TIMESTAMP}.log"
    
    # Detect timeout command
    local timeout_cmd=$(detect_timeout_command)
    
    # Run the test with timeout
    local start_time=$(date +%s)
    
    if [ -n "${timeout_cmd}" ]; then
        ${timeout_cmd} 10s "${EMILI_BIN}" "${TEST_INSTANCE}" ${config} -ro "${TIME_LIMIT}" \
            > "${stdout_file}" 2> "${stderr_file}"
        local exit_code=$?
    else
        # Run without timeout but with background process monitoring
        "${EMILI_BIN}" "${TEST_INSTANCE}" ${config} -ro "${TIME_LIMIT}" \
            > "${stdout_file}" 2> "${stderr_file}" &
        local pid=$!
        
        # Wait for process with manual timeout
        local count=0
        while kill -0 $pid 2>/dev/null && [ $count -lt 100 ]; do
            sleep 0.1
            ((count++))
        done
        
        if kill -0 $pid 2>/dev/null; then
            kill -9 $pid 2>/dev/null
            wait $pid 2>/dev/null
            local exit_code=124  # Timeout exit code
        else
            wait $pid
            local exit_code=$?
        fi
    fi
    
    local end_time=$(date +%s)
    local elapsed_time=$((($end_time - $start_time) * 1000))  # Convert to milliseconds
    
    # Analyze results
    local status="UNKNOWN"
    local status_color="${YELLOW}"
    
    if [ ${exit_code} -eq 0 ]; then
        status="PASS"
        status_color="${GREEN}"
        ((TESTS_PASSED++))
        
        # Extract solution quality if available
        if grep -q "Best solution:" "${stdout_file}"; then
            local best_solution=$(grep "Best solution:" "${stdout_file}" | tail -1 | sed 's/Best solution://g' | xargs)
            echo -e "  ${status_color}✓${NC} Completed in ${elapsed_time}ms (Solution: ${best_solution})"
        else
            echo -e "  ${status_color}✓${NC} Completed in ${elapsed_time}ms"
        fi
        
        # Clean up successful test logs to save space
        rm -f "${stdout_file}" "${stderr_file}"
        
    elif [ ${exit_code} -eq 124 ]; then
        status="TIMEOUT"
        status_color="${YELLOW}"
        ((TESTS_TIMEOUT++))
        ((TESTS_FAILED++))
        echo -e "  ${status_color}⏱${NC} Timed out after 10 seconds"
        print_warning "  Logs saved to: ${stderr_file}"
        
    elif [ ${exit_code} -eq 139 ] || [ ${exit_code} -eq 134 ]; then
        status="CRASH"
        status_color="${RED}"
        ((TESTS_CRASHED++))
        ((TESTS_FAILED++))
        echo -e "  ${status_color}✗${NC} Segmentation fault or crash (exit code: ${exit_code})"
        
        # Show first few lines of error if any
        if [ -s "${stderr_file}" ]; then
            print_error "  Error output:"
            head -3 "${stderr_file}" | sed 's/^/    /'
            print_warning "  Full logs saved to: ${stderr_file}"
        fi
        
    else
        status="FAIL"
        status_color="${RED}"
        ((TESTS_FAILED++))
        echo -e "  ${status_color}✗${NC} Failed with exit code ${exit_code}"
        
        # Show error output if any
        if [ -s "${stderr_file}" ]; then
            print_error "  Error output:"
            head -5 "${stderr_file}" | sed 's/^/    /'
            print_warning "  Full logs saved to: ${stderr_file}"
        fi
    fi
    
    echo ""  # Add blank line for readability
}

# Function to run all configurations (Phase 2.1)
run_all_configurations() {
    print_section "Phase 2.1: Multi-Configuration Testing"
    
    local config_count=$(grep -v '^#' "${CONFIG_FILE}" | grep -v '^$' | wc -l | xargs)
    print_info "Found ${config_count} configurations to test"
    echo ""
    
    local index=0
    while IFS= read -r config; do
        # Skip empty lines and comments
        if [[ -z "$config" ]] || [[ "$config" =~ ^# ]]; then
            continue
        fi
        
        ((index++))
        run_single_test "${config}" "Configuration ${index}" "${index}"
        
        # Progress indicator
        local progress=$((index * 100 / config_count))
        echo -e "${CYAN}Progress: ${index}/${config_count} (${progress}%)${NC}"
        echo "----------------------------------------"
        echo ""
    done < "${CONFIG_FILE}"
}

# Function to run quick test (first config only)
run_quick_test() {
    print_section "Quick Test Mode"
    
    # Read first configuration from file
    local first_config=$(grep -v '^#' "${CONFIG_FILE}" | grep -v '^$' | head -1)
    
    if [ -z "${first_config}" ]; then
        print_fail "No valid configuration found in ${CONFIG_FILE}"
        exit 1
    fi
    
    run_single_test "${first_config}" "Quick Test" "1"
}

# Function to run single configuration by index
run_single_configuration() {
    local target_index=$1
    print_section "Single Configuration Test (Index: ${target_index})"
    
    local index=0
    while IFS= read -r config; do
        # Skip empty lines and comments
        if [[ -z "$config" ]] || [[ "$config" =~ ^# ]]; then
            continue
        fi
        
        ((index++))
        if [ ${index} -eq ${target_index} ]; then
            run_single_test "${config}" "Configuration ${index}" "${index}"
            return
        fi
    done < "${CONFIG_FILE}"
    
    print_fail "Configuration index ${target_index} not found"
    exit 1
}

# Function to clean old logs
clean_old_logs() {
    if [ -d "${LOG_DIR}" ]; then
        # Remove logs older than 7 days
        find "${LOG_DIR}" -name "*.log" -mtime +7 -delete 2>/dev/null
        
        # Keep only last 100 log files
        local log_count=$(find "${LOG_DIR}" -name "*.log" 2>/dev/null | wc -l)
        if [ ${log_count} -gt 100 ]; then
            find "${LOG_DIR}" -name "*.log" -print0 2>/dev/null | \
                xargs -0 ls -t | tail -n +101 | xargs rm -f 2>/dev/null
        fi
    fi
}

# Function to print detailed summary
print_summary() {
    echo ""
    echo "======================================"
    echo "           Test Summary               "
    echo "======================================"
    
    print_info "Tests run:     ${TESTS_RUN}"
    
    if [ ${TESTS_PASSED} -gt 0 ]; then
        print_pass "Tests passed:  ${TESTS_PASSED}"
    else
        print_info "Tests passed:  ${TESTS_PASSED}"
    fi
    
    if [ ${TESTS_FAILED} -gt 0 ]; then
        print_fail "Tests failed:  ${TESTS_FAILED}"
        
        if [ ${TESTS_TIMEOUT} -gt 0 ]; then
            print_warning "  - Timeouts:  ${TESTS_TIMEOUT}"
        fi
        
        if [ ${TESTS_CRASHED} -gt 0 ]; then
            print_error "  - Crashes:   ${TESTS_CRASHED}"
        fi
    else
        print_info "Tests failed:  ${TESTS_FAILED}"
    fi
    
    # Calculate pass rate
    if [ ${TESTS_RUN} -gt 0 ]; then
        local pass_rate=$((TESTS_PASSED * 100 / TESTS_RUN))
        echo ""
        print_info "Pass rate:     ${pass_rate}%"
    fi
    
    # Show log directory if there are failures
    if [ ${TESTS_FAILED} -gt 0 ]; then
        echo ""
        print_warning "Failed test logs saved in: ${LOG_DIR}"
    fi
    
    echo ""
}

# Main test execution
main() {
    echo "======================================"
    echo "  EMILI Testing Framework - Phase 2  "
    echo "======================================"
    echo ""
    
    # Check prerequisites
    check_prerequisites
    echo ""
    
    # Clean old logs
    clean_old_logs
    
    # Check for timeout command availability
    local timeout_cmd=$(detect_timeout_command)
    if [ -z "${timeout_cmd}" ]; then
        print_warning "timeout command not found, using fallback timeout mechanism"
        echo ""
    fi
    
    # Run tests based on mode
    case "${TEST_MODE}" in
        quick)
            run_quick_test
            ;;
        single)
            run_single_configuration "${TEST_INDEX}"
            ;;
        all|*)
            run_all_configurations
            ;;
    esac
    
    # Print summary
    print_summary
    
    # Exit with appropriate code
    if [ ${TESTS_FAILED} -gt 0 ]; then
        if [ ${TESTS_CRASHED} -gt 0 ]; then
            print_fail "Critical: Some tests crashed!"
            exit 2
        else
            print_fail "Some tests failed!"
            exit 1
        fi
    else
        print_pass "All tests passed!"
        exit 0
    fi
}

# Run main function
main "$@"