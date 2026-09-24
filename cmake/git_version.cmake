# Writes ${BIN_DIR}/git_version.h from git_version.h.in with the current HEAD
# hash. Included at configure time and run again, via the git_version target,
# on every build so the 'commit :' line EMILI prints never goes stale.
# configure_file only touches the output when its content changes, so an
# unchanged hash does not trigger a recompile.
execute_process(COMMAND git rev-parse HEAD
                WORKING_DIRECTORY ${SRC_DIR}
                OUTPUT_VARIABLE GIT_COMMIT_HASH
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET)
if(NOT GIT_COMMIT_HASH)
    set(GIT_COMMIT_HASH "unknown")
endif()
configure_file(${SRC_DIR}/git_version.h.in ${BIN_DIR}/git_version.h @ONLY)
