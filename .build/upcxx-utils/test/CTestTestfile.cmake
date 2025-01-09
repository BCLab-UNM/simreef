# CMake generated Testfile for 
# Source directory: /users/mfricke/simreef/upcxx-utils/test
# Build directory: /users/mfricke/simreef/.build/upcxx-utils/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_version "upcxx-run" "-n" "2" "$CMAKE_CURRENT_BIN_DIR}/test_version")
add_test(test_upcxx_utils "upcxx-run" "-n" "2" "$CMAKE_CURRENT_BIN_DIR}/test_upcxx_utils")
add_test(test_log "upcxx-run" "-n" "2" "$CMAKE_CURRENT_BIN_DIR}/test_log")
add_test(test_timers "upcxx-run" "-n" "2" "$CMAKE_CURRENT_BIN_DIR}/test_timers")
add_test(test_progress_bar "upcxx-run" "-n" "2" "$CMAKE_CURRENT_BIN_DIR}/test_progress_bar")
add_test(test_aggr_store "upcxx-run" "-n" "2" "$CMAKE_CURRENT_BIN_DIR}/test_aggr_store")
