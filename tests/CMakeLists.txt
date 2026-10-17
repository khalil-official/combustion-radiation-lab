# Numerical property and regression tests for radlib.
#
# The tests link the shared radlib library; on Windows the test executables
# are placed next to radlib.dll so the loader finds it without installation.

add_executable(test_rad_properties test_rad_properties.cc)
target_link_libraries(test_rad_properties PRIVATE radlib)
set_target_properties(test_rad_properties PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY $<TARGET_FILE_DIR:radlib>)
add_test(NAME rad_properties COMMAND test_rad_properties)

# Error-path check: an out-of-range band must terminate with a failure exit
# code (the library historically called exit(0), reporting failure as
# success). WILL_FAIL passes only when the command exits nonzero.
add_executable(test_invalid_band test_invalid_band.cc)
target_link_libraries(test_invalid_band PRIVATE radlib)
set_target_properties(test_invalid_band PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY $<TARGET_FILE_DIR:radlib>)
add_test(NAME invalid_band_reports_failure COMMAND test_invalid_band)
set_tests_properties(invalid_band_reports_failure PROPERTIES WILL_FAIL TRUE)
