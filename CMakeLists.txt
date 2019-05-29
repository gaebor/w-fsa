cmake_minimum_required (VERSION 3.1)
project (wfsa)

link_directories($ENV{MKLROOT}/lib/intel64)

set (CMAKE_CXX_STANDARD 11)

if(UNIX)
    add_compile_options("-m64" "-Wall")
elseif(WIN32)
    add_compile_options("/W4" "/favor:INTEL64")
endif()

file(GLOB WFSA_SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB WFSA_INC ${PROJECT_SOURCE_DIR}/inc/*.h)

add_executable(wfsa ${WFSA_SRC} ${WFSA_INC})

include_directories("$ENV{MKLROOT}/include" "${PROJECT_SOURCE_DIR}/inc")

# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
add_definitions(-DMKL_ILP64)

target_link_libraries(wfsa mkl_intel_ilp64 mkl_sequential mkl_core)
# target_link_libraries(wfsa mkl_rt)

if (UNIX)
    target_link_libraries(wfsa pthread)
endif()

include(CTest)

add_test(testhelp wfsa -h)

add_test(test1 wfsa -a ../data/test.wfsa -c ../data/test.corpus  -n -eval -e 20 -i 31)
add_test(test2 wfsa -a ../data/test2.wfsa -c ../data/test.corpus -n -eval -e 20 -i 31)
add_test(test3 wfsa -a ../data/test3.wfsa -c ../data/test.corpus -n -eval -e 20 -i 31)
add_test(test4 wfsa -a ../data/test4.wfsa -c ../data/test.corpus -n -eval -e 20 -i 31)
add_test(test_loop wfsa -a ../data/test.loop.wfsa -c ../data/test.corpus -n -eval -e 20 -i 31)
add_test(test_talk wfsa -a ../data/talk.wfsa -c ../data/talk.corpus -n -eval -e 20 -i 31)
add_test(test_talk2 wfsa -a ../data/talk.wfsa -c ../data/test.corpus -n -eval -e 20 -i 31)

add_test(test1v wfsa -a ../data/test.wfsa -c ../data/test.corpus  -p -n -eval -e 20 -i 31)
add_test(test2v wfsa -a ../data/test2.wfsa -c ../data/test.corpus -p -n -eval -e 20 -i 31)
add_test(test3v wfsa -a ../data/test3.wfsa -c ../data/test.corpus -p -n -eval -e 20 -i 31)
add_test(test4v wfsa -a ../data/test4.wfsa -c ../data/test.corpus -p -n -eval -e 20 -i 31)
add_test(test_loopv wfsa -a ../data/test.loop.wfsa -c ../data/test.corpus -p -n -eval -e 20 -i 31)
add_test(test_talkv wfsa -a ../data/talk.wfsa -c ../data/talk.corpus -p -n -eval -e 20 -i 31)
add_test(test_talk2v wfsa -a ../data/talk.wfsa -c ../data/test.corpus -p -n -eval -e 20 -i 31)

set_tests_properties(test1 PROPERTIES WILL_FAIL TRUE)
set_tests_properties(test1v PROPERTIES WILL_FAIL TRUE)
