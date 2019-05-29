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
