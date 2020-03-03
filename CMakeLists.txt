cmake_minimum_required (VERSION 3.1)
project (wfsa)

set (CMAKE_CXX_STANDARD 11)

if(UNIX)
    add_compile_options( -m64 -Wall )
elseif(WIN32)
    add_compile_options( /W4 /favor:INTEL64 )
endif()

set (srcdir ${PROJECT_SOURCE_DIR}/src)
set (incdir ${PROJECT_SOURCE_DIR}/inc)

# add_executable(wfsa ${srcdir}/main.cpp ${incdir}/ArgParser.h
# 					${srcdir}/Corpus.cpp ${incdir}/Corpus.h
# 					${srcdir}/Utils.cpp ${incdir}/Utils.h
# 					${srcdir}/Isatty.cpp ${incdir}/Isatty.h
# 					${srcdir}/Fsa.cpp ${incdir}/Fsa.h)

add_executable(attol ${srcdir}/attol.cpp
					${srcdir}/Utils.cpp ${incdir}/Utils.h
					${srcdir}/FlagDiacritics.cpp ${incdir}/FlagDiacritics.h
					${srcdir}/Transducer.cpp ${incdir}/Transducer.h)

include_directories("${incdir}")
