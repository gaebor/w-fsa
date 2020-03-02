#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <string>

#include "ArgParser.h"
#include "Transducer.h"
#include "Utils.h"

int main(int argc, const char** argv)
{
    std::string transducer_filename;
    {
        arg::Parser<> parser("AT&T Optimized Lookup", { "-h", "--help" }, std::cout, std::cerr, "", 80);

        parser.AddArg(transducer_filename, {}, "AT&T format transducer filename", "filename");
        
        parser.Do(argc, argv);
    }
try{
    Transducer t;
    std::ifstream ifs(transducer_filename);
    if (ifs)
        t.Read(ifs);
    else
    {
        std::cerr << "Cannot read \"" << transducer_filename << "\"!" << std::endl;
        return 1;
    }

    std::cerr << "Number of states (including start): " << t.GetNumberOfStates() << std::endl;
    std::cerr << "Number of transitions (including finishing from final states): " << t.GetNumberOfTransitions() << std::endl;
    std::cerr << "Memory usage (approximate): " << t.GetAllocatedMemory() << " bytes" << std::endl;

    Transducer::ResultHandler resulthandler = [](const Transducer::Path& path)
    {
        for (auto i : path)
        {
            printf("%u ", i);
        }
        putc('\n', stdout);
    };

    std::string word;
    while (std::cin >> word)
    {
        t.Lookup(word.c_str(), resulthandler);
        putc('\n', stdout);
    }

    return 0;
}
catch (const MyError& e)
{
    std::cerr << e.what() << std::endl;
    return 1;
}
}
