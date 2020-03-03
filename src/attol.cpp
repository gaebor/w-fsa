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
    bool binary = false;
    std::string dump;
    {
        arg::Parser<> parser("AT&T Optimized Lookup", { "-h", "--help" }, std::cout, std::cerr, "", 80);

        parser.AddArg(transducer_filename, {}, "AT&T (text) format transducer filename", "filename");
        parser.AddFlag(binary, { "-b", "--binary" }, "read ATTOL binary format instead of text");
        parser.AddArg(dump, { "-d", "--dump" }, "after reading, dump the transducer in ATTOL binary format\ndon't perform actual lookup", "filename");

        parser.Do(argc, argv);
    }
try{
    Transducer t;
    if (binary)
    {
        if (FILE* f = fopen(transducer_filename.c_str(), "rb"))
            t.ReadBinary(f);
        else
        {
            std::cerr << "Cannot open to read \"" << transducer_filename << "\"!" << std::endl;
            return 1;
        }
    } else
    {
        std::ifstream ifs(transducer_filename);
        if (ifs)
            t.Read(ifs);
        else
        {
            std::cerr << "Cannot open to read \"" << transducer_filename << "\"!" << std::endl;
            return 1;
        }
        std::cerr << "Number of states (including start): " << t.GetNumberOfStates() << std::endl;
        std::cerr << "Number of transitions (including finishing from final states): " << t.GetNumberOfTransitions() << std::endl;
    }

    if (!dump.empty())
    {
        if (FILE* f = fopen(dump.c_str(), "wb"))
        {
            t.DumpBinary(f);
            return 0;
        }
        else
        {
            std::cerr << "Cannot open to write \"" << dump << "\"!" << std::endl;
            return 1;
        }
    }
    bool has_analysis = false;
    std::string word;

    Transducer::ResultHandler resulthandler = [&](const Transducer::Path& path)
    {
        fputs(word.c_str(), stdout);
        putc(' ', stdout);
        for (auto i : path)
        {
            printf("%u ", i);
        }
        putc('\n', stdout);
        has_analysis = true;
    };

    while (std::cin >> word)
    {
        has_analysis = false;
        t.Lookup(word.c_str(), resulthandler);
        if (!has_analysis)   
        {
            fputs(word.c_str(), stdout);
            puts(" ?");
        }
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
