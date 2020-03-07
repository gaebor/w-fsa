#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <string>

#include "ArgParser.h"
#include "Transducer.h"
#include "Utils.h"
#include "Isatty.h"

int main(int argc, const char** argv)
{
    // @see http://juditacs.github.io/2015/11/26/wordcount.html
    if (!IsTty(stdout))
        std::ios_base::sync_with_stdio(false);

    std::string transducer_filename, input_filename, output_filename;
    bool binary = false;
    std::string dump;
    {
        arg::Parser<> parser("AT&T Optimized Lookup", { "-h", "--help" }, std::cout, std::cerr, "", 80);

        parser.AddArg(transducer_filename, {}, "AT&T (text) format transducer filename", "filename");
        parser.AddFlag(binary, { "-b", "--binary" }, "read ATTOL binary format instead of text");
        parser.AddArg(dump, { "-d", "--dump" }, "after reading, dump the transducer in ATTOL binary format\ndon't perform actual lookup", "filename");
        parser.AddArg(input_filename, { "-i", "--input" }, "input file to analyze, stdin if empty", "filename");
        parser.AddArg(output_filename, { "-o", "--output" }, "output file, stdout if empty", "filename");
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
    std::string word;

    FILE* input = input_filename.empty() ? stdin : fopen(input_filename.c_str(), "r");
    if (!input)
        throw MyError("Unable to open input file \"", input_filename, "\"!");
    
    FILE* output = output_filename.empty() ? stdout : fopen(output_filename.c_str(), "w");
    if (!output)
        throw MyError("Unable to open output file \"", output_filename, "\"!");

    Transducer::ResultHandler resulthandler = [&](const Transducer::Path& path)
    {
        fputs(word.c_str(), output);
        putc(' ', output);
        for (auto i : path)
        {
            fprintf(output, "%u ", i);
        }
        putc('\n', output);
    };

    int c = 0;
    while (c != EOF)
    {
        word.clear();
        while ((c = fgetc(input)) != EOF && c != '\n' && c != '\0')
            word.push_back(static_cast<char>(c));

        if (!t.Lookup(word.c_str(), resulthandler))
        {
            fputs(word.c_str(), output);
            fputs(" ?", output);
        }
        putc('\n', output);
    }

    return 0;
}
catch (const MyError& e)
{
    std::cerr << e.what() << std::endl;
    return 1;
}
}
