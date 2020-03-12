#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>

#include "ArgParser.h"
#include "atto/Transducer.h"
#include "Utils.h"
#include "Isatty.h"

int main(int argc, const char** argv)
{
    // @see http://juditacs.github.io/2015/11/26/wordcount.html
    if (!IsTty(stdout))
        std::ios_base::sync_with_stdio(false);

    std::string transducer_filename, input_filename, output_filename;
    double time_limit = 0.0;
    size_t max_depth = 0, max_results = 0;
    int print = 1;

    int flag_strategy = atto::Transducer::OBEY;
    std::string dump;
    {
        arg::Parser<> parser("AT&T Optimized Lookup", { "-h", "--help" }, std::cout, std::cerr, "", 80);

        parser.AddArg(transducer_filename, {}, 
                        "AT&T (text) format transducer filename", "filename");
        
        parser.AddArg(dump, { "-w", "--write", "--dump" }, 
                        "after reading, dump the transducer in ATTOL binary format\n"
                        "don't perform actual lookup", "filename");
        
        parser.AddArg(input_filename, { "-i", "--input" },
                        "input file to analyze, stdin if empty", "filename");
        parser.AddArg(output_filename, { "-o", "--output" },
                        "output file, stdout if empty", "filename");
        
        parser.AddArg(time_limit, { "-t", "--time" }, 
                        "time limit (in seconds) when not to search further\n"
                        "unlimited if set to 0");
        parser.AddArg(max_results, { "-n" },
                        "max number of results for one word\n"
                        "unlimited if set to 0");
        parser.AddArg(max_depth, { "-d", "--depth" },
                        "maximum depth to go down during lookup\n"
                        "unlimited if set to 0");
        
        parser.AddArg(flag_strategy, { "-f", "--flag" }, 
            ToStr("how to treat the flag diacritics\n",
                atto::Transducer::IGNORE, ": ignore (off)\n",
                atto::Transducer::OBEY, ": obey\n", 
                atto::Transducer::NEGATIVE, ": negative, return only those paths that were invalid flag-wise but correct analysis otherwise."),
            "", std::vector<int>({ atto::Transducer::OBEY, atto::Transducer::IGNORE, atto::Transducer::NEGATIVE }));
        
        parser.AddArg(print, { "-p", "--print" },
            ToStr("What to print about the analyses\n",
                0, ": output tape result\n",
                1, ": output tape result with weights\n",
                2, ": transition IDs along the path\n",
                3, ": transition IDs along the path with weights"),
            "", std::vector<int>({ 0,1,2,3 }));
  
        parser.Do(argc, argv);
    }
try{ 
    atto::Transducer t;
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

    std::string word;

    FILE* input = input_filename.empty() ? stdin : fopen(input_filename.c_str(), "r");
    if (!input)
        throw MyError("Unable to open input file \"", input_filename, "\"!");
    
    FILE* output = output_filename.empty() ? stdout : fopen(output_filename.c_str(), "w");
    if (!output)
        throw MyError("Unable to open output file \"", output_filename, "\"!");
    
    t.max_depth = max_depth;
    t.max_results = max_results;
    t.time_limit = time_limit;

    std::vector<float> weights;
    bool has_analysis;
    switch (print)
    {
    case 1:
        t.resulthandler = [&](const atto::Transducer::Path& path)
        {
            float weight = 0;
            fputs(word.c_str(), output);
            putc('\t', output);
            for (auto i : path)
            {
                if (!atto::FlagDiacritics<>::IsIt(std::get<4>(i)))
                    fputs(std::get<4>(i), output);
                weight += std::get<1>(i);
            }
            fprintf(output, "\t%g\n", weight);
            has_analysis = true;
        };
        break;
    case 2:
        t.resulthandler = [&](const atto::Transducer::Path& path)
        {
            fputs(word.c_str(), output);
            for (auto i : path)
            {
                fprintf(output, " %u", std::get<0>(i) + 1);
            }
            putc('\n', output);
            has_analysis = true;
        };
        break;
    case 3:
        t.resulthandler = [&](const atto::Transducer::Path& path)
        {
            float weight = 0;
            fputs(word.c_str(), output);
            for (auto i : path)
            {
                fprintf(output, " %u", std::get<0>(i) + 1);
                weight += std::get<1>(i);
            }
            fprintf(output, "\t%g\n", weight);
            has_analysis = true;
        };
        break;
    default:
        t.resulthandler = [&](const atto::Transducer::Path& path)
        {
            fputs(word.c_str(), output);
            putc('\t', output);
            for (auto i : path)
            {
                if (!atto::FlagDiacritics<>::IsIt(std::get<4>(i)))
                    fputs(std::get<4>(i), output);
            }
            putc('\n', output);
            has_analysis = true;
        };
        break;
    }

    int c = ~EOF;
    while (c != EOF)
    {
        word.clear();
        while ((c = fgetc(input)) != EOF && c != '\n' && c != '\0')
            word.push_back(static_cast<char>(c));
        has_analysis = false;
        t.Lookup(word.c_str(), atto::Transducer::FlagStrategy(flag_strategy));
        if (!has_analysis)
        {
            fprintf(output, "%s\t%s+\n", word.c_str(), word.c_str());
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
