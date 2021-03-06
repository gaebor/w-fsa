#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <memory>
#include <float.h>

#include "mkl.h"

#include "Test.h"

#include "Fsa.h"
#include "HessianLearner.h"
#include "QuasiNewtonLearner.h"
#include "Corpus.h"
#include "Utils.h"
#include "Recognize.h"

#include "ArgParser.h"

int main(int argc, const char* argv[])
{
    std::string automaton_filename, corpus_filename;
    std::string matrices_filename;
    std::string output_filename;
    std::string optimizer = "Hessian";
    int epochs = 20, initflags = 0;
    bool normalize = false, print = false, suppress = false, evaluate = false;
    bool print_recognize = false;
    double eta = 1.0;
    int recognize = 0;
    bool initx = false;
    double tolerance = 1e-6;
    int threads = 0;
    const int max_threads = mkl_get_max_threads();
    size_t test_n = 0, test_t = 0;
#if defined MKL_ILP64 && defined mkl_set_interface_layer
    if (mkl_set_interface_layer(MKL_INTERFACE_ILP64) == -1)
    {
       std::cerr << "Unable to specify ILP64 interface layer!" << std::endl;
       return 1;
    }
#endif
    {
        MKLVersion Version;
        mkl_get_version(&Version);
        std::ostringstream oss;
        oss << "\n --- Put the W in FSA --- \n\n"
            "Learning weights of a finite state automaton\n"
            "Author: Gabor Borbely, Contact: borbely@math.bme.hu, License: MIT\n\n" <<
            "Using Intel(R) Math Kernel Library Version " << Version.MajorVersion <<
            '.' << Version.MinorVersion << '.' << Version.UpdateVersion <<
            ' ' << Version.ProductStatus << " Build " << Version.Build <<
            "\nfor " << Version.Processor << "\n";
        oss << "MKL_INT = " << sizeof(MKL_INT) << "Bytes\n";
        if (max_threads == 1)
            oss << "Only 1 thread";
        else
            oss << "Maximum " << max_threads << " threads";
        oss << " supported!\n\n"
            "Environment variables to control MKL:\n"
            "MKL_THREADING_LAYER=[INTEL|GNU|TBB|SEQUENTIAL]\n"
            "MKL_NUM_THREADS=int\n"
            "MKL_VML_MODE=[VML_HA|VML_LA|VML_EP] meaning 'high accuracy', 'low accuracy' and 'enhanced performance' respectively";

        arg::Parser<> parser(oss.str(), { "-h", "--help" }, std::cout, std::cerr, "", 80);

        parser.AddArg(corpus_filename, { "-c", "--corpus" }, "corpus to load");
        parser.AddArg(automaton_filename, { "-a", "--automaton", "--load" }, "FSA to load");
        parser.AddArg(output_filename, { "-o", "--output" }, "Write learned WFSA in this file, or stdout if empty");
        parser.AddArg(epochs, { "-e", "--epoch", "--epochs" }, "number of maximum optimization epochs");
        parser.AddArg(eta, { "-l", "--learning", "--eta" }, "learning rate");
        parser.AddArg(tolerance, { "-tol", "--tol", "--tolerance" }, "tolerance when to stop");
        parser.AddFlag(evaluate, { "-eval", "--eval", "--evaluate" }, "evaluate model after optimization");
        parser.AddFlag(normalize, { "-n", "--normalize" }, "normalize the automaton after optimization");
        parser.AddFlag(print, { "-p", "--print" }, "prints extra info to stderr");
        parser.AddFlag(suppress, { "-s", "--suppress" }, "suppresses printing of learned FSA to stdout");
        parser.AddArg(threads, { "-t", "--thread", "--threads" }, "sets the number of MKL threads,\nif zero then leave it to the environment variable");
        parser.AddArg(recognize, { "-r", "--recognize" }, "recognize algorithm, 0: Breadth-first search, 1: Depth-first search", "", std::vector<int>({ 0, 1 }));
        parser.AddArg(optimizer, { "-opt", "--optimizer" }, "", "", std::vector<std::string>({ "Hessian", "QuasiNewton" }));
        parser.AddFlag(print_recognize, { "-pr", "--print-recognize" }, "print recognized paths to stderr");
        parser.AddArg(matrices_filename, { "-m", "--matrices", "--matrix" }, "name for matrix files\n"
            "If filename starts with \"<\" then loads, if \">\" then saves.\n"
            "loading from matrix files implies uniform initialization, unless initial x vector is provided with \"-x\"");
        parser.AddFlag(initx, { "-x", "--initx", "--initial" }, "if set, reads initial x vector from stdin");
        parser.AddArg(initflags, { "-i", "--init" }, "bitfield to control initialization before optimization\n"
            "1: make weights total uniform\n"
            "2: normalize automaton\n"
            "4: initialize Lagrange multipliers\n"
            "8: use Hessian of objective, otherwise use only Hessian of constraints\n"
            "16: computes a permutation that minimizes the fill-in during the factorization phase\n"
            "32: exponential Lagrange multipliers\n"
        );
        parser.AddArg(test_n, { "-testn", "-test_n", "--testn", "--test_n" }, "bitfield which vector sizes to test, calculated in powers of 10:\n"
            "0: no tests, 1: 1, 2: 10, 4: 100, 2^n -> 10^n \n"
            "for the tests the '--epoch' argument is used to specify the number of repetitions."
        );
        parser.AddArg(test_t, { "-testt", "-test_t", "--testt", "--test_t" }, "bitfield for which tests to run\n"
            "if positive, then performs the tests and exits.\n"
            "1\tstd::exp\n"
            "2\tvdExp\n"
            "4\tstd::accumulate\n"
            "8\tcblas_ddot with ones vector\n"
            "16\tstd::fill with zeros\n"
            "32\tcblas_dscal with zero"
        );

        parser.Do(argc, argv);
    }
    if (threads > 0)
    {
       mkl_set_num_threads(threads);
    }
    if (test_t > 0)
    {
        test(std::cout, epochs, test_n, test_t);
        return 0;
    }
    try {

    Fsa fsa;
    std::unique_ptr<Learner> learner(optimizer == "Hessian" ? (Learner*)(new HessianLearner()) : (Learner*)(new QuasiNewtonLearner()));
    if (!matrices_filename.empty() && matrices_filename.front() == '<')
    {   // skip Corpus and Fsa
        std::cerr << "Loading matrices \"" << matrices_filename.substr(1) << "\" ... ";
        std::cerr.flush();
        if (!learner->LoadMatrices(matrices_filename.substr(1)))
            throw LearnerError("Unable to load Learner from \"", matrices_filename.substr(1), "\"");
        std::cerr << "done" << std::endl;
        std::cerr << "Info:" <<
            "\n\tparameters: " << learner->GetNumberOfParameters() <<
            "\n\tconstraints: " << learner->GetNumberOfConstraints() <<
            "\n\tstrings: " << learner->GetNumberOfStrings() <<
            "\n\tpaths: " << learner->GetNumberOfPaths() <<
            "\n\tcommon support: " << learner->GetCommonSupport() <<
            "\n\tunique paths: " << (learner->HasUniquePaths() ? "true" : "false") <<
            std::endl;
    }else
    {
        std::cerr << "Corpus: "; std::cerr.flush();
        Corpus corpus;
        if (FILE* f = fopen(corpus_filename.c_str(), "r"))
        {
            corpus.Read(f);
            fclose(f);
        }
        else
        {
            std::cerr << "\nUnable to open \"" << corpus_filename << "\"!" << std::endl;
            return 1;
        }
        std::cerr << "\n\tsize: " << corpus.size() <<
            "\n\tsum: " << corpus.Sum();
        corpus.Renormalize();
        std::cerr << ", renormalized to " << corpus.Sum() << std::endl;

        std::cerr << "Automaton: "; std::cerr.flush();

        if (FILE* f = fopen(automaton_filename.c_str(), "r"))
        {
            fsa.Read(f);
            fclose(f);
        }
        else
        {
            std::cerr << "\nUnable to open \"" << automaton_filename << "\"!" << std::endl;
            return 1;
        }

        std::cerr << "\n\tstates: " << fsa.GetNumberOfStates() <<
            "\n\ttransitions: " << fsa.GetNumberOfTransitions() <<
            "\n\temissions: " << fsa.GetNumberOfEmissions() <<
            "\n\tparameters: " << fsa.GetNumberOfParameters() <<
            "\n\tconstraints: " << fsa.GetNumberOfConstraints() <<
            "\n\tfree parameters: " << fsa.GetNumberOfFreeParameters() <<
            std::endl;

        if (print || print_recognize)
        {
            typedef std::pair<std::string, std::string> Path;
            Recognizer<Path> recognizer(fsa.GetEndState(),
                [](Path& history, const Fsa::Transitions::value_type& transition,
                    const Fsa::Emissions::value_type& emission) -> Path&
            {
                history.first += emission.str;
                history.second += " -> ";
                history.second += transition.next->first;
                if (emission.str[0])
                {
                    history.second += '"';
                    history.second += emission.str;
                    history.second += '"';
                }
                return history;
            }, [](const Path& path)
            {
                fprintf(stderr, "%s: %s\n", path.first.c_str(), path.second.c_str());
            });
            const auto& start_state = fsa.GetTransitionMtx().at(fsa.GetStartState());
            const auto empty_history = Path("", fsa.GetStartState());
            for (const auto& word : corpus)
                recognizer.Recognize(word.first.c_str(), start_state, empty_history, recognize == 0);
        }

        std::cerr << "Recognize: "; std::cerr.flush();
        learner->BuildFrom(fsa, corpus, recognize == 0);

        std::cerr << "\n\tstrings: " << learner->GetNumberOfStrings() <<
            "\n\tpaths: " << learner->GetNumberOfPaths() <<
            "\n\tcommon support: " << learner->GetCommonSupport() <<
            "\n\tunique paths: " << (learner->HasUniquePaths() ? "true" : "false") <<
            "\nAfter trimming:\n\tparameters: " << learner->GetNumberOfParameters() <<
            "\n\tconstraints: " << learner->GetNumberOfConstraints() <<
            std::endl;
    }

    if (learner->GetNumberOfParameters() == 0)
    {
        //TODO The automaton may generate one string deterministically!
        std::cerr << "Empty automaton!" << std::endl;
        return 1;
    }
    if (learner->GetNumberOfStrings() == 0)
    {
        //TODO aux model only!
        std::cerr << "Automaton cannot generate any of the strings!" << std::endl;
        return 1;
    }
    learner->Finalize();

    if (print)
    {
        fputs("C:\n", stderr);
        learner->PrintC(stderr);
        fputs("M:\n", stderr);
        learner->PrintM(stderr);
        fputs("P:\n", stderr);
        learner->PrintP(stderr);
    }

    if (!matrices_filename.empty() && matrices_filename.front() == '>')
    {
        std::cerr << "Saving matrices \"" << matrices_filename.substr(1) << "\" ... "; std::cerr.flush();
        if (!learner->SaveMatrices(matrices_filename.substr(1)))
            std::cerr << "Failed!";
        else
            std::cerr << "Done";
        std::cerr << std::endl;
    }

    std::cerr << "Initialize ... "; std::cerr.flush();
    if (initx)
    {
        std::vector<double> x;
        while (std::cin && x.size() < learner->GetNumberOfParameters())
        {
            x.emplace_back();
            std::cin >> x.back();
        };
        if (!std::cin)
            throw MyError("Cannot read initial x value!");
        learner->Init(initflags, x.data());
    }else
        learner->Init(initflags);

    std::cerr << "done" << std::endl;
    //std::cerr << "augmented Hessian:\n\trows: " << learner.GetNumberOfParameters() + learner.GetNumberOfConstraints() <<
    //    "\n\tnnz: " << learner.GetNnzHessian() <<
    //    "\n\tfill: " << learner.GetHessianFillRatio() << std::endl;
    //
    const auto width = (int)std::ceil(std::log10(epochs + 1));
    if (epochs > 0)
    {
        // Clock<> clock;
        std::cerr << "Optimization:" << std::endl;
        for (int e = 1; e <= epochs; ++e)
        {
            if (e % 20 == 1)
                    std::cerr << "epoch\t" << learner->GetOptimizationHeader() << std::endl;
            learner->OptimizationStep(eta, print);
            fprintf(stderr, "%0*d\t", width, e);
            // PrintFixedWidth(stderr, clock.Tock()/e, 9);
            // fputs("\t", stderr);
            const auto info = learner->GetOptimizationInfo();
            for (double x : info)
            {
                PrintFixedWidth(stderr, x, 9);
                fputs(" ", stderr);
            }

            std::cerr << std::endl;
            for (double x : info)
            {
                if (!std::isfinite(x))
                {
                    std::cerr << std::endl;
                    throw LearnerError(x, " detected at epoch ", e);
                }
            }

            if (learner->HaltCondition(tolerance))
                break;
        }
    }

    if (normalize)
        learner->Renormalize();

    if (evaluate)
    {
        const auto results = learner->GetOptimizationResult(print);
        std::cerr.precision(DBL_DIG);
        std::cerr << "Result:";
        for (double x : results)
            std::cerr << ' ' << x;
        std::cerr << std::endl;
    }
    if (!suppress)
    {
        FILE* outf = output_filename.empty() ? stdout : fopen(output_filename.c_str(), "w");
        if (outf)
        {
            if (fsa.GetNumberOfParameters() > 0)
            {// FSA has been loaded
                learner->RewriteWeights(fsa);
                fsa.Dump(outf);
            }
            else
            {// loaded from matrices
                const MKL_INT n = learner->GetNumberOfParameters();
                const MKL_INT k = learner->GetNumberOfConstraints();
                if (n > 0)
                    fprintf(outf, "%g", *(learner->GetWeights()));
                for (auto* x = learner->GetWeights() + 1; x < learner->GetWeights() + n; ++x)
                    fprintf(outf, " %g", *x);
                std::cout << std::endl;
                if (k > 0)
                    fprintf(outf, "%g", *(learner->GetWeights() + n));
                for (auto* x = learner->GetWeights() + n + 1; x < learner->GetWeights() + n + k; ++x)
                    fprintf(outf, " %g", *x);
                std::cout << std::endl;
            }
            if (outf != stdout)
                fclose(outf);
        }else
            throw MyError("Unable to open output file \"", output_filename, "\" for writing!");
    }
    }
    catch (std::exception& e)
    {
        fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    return 0;
}
