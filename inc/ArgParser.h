#pragma once

#include <typeinfo>
#include <string>
#include <cstring>
#include <vector>
#include <memory>
#include <unordered_set>
#include <iostream>
#include <sstream>
#include <algorithm>

namespace arg
{
    template<class String>
    struct StdStreams;

    template<class Allocator>
    struct StdStreams<std::basic_string<char, std::char_traits<char>, Allocator>>
    {
        static std::basic_istream<char, std::char_traits<char>>& cin;
        static std::basic_ostream<char, std::char_traits<char>>& cout;
        static std::basic_ostream<char, std::char_traits<char>>& cerr;
    };

    template<class Allocator>
    struct StdStreams<std::basic_string<wchar_t, std::char_traits<wchar_t>, Allocator>>
    {
        static std::basic_istream<wchar_t, std::char_traits<wchar_t>>& cin;
        static std::basic_ostream<wchar_t, std::char_traits<wchar_t>>& cout;
        static std::basic_ostream<wchar_t, std::char_traits<wchar_t>>& cerr;
    };

    template<class Allocator>
    std::basic_istream<char, std::char_traits<char>>& StdStreams<std::basic_string<char, std::char_traits<char>, Allocator>>::cin = std::cin;

    template<class Allocator>
    std::basic_ostream<char, std::char_traits<char>>& StdStreams<std::basic_string<char, std::char_traits<char>, Allocator>>::cout = std::cout;

    template<class Allocator>
    std::basic_ostream<char, std::char_traits<char>>& StdStreams<std::basic_string<char, std::char_traits<char>, Allocator>>::cerr= std::cerr;

    template<class Allocator>
    std::basic_istream<wchar_t, std::char_traits<wchar_t>>& StdStreams<std::basic_string<wchar_t, std::char_traits<wchar_t>, Allocator>>::cin = std::wcin;

    template<class Allocator>
    std::basic_ostream<wchar_t, std::char_traits<wchar_t>>& StdStreams<std::basic_string<wchar_t, std::char_traits<wchar_t>, Allocator>>::cout = std::wcout;

    template<class Allocator>
    std::basic_ostream<wchar_t, std::char_traits<wchar_t>>& StdStreams<std::basic_string<wchar_t, std::char_traits<wchar_t>, Allocator>>::cerr = std::wcerr;


    template<
        class Chr = char,
        class Traits = std::char_traits<Chr>,
        class Allocator = std::allocator<Chr>>
    void Print(std::basic_ostream<Chr, Traits>& os,
                const std::basic_string<Chr, Traits, Allocator>& str, Chr separator = Chr('\n'),
                const std::basic_string<Chr, Traits, Allocator>& preface = std::basic_string<Chr, Traits, Allocator>(), int width = 80)
    {
        int w = 0;
        for (auto chr = str.begin(); chr != str.end(); )
        {
            if (*chr == separator)
            {
                os << separator << preface;
                w = 0;
                ++chr;
            }
            else if (w >= width)
            {
                os << separator << preface;
                w = 0;
            }
            else
            {
                os << *chr;
                ++w;
                ++chr;
            }
        }
    }

    template<class String, class Chr>
    String Convert(const Chr* c_str)
    {
        const String s(c_str, c_str + std::char_traits<Chr>::length(c_str));
        return s;
    }

    template<class Ty, class String>
    bool ReadVal(Ty& val, const String& arg)
    {
        std::basic_istringstream<typename String::value_type, typename String::traits_type, typename String::allocator_type> iss;
        iss.str(arg);
        iss >> val;
        //TODO check whether the string was fully consumed!
        return !iss.fail();
    }

    template<class Chr, class Traits, class Allocator, class String2>
    bool ReadVal(std::basic_string<Chr, Traits, Allocator>& val, const String2& arg)
    {
        val.assign(arg.begin(), arg.end());
        return true;
    }

    template<class String = std::string>
    struct Argument
    {
        typedef typename String::value_type Chr;
        typedef std::basic_ostream<Chr, typename String::traits_type> OStream;

        Argument(const std::initializer_list<const Chr*>& args,
                const String& info,
                const String& meta = String())
            : options(args.begin(), args.end()), _info(info), _meta(meta)
        {
        }
        virtual ~Argument() {}
        //! writes the corresponding help
        virtual void WriteLong(OStream& os, int width = 80) const
        {
            os << Chr('\t');
            if (this->options.empty())
            {
                os << _meta;
            }
            else
            {
                for (auto option : this->options)
                    os << option << Chr(' ');
                os << Chr('\'') << _meta << Chr('\'');
            }
            os << Convert<String>(" default: ");
            PrintVal(os);

            String prefix; prefix.push_back(Chr('\t')); prefix.push_back(Chr('\t'));

            if (!this->_info.empty())
            {
                os << Chr('\n') << Chr('\t') << Chr('\t');
                Print(os, this->_info, Chr('\n'), prefix, width);
            }
            PrintChoices(os);

            os << std::endl;
        }
        virtual void WriteShort(OStream& os)const
        {
            if (this->options.empty())
                os << _meta;
            else
                os << this->options[0] << Chr(' ') << Chr('\'') << _meta << Chr('\'');
        }
        //! returns the number of arguments consumed. If 0 then the read was not successful.
        virtual int Read(int argc, const Chr** argv)const = 0;

        const std::vector<String> options;
        const String _info;
        const String _meta;

        bool Match(const Chr* arg) const
        {
            if (options.empty())
                return true;
            for (const String& option : options)
                if (option == arg)
                    return true;
            return false;
        }
        virtual void PrintChoices(OStream& os) const = 0;
        virtual void PrintVal(OStream& os) const = 0;
    };

    template<class Ty, class String>
    struct TypedArgument : Argument<String>
    {
        using typename Argument<String>::Chr;
        using typename Argument<String>::OStream;

        TypedArgument(Ty& def_val, 
            const std::initializer_list<const Chr*>& args = {},
            const String& info = String(), const String& meta = String(),
            const std::initializer_list<Ty>& choices = {})
            :   Argument<String>(args, info, meta.empty() ? Convert<String>(typeid(_val).name()) : meta),
                _val(def_val), _choices(choices.begin(), choices.end())
        {
        }
        virtual ~TypedArgument() {}

        //! returns the number of arguments consumed. If 0 then the read was not successful.
        virtual int Read(int argc, const Chr** argv)const
        {
            if (argc > 0)
            {
                if (this->options.empty())
                {
                    return ReadVal<Ty, String>(_val, argv[0]) ? 1 : 0;
                }else if (this->Match(argv[0]) && argc > 1)
                {
                    if (ReadVal<Ty, String>(_val, argv[1]) && (_choices.empty() || _choices.find(_val) != _choices.end()))
                        return 2;
                    StdStreams<String>::cerr << "At option \"" << argv[0] << "\" the argument \"" << argv[1] << "\" is not valid!" << std::endl;
                    exit(1);
                }
            }
            return 0;
        }

        virtual void PrintChoices(OStream& os)const
        {
            if (!_choices.empty())
            {
                os << Convert<String>("\n\t\tpossible values:");
                for (const auto& choice : _choices)
                    os << Chr(' ') << choice;
            }
        }
        virtual void PrintVal(OStream& os)const
        {
            os << _val;
        }
        Ty& _val;
        const std::unordered_set<Ty> _choices;
    };

    template<class Chr2, class String, class Traits2, class Allocator2>
    struct TypedArgument<std::basic_string<Chr2, Traits2, Allocator2>, String> : Argument<String>
    {
        using typename Argument<String>::Chr;
        using typename Argument<String>::OStream;
        typedef std::basic_string<Chr2, Traits2, Allocator2> Ty;

        TypedArgument(Ty& def_val,
            const std::initializer_list<const Chr*>& args = {},
            const String& info = String(), const String& meta = String(),
            const std::initializer_list<Ty>& choices = {})
            : Argument<String>(args, info, meta.empty() ? Convert<String>("string") : meta),
            _val(def_val), _choices(choices.begin(), choices.end())
        {
        }
        virtual ~TypedArgument() {}

        //! returns the number of arguments consumed. If 0 then the read was not successful.
        virtual int Read(int argc, const Chr** argv)const
        {
            if (argc > 0)
            {
                if (this->options.empty())
                {
                    _val.assign(argv[0], argv[0] + String::traits_type::length(argv[0]));
                    if (_choices.empty() || _choices.find(_val) != _choices.end())
                        return 1;
                    
                    StdStreams<String>::cerr << Convert<String>("At option \"");
                    this->WriteShort(StdStreams<String>::cerr);
                    StdStreams<String>::cerr << Convert<String>("\" the argument \"") << argv[0] << Convert<String>("\" is not valid!");
                    this->PrintChoices(StdStreams<String>::cerr);
                    exit(1);
                }
                else if (this->Match(argv[0]) && argc > 1)
                {
                    _val.assign(argv[1], argv[1] + String::traits_type::length(argv[1]));
                    if (_choices.empty() || _choices.find(_val) != _choices.end())
                        return 2;

                    StdStreams<String>::cerr << Convert<String>("At option \"");
                    this->WriteShort(StdStreams<String>::cerr);
                    StdStreams<String>::cerr << Convert<String>("\" the argument \"") << argv[1] << Convert<String>("\" is not valid!");
                    this->PrintChoices(StdStreams<String>::cerr);
                    exit(1);
                }
            }
            return 0;
        }

        virtual void PrintChoices(OStream& os)const
        {
            if (!_choices.empty())
            {
                auto buffer = Convert<String>("\n\t\tpossible values:");
                os << buffer;
                for (const auto& choice : _choices)
                {
                    buffer.assign(choice.begin(), choice.end());
                    os << Chr(' ') << Chr('"') << buffer << Chr('"');
                }
            }
        }
        virtual void PrintVal(OStream& os)const
        {
            const String buffer(_val.begin(), _val.end());
            os << Chr('"') << buffer << Chr('"');
        }
        Ty& _val;
        const std::unordered_set<Ty> _choices;
    };

    template<class String>
    struct SetFlag : Argument<String>
    {
        using typename Argument<String>::Chr;
        using typename Argument<String>::OStream;

        SetFlag(bool& def_val,
            const std::initializer_list<const Chr*>& args = {},
            const String& info = String(), bool reset=false, const String& meta = String())
            : Argument<String>(args, info, meta.empty() ? Convert<String>(typeid(_val).name()) : meta),
                _val(def_val), _reset(reset)
        {
        }
        virtual ~SetFlag() {}
        virtual void PrintVal(OStream& os)const
        {
            os << (_val ? "true" : "false");
        }
        //! returns the number of arguments consumed. If 0 then the read was not successful.
        virtual int Read(int argc, const Chr** argv)const
        {
            if (argc > 0)
            {
                if (this->Match(argv[0]))
                {
                    _val = _reset ? false : true;
                    return 1;
                }
            }
            return 0;
        }
        virtual void PrintChoices(OStream& )const
        {
        }
        virtual void WriteShort(OStream& os)const
        {
            os << this->options[0];
        }
        bool& _val;
        const bool _reset;
    };

    //! TODO: line break, optional, default and intentional, too
    //! TODO range checking functions!
    template<class String = std::string>
    class Parser
    {
        typedef typename String::value_type Chr;
        typedef std::basic_ostream<Chr, typename String::traits_type> OStream;
    public:
        Parser(const String& header,
                const std::initializer_list<const Chr*>& helps,
                const String& footer = String(),
                int width = 80, bool strict=false)
            :   positional_arguments(), optional_arguments(),
                help_options(helps.begin(), helps.end()), _options(),
                program_name(), header(header), footer(footer),
                width(width), strict(strict)
        {
        }
        ~Parser() {}

        void Do(int argc, const Chr** argv)
        {
            auto positional = positional_arguments.begin();
            program_name = argv[0];
            for (++argv, --argc; argc > 0;)
            {
                for (const String& help : help_options)
                {
                    if (help == *argv)
                    {
                        HelpLong(StdStreams<String>::cout);
                        exit(0);
                    }
                }
                bool found = false;
                for (const auto& argument : optional_arguments)
                {
                    const auto read = argument->Read(argc, argv);
                    if (read > 0)
                    {
                        argc -= read;
                        argv += read;
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    if (positional != positional_arguments.end())
                    {
                        const auto read = (*positional)->Read(argc, argv);
                        argc -= read;
                        argv += read;
                        ++positional;
                    }
                    else
                    {
                        StdStreams<String>::cerr << Convert<String>("Unknown argument: \"") << *argv << Chr('\"') << std::endl;
                        if (strict)
                        {
                            StdStreams<String>::cerr << std::endl;
                            HelpShort(StdStreams<String>::cerr);
                            exit(1);
                        }
                        ++argv;
                        --argc;
                    }
                }
            }
            if (positional != positional_arguments.end())
            {
                StdStreams<String>::cerr << Convert<String>("Positional argument \"");
                (*positional)->WriteShort(StdStreams<String>::cerr);
                StdStreams<String>::cerr << Convert<String>("\" was not provided!\n") << std::endl;
                HelpShort(StdStreams<String>::cerr);
                exit(1);
            }
        }
        void HelpShort(OStream& os = std::cout) const
        {
            os << Convert<String>("HELP:\n") << program_name;
            for (const auto& h : help_options)
                os << Chr(' ') << h;

            os << Convert<String>("\n\nUSAGE:\n") << program_name;
            for (const auto& argument : positional_arguments)
            {
                os << Chr(' ');
                argument->WriteShort(os);
            }
            for (const auto& argument : optional_arguments)
            {
                os << Chr(' ') << Chr('[');
                argument->WriteShort(os);
                os << Chr(']');
            }
            os << std::endl;
        }
        void HelpLong(OStream& os = std::cout) const
        {
            if (!header.empty())
            {
                Print(os, header, Chr('\n'), String(), width);
                os << Chr('\n');
            }
            os << Convert<String>("\nPOSITIONALS:\n");
            for (const auto& argument : positional_arguments)
                argument->WriteLong(os);

            os << Convert<String>("\nOPTIONS:\n");
            for (const auto& argument : optional_arguments)
                argument->WriteLong(os);

            if (!footer.empty())
            {
                os << Chr('\n');
                Print(os, footer, Chr('\n'), String(), width);
                os << Chr('\n');
            }
        }
        template <class Ty>
        void AddArg(Ty& value,
            const std::initializer_list<const Chr*>& args = {},
            const String& info = String(), const String& meta = String(),
            const std::initializer_list<Ty>& choices = {})
        {
            if (args.size() == 0)
            {
                positional_arguments.emplace_back(new TypedArgument<Ty, String>(value, args, info, meta, choices));
            }
            else
            {
                optional_arguments.emplace_back(new TypedArgument<Ty, String>(value, args, info, meta, choices));
                if (choices.size() > 0 &&
                    std::find(choices.begin(), choices.end(), value) == choices.end())
                {   // has restrictions but the initial value does not match.
                    // the initial value is not within choices
                    StdStreams<String>::cerr << Convert<String>("Default value ");
                    optional_arguments.back()->PrintVal(StdStreams<String>::cerr);
                    StdStreams<String>::cerr << Convert<String>(" is not allowed for '") << *args.begin() << Chr('\'');
                    optional_arguments.back()->PrintChoices(StdStreams<String>::cerr);
                    StdStreams<String>::cerr << std::endl;
                    exit(1);
                }
                CheckLastArgumentOptions();
            }
        }
        void AddFlag(bool& value,
            const std::initializer_list<const Chr*>& args = {},
            const String& info = String(), bool reset=false, const String& meta = String())
        {
            if (args.size() == 0)
            {
                StdStreams<String>::cerr << Convert<String>("Boolean (flag) argument cannot be positional!") << std::endl;
                exit(1);
            }
            optional_arguments.emplace_back(new SetFlag<String>(value, args, info, reset, meta));
            CheckLastArgumentOptions();
        }

    private:
        void CheckLastArgumentOptions()
        {
            for (const auto& option : optional_arguments.back()->options)
            {
                if (!_options.insert(option).second)
                {
                    StdStreams<String>::cerr << Convert<String>("Duplicate option: \"") << option << Chr('\"') << std::endl;
                    exit(1);
                }
            }
        }
        std::vector<std::unique_ptr<Argument<String>>> positional_arguments;
        std::vector<std::unique_ptr<Argument<String>>> optional_arguments;
        std::vector<String> help_options;
        std::unordered_set<String> _options;
        String program_name;
        String header, footer;
        int width;
        bool strict;
    };
}
