#pragma once

#include <typeinfo>
#include <type_traits>
#include <string>
#include <cstring>
#include <vector>
#include <memory>
#include <unordered_set>
#include <set>
#include <sstream>
#include <ostream>
#include <algorithm>
#include <functional>
#include <utility>

namespace arg
{

    template<class String, class OStream>
    void PrintWidth(OStream& os,
                const String& str, typename String::value_type separator /*= typename String::value_type('\n')*/,
                const String& preface /*= String()*/, int width /*= 80*/)
    {
        typedef typename String::value_type Chr;
        int w = 0;
        width = std::max(1, width - (int)preface.size());
        for (auto chr = str.begin(); chr != str.end(); )
        {
            if (*chr == separator)
            {
                os << separator << preface;
                w = 0;
                ++chr;
            }
            else if (w >= width && *chr == Chr(' '))
            {
                os << separator << preface;
                w = 0;
                ++chr;
            }
            else
            {
                os << *chr;
                ++w;
                ++chr;
            }
        }
    }

    template<class Ty, class String>
    String GetTypeName(const Ty& val, const String& meta)
    {
        if (meta.empty())
        {
            std::basic_ostringstream<typename String::value_type> buffer;
            buffer << typeid(val).name();
            return buffer.str();
        }
        else
            return meta;
    }

    template<class Ty, class Chr>
    bool ReadVal(Ty& val, const Chr* arg)
    {
        typedef std::basic_string<Chr> String;
        std::basic_istringstream<typename String::value_type, typename String::traits_type, typename String::allocator_type> iss;
        iss.str(arg);
        iss >> val;
        return !iss.fail() && iss.eof();
    }
    
    template<class Chr, class Traits, class Allocator, class Chr2>
    bool ReadVal(std::basic_string<Chr, Traits, Allocator>& val, const Chr2* arg)
    {
        val.assign(arg, arg + std::char_traits<Chr2>::length(arg));
        return true;
    }
    
    template<class Chr>
    bool ReadVal(bool& val, const Chr* arg)
    {
        typedef std::basic_string<Chr> String;
        std::basic_istringstream<typename String::value_type, typename String::traits_type, typename String::allocator_type> iss;
        iss.str(arg);
        iss >> val;
        if (iss.fail() || !iss.eof())
        {
            iss.clear();
            iss.str(arg);
            iss >> std::boolalpha >> val;
            return !iss.fail();
        }else
            return true;
    }
    
    template<class Ty, class OStream>
    void PrintVal(const Ty& val, OStream& os)
    {
        os << val;
    }

    template<class Chr, class Traits, class Allocator, class OStream>
    void PrintVal(const std::basic_string<Chr, Traits, Allocator>& val, OStream& os)
    {
        os << '"' << val << '"';
    }
    
    template<class OStream>
    void PrintVal(const bool& val, OStream& os)
    {
        const auto flag_state = os.flags();
        os << std::boolalpha << val;
        os.flags(flag_state);
    }

    template<class Ty, class OStream>
    struct Checker
    {
        Checker() {}
        virtual ~Checker(){};
        virtual bool Do(const Ty&) const
        {
            return true;
        }
        
        virtual void Print(OStream&, int)const {}
        
        static const Checker<Ty, OStream>* New() {return new Checker<Ty, OStream>();}
        static const Checker<Ty, OStream>* New(Ty a);
        static const Checker<Ty, OStream>* New(Ty a, Ty b);
        template<class ArgList>
        static const Checker<Ty, OStream>* New(const ArgList& list);
        static const Checker<Ty, OStream>* New(std::function<bool(const Ty&)> cond);
    };
    
    template<class Ty, class OStream>
    struct MinMaxChecker : Checker<Ty, OStream>
    {
        MinMaxChecker(const Ty& a, const Ty& b)
            : _min(a), _max(b)
        {}
        virtual bool Do(const Ty& val) const
        {
            return _min <= val && (_min >= _max || val <= _max);
        }
        virtual void Print(OStream& os, int width)const
        {
            typedef typename OStream::char_type Chr;
            typedef std::basic_string<Chr> String;
            std::basic_ostringstream<Chr, typename OStream::traits_type> buffer;
            
            buffer << "\t\t";
            const String preface = buffer.str();
            buffer.str(String());

            buffer << "\nmin: ";
            PrintVal(_min, buffer);
            if (_max > _min)
            {   
                buffer << ", max: ";
                PrintVal(_max, buffer);
            }
            PrintWidth(os, buffer.str(), Chr('\n'), preface, width);
        }
        const Ty _min;
        const Ty _max;
    };
    template<class Ty, class OStream>
    const Checker<Ty, OStream>* Checker<Ty, OStream>::New(Ty a, Ty b)
    {
        return new MinMaxChecker<Ty, OStream>(a, b);
    }
    template<class Ty, class OStream>
    const Checker<Ty, OStream>* Checker<Ty, OStream>::New(Ty a)
    {
        return new MinMaxChecker<Ty, OStream>(a, a);
    }
    
    template<class Ty, class OStream>
    struct ChoiceChecker : Checker<Ty, OStream>
    {
        template<class ArgList>
        ChoiceChecker(const ArgList& list)
            : choices(list.begin(), list.end())
        {}
        virtual bool Do(const Ty& val) const
        {
            return choices.empty() || (choices.find(val) != choices.end());
        }
        virtual void Print(OStream& os, int width)const
        {
            typedef typename OStream::char_type Chr;
            typedef std::basic_string<Chr> String;

            if (!choices.empty())
            {
                std::basic_ostringstream<Chr, typename OStream::traits_type> buffer;
                buffer << "\t\t";
                const String preface = buffer.str();
                buffer.str(String());

                buffer << "\npossible values:";
                for (const auto& choice : choices)
                {
                    buffer << ' ';
                    PrintVal(choice, buffer);
                }
                PrintWidth(os, buffer.str(), Chr('\n'), preface, width);
            }
        }
        const std::set<Ty> choices;
    };
    template<class Ty, class OStream>
    template<class ArgList>
    const Checker<Ty, OStream>* Checker<Ty, OStream>::New(const ArgList& list)
    {
        return new ChoiceChecker<Ty, OStream>(list);
    }
    
    template<class Ty, class OStream>
    struct ConditionChecker : Checker<Ty, OStream>
    {
        ConditionChecker(std::function<bool(const Ty&)> cond)
            : condition(cond)
        {}
        virtual bool Do(const Ty& val) const
        {
            return condition(val);
        }
        const std::function<bool(const Ty&)> condition;
    };
    template<class Ty, class OStream>
    const Checker<Ty, OStream>* Checker<Ty, OStream>::New(std::function<bool(const Ty&)> cond)
    {
        return new ConditionChecker<Ty, OStream>(cond);
    }
    
    template<class String, class OStream>
    struct Argument
    {
        typedef typename String::value_type Chr;

        template<class ArgList>
        Argument(const ArgList& args, const String& info, const String& meta)
            : options(args.begin(), args.end()), _info(info), _meta(meta)
        {
        }
        virtual ~Argument() {}
        //! writes the corresponding help
        virtual void WriteLong(OStream& os, int width = 80) const
        {
            os << '\t';
            if (this->options.empty())
            {
                WriteShort(os);
            }
            else
            {
                for (auto option : this->options)
                    os << option << ' ';
                os << '\'' << _meta << '\'';
            }
            os << " default: "; Print(os);
            
            std::basic_ostringstream<Chr> buffer;
            buffer << "\t\t";
            const String prefix = buffer.str();
            
            if (!this->_info.empty())
            {
                os << "\n\t\t";
                PrintWidth(os, this->_info, Chr('\n'), prefix, width);
            }
            PrintRestrictions(os, width);

            os << std::endl;
        }
        virtual void WriteShort(OStream& os)const
        {
            if (!this->options.empty())
                os << this->options[0] << ' ';
            
            os << '\'' << _meta << '\'';
        }
        //! returns the number of arguments consumed. If 0 then the read was not successful.
        /*! writes cerr only if something went wrong.
        */
        virtual int Read(int argc, const Chr** argv, OStream& cerr)const = 0;
        virtual void Print(OStream& os)const = 0;
        virtual void PrintRestrictions(OStream& os, int width) const = 0;

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
    };

    template<class Ty, class String, class OStream>
    struct TypedArgument : Argument<String, OStream>
    {
        using typename Argument<String, OStream>::Chr;

        template<class ArgList, class... _Valty>
        TypedArgument(Ty& def_val,
            const ArgList& args, const String& info, const String& meta,
            _Valty... _Val)
            : Argument<String, OStream>(args, info, GetTypeName(def_val, meta)),
            _val(def_val), checker(Checker<Ty, OStream>::New(_Val...))
        {
        }
        virtual ~TypedArgument() {}

        virtual int Read(int argc, const Chr** argv, OStream& cerr)const
        {
            int ret = 0;
            if (!this->options.empty())
            {
                if (this->Match(argv[0]))
                {
                    ++argv;
                    --argc;
                    ++ret;
                }else
                    return 0;
            }
            
            if (argc > 0)
            {
                Ty val;
                if (ReadVal(val, argv[0]) && checker->Do(val))
                {
                    ++ret;
                    std::swap(_val, val);
                }
                else
                {
                    cerr << "At option ";
                    this->WriteShort(cerr);
                    cerr << " the argument \"" << argv[0] << "\" is not valid!";
                    this->PrintRestrictions(cerr, 80);
                    cerr << std::endl;
                    return -1;
                }
            } else
            {
                cerr << "Option ";
                this->WriteShort(cerr);
                cerr << " requires an argument!" << std::endl;
                return -1;
            }
            return ret;
        }
        
        virtual void Print(OStream& os)const
        {
            PrintVal(_val, os);
        }

        virtual void PrintRestrictions(OStream& os, int width)const
        {
            checker->Print(os, width);
        }
        Ty& _val;
        std::unique_ptr<const Checker<Ty, OStream>> checker;
    };

    template<class String, class OStream>
    struct SetFlag : Argument<String, OStream>
    {
        using typename Argument<String, OStream>::Chr;

        template<class ArgList>
        SetFlag(bool& def_val,
            const ArgList& args,
            const String& info, bool reset, const String& meta)
            : Argument<String, OStream>(args, info, GetTypeName(def_val, meta)),
                _val(def_val), _reset(reset)
        {
        }
        virtual ~SetFlag() {}
        //! returns the number of arguments consumed. If 0 then the read was not successful.
        virtual int Read(int argc, const Chr** argv, OStream& )const
        {
            if (argc > 0 && this->Match(argv[0]))
            {
                _val = _reset ? false : true;
                return 1;
            }
            return 0;
        }
        virtual void Print(OStream& os)const
        {
            PrintVal(_val, os);
        }
        virtual void PrintRestrictions(OStream&, int width)const {}
        virtual void WriteShort(OStream& os)const
        {
            os << this->options[0];
        }
        bool& _val;
        const bool _reset;
    };

    typedef std::initializer_list<const char*> list;
    typedef std::initializer_list<const wchar_t*> wlist;

    template<class String = std::string, bool strict = true, class OStream = std::ostream>
    class Parser
    {
        typedef typename String::value_type Chr;
        static_assert(std::is_base_of<std::basic_ostream<Chr, typename String::traits_type>, OStream>::value, "OStream type should be child of std::basic_ostream");
        bool fail()const
        {
            if (strict)
                exit(1);
            else
                return false;
            return true;
        }
        bool succeed()const
        {
            if (strict)
                exit(0);
            return true;
        }

    public:
        template<class ArgList = std::initializer_list<String>>
        Parser(const String& header,
                const ArgList& helps,
                OStream& cout, OStream& cerr,
                const String& footer = String(),
                int width = 80)
            :   positional_arguments(), optional_arguments(),
                cout(cout), cerr(cerr),
                help_options(helps.begin(), helps.end()), _options(),
                program_name(), header(header), footer(footer),
                width(width)
        {
        }
        ~Parser() {}

        bool Do(int argc, const Chr** argv)
        {
            auto positional = positional_arguments.begin();
            program_name = argv[0];
            for (++argv, --argc; argc > 0;)
            {
                for (const String& help : help_options)
                {
                    if (help == *argv)
                    {
                        HelpLong(cout);
                        return succeed();
                    }
                }
                bool found = false;
                for (const auto& argument : optional_arguments)
                {
                    const auto read = argument->Read(argc, argv, cerr);
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
                        const auto read = (*positional)->Read(argc, argv, cerr);
                        if (read < 0)
                        {
                            HelpShort(cerr << std::endl);
                            return fail();
                        }
                        argc -= read;
                        argv += read;
                        ++positional;
                    }
                    else
                    {
                        cerr << "Unknown argument: \"" << *argv << '\"' << std::endl;
                        HelpShort(cerr << std::endl);
                        fail(); // only if 'strict'
                        ++argv;
                        --argc;
                    }
                }
            }
            if (positional != positional_arguments.end())
            {
                cerr << "Positional argument ";
                (*positional)->WriteShort(cerr);
                cerr << " was not provided!\n" << std::endl;
                HelpShort(cerr);
                return fail();
            }
            return true;
        }
        void HelpShort(OStream& os) const
        {
            os << "HELP:\n" << program_name;
            for (const auto& h : help_options)
                os << ' ' << h;

            os << "\n\nUSAGE:\n" << program_name;
            for (const auto& argument : positional_arguments)
            {
                os << ' ';
                argument->WriteShort(os);
            }
            for (const auto& argument : optional_arguments)
            {
                os << " [";
                argument->WriteShort(os);
                os << ']';
            }
            os << std::endl;
        }
        void HelpLong(OStream& os) const
        {
            if (!header.empty())
            {
                PrintWidth(os, header, Chr('\n'), String(), width);
                os << Chr('\n');
            }
            os << "\nPOSITIONALS:\n";
            for (const auto& argument : positional_arguments)
                argument->WriteLong(os, width);

            os << "\nOPTIONALS:\n";
            for (const auto& argument : optional_arguments)
                argument->WriteLong(os, width);

            if (!footer.empty())
            {
                os << '\n';
                PrintWidth(os, footer, Chr('\n'), String(), width);
                os << '\n';
            }
        }
        
        template <class Ty, class ArgList = std::initializer_list<String>, class ...Vargs>
        bool AddArg(Ty& value,
            const ArgList& args = ArgList(),
            const String& info = String(), const String& meta = String(),
            Vargs ...vargs)
        {
            typedef TypedArgument<Ty, String, OStream> ArgType;
            std::unique_ptr<ArgType> arg(new ArgType(value, args, info, meta, vargs...));
            if (args.size() == 0)
            {
                positional_arguments.emplace_back(std::move(arg));
                return true;
            }
            else
            {
                if (!arg->checker->Do(value))
                {
                    cerr << "Default value ";
                    arg->Print(cerr);
                    cerr << " is not allowed for '" << *args.begin() << '\'';
                    arg->PrintRestrictions(cerr, width);
                    cerr << std::endl;
                    return fail();
                }
                return AddCheckedArgument(arg);
            }
        }

        template<class ArgList = std::initializer_list<String>>
        bool AddFlag(bool& value,
            const ArgList& args = ArgList(),
            const String& info = String(), bool reset=false, const String& meta = String())
        {
            if (args.size() == 0)
            {
                cerr << "Boolean (flag) argument cannot be positional!" << std::endl;
                return fail();
            }
            std::unique_ptr<Argument<String, OStream>> arg(new SetFlag<String, OStream>(value, args, info, reset, meta));
            return AddCheckedArgument(arg);
        }

    private:
        template<class Ptr>
        bool AddCheckedArgument(Ptr& arg)
        {
            for (const auto& option : arg->options)
            {
                if (_options.find(option) != _options.end())
                {
                    cerr << "Duplicate option: \"" << option << '\"' << std::endl;
                    return fail();
                }
            }
            _options.insert(arg->options.begin(), arg->options.end());
            optional_arguments.emplace_back(std::move(arg));
            return true;
        }
        std::vector<std::unique_ptr<Argument<String, OStream>>> positional_arguments;
        std::vector<std::unique_ptr<Argument<String, OStream>>> optional_arguments;
        OStream &cout, &cerr;
        std::vector<String> help_options;
        std::unordered_set<String> _options;
        String program_name;
        String header, footer;
        int width;
    };
}
