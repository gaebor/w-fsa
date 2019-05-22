#pragma once

#include <typeinfo>
#include <string>
#include <vector>
#include <memory>
#include <unordered_set>
#include <iostream>
#include <sstream>

namespace arg
{
    template<
        class Chr = char,
        class Traits = std::char_traits<Chr>,
        class Allocator = std::allocator<Chr>>
    void Print(std::basic_ostream<Chr, Traits>& os,
                const std::basic_string<Chr, Traits, Allocator>& str, Chr separator = Chr(10),
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

    template<
        class Chr,
        class Ty,
        class Traits = std::char_traits<Chr>,
        class Allocator = std::allocator<Chr>>
    bool ReadVal(Ty& val, const std::basic_string<Chr, Traits, Allocator>& arg)
    {
        std::basic_istringstream<Chr, Traits, Allocator> iss;
        iss.str(arg);
        iss >> val;
        return !iss.fail();
    }

    template<
        class Chr = char,
        class Traits = std::char_traits<Chr>,
        class Allocator = std::allocator<Chr>>
    struct Argument
    {
        typedef std::basic_string<Chr, Traits, Allocator> String;
        typedef std::basic_ostream<Chr, Traits> OStream;
        Argument(const std::initializer_list<const Chr*>& args, const String& info)
            : options(args.begin(), args.end()), _info(info)
        {
        }
        virtual ~Argument() {}
        //! writes the corresponding help
        virtual void WriteLong(OStream& os, int width = 80)const = 0;
        virtual void WriteShort(OStream& os)const = 0;
        //! returns the number of arguments consumed. If 0 then the read was not successful.
        virtual int Read(int argc, const Chr** argv)const = 0;

        const std::vector<String> options;
        const String _info;
    protected:
        bool Match(const Chr* arg) const
        {
            if (options.empty())
                return true;
            for (const auto& option : options)
                if (option == arg)
                    return true;
            return false;
        }
    };

    template<class Ty,
        class Chr = char,
        class Traits = std::char_traits<Chr>,
        class Allocator = std::allocator<Chr>>
    struct TypedArgument : Argument<Chr, Traits, Allocator>
    {
        using typename Argument<Chr, Traits, Allocator>::String;
        using typename Argument<Chr, Traits, Allocator>::OStream;

        TypedArgument(Ty& def_val, 
            const std::initializer_list<const Chr*>& args = {},
            const String& info = String(), const std::string& meta = "",
            const std::initializer_list<Ty>& choices = {})
            :   Argument<Chr, Traits, Allocator>(args, info), _val(def_val),
                _meta(meta.empty() ? typeid(_val).name() : meta),
                _choices(choices.begin(), choices.end())
        {
        }
        virtual ~TypedArgument() {}
        virtual void WriteShort(OStream& os)const
        {
            os << " [" << this->options[0] << " '" << String(_meta.begin(), _meta.end()) << "']";
        }
        //! writes the corresponding help
        virtual void WriteLong(OStream& os, int width = 80)const
        {
            os << "\t";
            for (auto option : this->options)
                os << option << " ";
            os << "'" << String(_meta.begin(), _meta.end()) << "' default: " << _val;
            String prefix; prefix.push_back(Chr(9)); prefix.push_back(Chr(9));
            if (!this->_info.empty())
            {
                os << "\n\t\t";
                Print(os, this->_info, Chr(10), prefix, width);
            }
            
            if (!_choices.empty())
            {
                os << "\n\t\tpossible values:";
                for (const auto& choice : _choices)
                    os << ' ' << choice;
            }
            os << std::endl;
        }
        //! returns the number of arguments consumed. If 0 then the read was not successful.
        virtual int Read(int argc, const Chr** argv)const
        {
            if (argc > 0)
            {
                if (this->options.empty())
                {
                    return ReadVal<Chr, Ty, Traits, Allocator>(_val, argv[0]) ? 1 : 0;
                }else if (this->Match(argv[0]) && argc > 1)
                {
                    if (ReadVal<Chr, Ty, Traits, Allocator>(_val, argv[1]) && (_choices.empty() || _choices.find(_val) != _choices.end()))
                        return 2;
                    std::cerr << "At option \"" << argv[0] << "\" the argument \"" << argv[1] << "\" is not valid!" << std::endl;
                    exit(1);
                }
            }
            return 0;
        }
    protected:
        Ty & _val;
        const std::string _meta;
        std::unordered_set<Ty> _choices;
    };

    template<class Chr, class Traits, class Allocator>
    struct TypedArgument<std::basic_string<Chr, Traits, Allocator>> : Argument<Chr, Traits, Allocator>
    {
    private:
        typedef std::basic_string<Chr, Traits, Allocator> Ty;
        using typename Argument<Chr, Traits, Allocator>::String;
        using typename Argument<Chr, Traits, Allocator>::OStream;

    public:
        TypedArgument(Ty& def_val,
            const std::initializer_list<const char*>& args = {},
            const std::string& info = "", const std::string& meta = "",
            const std::initializer_list<Ty>& choices = {})
            : Argument<Chr, Traits, Allocator>(args, info), _val(def_val),
            _meta(meta.empty() ? typeid(_val).name() : meta),
            _choices(choices.begin(), choices.end())
        {
        }
        virtual ~TypedArgument() {}
        virtual void WriteShort(std::ostream& os)const
        {
            os << " [" << this->options[0] << " '" << _meta << "']";
        }
        //! writes the corresponding help
        virtual void WriteLong(std::ostream& os, int width = 80)const
        {
            os << "\t";
            for (auto option : this->options)
                os << option << " ";
            os << "'" << String(_meta.begin(), _meta.end()) << "' default: \"" << _val << "\"";

            String prefix; prefix.push_back(Chr(9)); prefix.push_back(Chr(9));

            if (!this->_info.empty())
            {
                os << "\n\t\t";
                Print(os, this->_info, Chr(10), prefix, width);
            }

            if (!_choices.empty())
            {
                os << "\n\t\tpossible values:";
                for (const auto& choice : _choices)
                    os << " \"" << choice << "\"";
            }
            os << std::endl;
        }
        //! returns the number of arguments consumed. If 0 then the read was not successful.
        virtual int Read(int argc, const Chr** argv)const
        {
            if (argc > 0)
            {
                if (this->options.empty())
                {
                    return ReadVal<Chr, String, Traits, Allocator>(_val, argv[0]) ? 1 : 0;
                }
                else if (this->Match(argv[0]) && argc > 1)
                {
                    if (ReadVal<Chr, String, Traits, Allocator>(_val, argv[1]) && (_choices.empty() || _choices.find(_val) != _choices.end()))
                        return 2;
                    std::wcerr << "At option \"" << argv[0] << "\" the argument \"" << argv[1] << "\" is not valid!" << std::endl;
                    exit(1);
                }
            }
            return 0;
        }
    protected:
        Ty & _val;
        const std::string _meta;
        std::unordered_set<Ty> _choices;
    };

    template<class Chr = char,
        class Traits = std::char_traits<Chr>,
        class Allocator = std::allocator<Chr>>
    struct SetFlag : Argument<Chr, Traits, Allocator>
    {
        using typename Argument<Chr, Traits, Allocator>::String;
        using typename Argument<Chr, Traits, Allocator>::OStream;

        SetFlag(bool& def_val,
            const std::initializer_list<const Chr*>& args = {},
            const String& info = String(), bool reset=false, const std::string& meta = "")
            : Argument<Chr, Traits, Allocator>(args, info), _val(def_val), _reset(reset),
            _meta(meta.empty() ? typeid(_val).name() : meta)
        {
        }
        virtual ~SetFlag() {}
        virtual void WriteShort(OStream& os)const
        {
            os << " [" << this->options[0] << "]";
        }
        //! writes the corresponding help
        virtual void WriteLong(OStream& os, int width = 80)const
        {
            os << "\t";
            for (auto option : this->options)
                os << option << " ";
            os << "'" << String(_meta.begin(), _meta.end()) << "' default: " << (_val ? "true" : "false");

            String prefix; prefix.push_back(Chr(9)); prefix.push_back(Chr(9));

            if (!this->_info.empty())
            {
                os << "\n\t\t";
                Print(os, this->_info, Chr(10), prefix, width);
            }

            os << std::endl;
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
    protected:
        bool& _val;
        bool _reset;
        const std::string _meta;
    };

    //! TODO: line break, optional, default and intentional, too
    //! TODO range checking functions!
    template<class Chr = char,
        class Traits = std::char_traits<Chr>,
        class Allocator = std::allocator<Chr>>
    class Parser
    {
        typedef std::basic_string<Chr, Traits, Allocator> String;
        typedef std::basic_ostream<Chr, Traits> OStream;
    public:
        Parser(const String& info,
                const std::initializer_list<const Chr*>& helps,
                int width = 80)
            :   arguments(), help_options(helps.begin(), helps.end()), _options(),
                program_name(), header(info),
                width(width)
        {
        }
        ~Parser() {}

        void Do(int argc, const Chr** argv)
        {
            program_name = argv[0];
            for (++argv, --argc; argc > 0;)
            {
                for (const auto& help : help_options)
                {
                    if (help == *argv)
                    {
                        Help(std::cout);
                        exit(0);
                    }
                }
                bool found = false;
                for (const auto& argument : arguments)
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
                    std::cerr << "Unknown argument: \"" << *argv << "\"" << std::endl;
                    ++argv;
                    --argc;
                }
            }
        }

        void Help(OStream& os = std::cout) const
        {
            if (!header.empty())
            {
                Print(os, header, Chr(10), String(), width);
                os << "\n\n";
            }
            os << "USAGE:\n\t" << program_name;
            for (const auto& argument : arguments)
                argument->WriteShort(os);
            os << "\n";
            os << "OPTIONS:\n";
            for (const auto& argument : arguments)
                argument->WriteLong(os);
        }
        template <class Ty>
        void AddArg(Ty& value,
            const std::initializer_list<const Chr*>& args = {},
            const String& info = String(), const std::string& meta = "",
            const std::initializer_list<Ty>& choices = {})
        {
            arguments.emplace_back(new TypedArgument<Ty, Chr, Traits, Allocator>(value, args, info, meta, choices));

            for (const auto& option : arguments.back()->options)
                if (!_options.insert(option).second)
                {
                    std::cerr << "Duplicate option: \"" << option  << "\"" << std::endl;
                    exit(1);
                }
        }
        void AddFlag(bool& value,
            const std::initializer_list<const Chr*>& args = {},
            const String& info = String(), bool reset=false, const std::string& meta = "")
        {
            arguments.emplace_back(new SetFlag<Chr, Traits, Allocator>(value, args, info, reset, meta));

            for (const auto& option : arguments.back()->options)
                if (!_options.insert(option).second)
                {
                    std::cerr << "Duplicate option: \"" << option << "\"" << std::endl;
                    exit(1);
                }
        }

    private:
        std::vector<std::unique_ptr<Argument<Chr, Traits, Allocator>>> arguments;
        std::vector<String> help_options;
        std::unordered_set<String> _options;
        String program_name;
        String header;
        int width;
    };
}
