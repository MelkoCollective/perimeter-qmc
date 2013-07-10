// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    21.04.2012 13:19:09 CEST
// File:    bash_parameter_msk.hpp

#ifndef __BASH_PARAMETER_HEADER_MSK
#define __BASH_PARAMETER_HEADER_MSK

/* minimal code

#include <bash_parameter_msk.hpp>
#include "bash_parameter_msk.hpp"
//parameter is defined in the header
addon::parameter.set(energy, "energy");
addon::parameter.set(magn, "magn", 100); //set with default
addon::parameter.set("equilib"); //set iternal
addon::parameter.read(argc, argv);
addon::parameter.get();

*/

#include "color.hpp"
#include "serialize/archive_enum.hpp"

#include <map>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <stdexcept>
#include <memory>
#include <algorithm>

//timer2_msk.hpp documents addon
namespace addon
{
    namespace detail {
        //------------------- something similar to boost any -------------------
        class any_base {
        public:
            virtual void * get() = 0;
        };
        template<typename T>
        class any_der: public any_base {
        public:
            any_der(T const & t): t_(t) {}
            void * get() {
                return (void *)(&t_);
            }
            
        private:
            T t_;
        };
        //------------------- runtime checks -------------------
        template<typename T>
        bool convertable_to(const std::string& str) {
            std::istringstream iss(str);
            T obj;
            iss >> std::ws >> obj >> std::ws;

            if(!iss.eof())
                return false;

            return true; 
        }
        template< typename T >
        inline T convert(const std::string& str) {
            std::istringstream iss(str);
            T obj;

            iss >> std::ws >> obj >> std::ws;

            if(!iss.eof())
                throw std::runtime_error("not convertable!");

            return obj; 
        }

        //------------------- compiletime checks -------------------
        template<typename T, typename U>
        struct is_convertable {
            static char check(U const &);
            static double check(...);
            static T t;
            enum{value = (sizeof(char) == sizeof(check(t)))};
        };
        
        //------------------- chameleon type -------------------
        class chameleon_type {
            typedef std::string type_type;
            
            //------------------- init -------------------
            template<typename T>
            static type_type tn() {
                return type_type(typeid(T).name());
            }
            template<typename T, typename U, bool convertable>
            struct init {
                init<T, U, convertable>(T const & t, chameleon_type * b) {
                    b->type = tn<U>();
                    b->any = std::shared_ptr<any_base>(new any_der<U>(t));
                }
            };
            template<typename T, typename U>
            struct init<T, U, false> {
                init<T, U, false>(T const & t, chameleon_type * b) {
                }
            };
            template<typename T>
            void set(T const & t) {
                any.reset();
                #define init_op(U) init<T, U, is_convertable<T, U>::value>(t, this);
                
                init_op(double)
                init_op(std::string)
                
                #undef init_op
                if(!any) {
                    std::stringstream ss;
                    ss << "cannot convert type " << typeid(T).name() << " to double or std::string";
                    throw std::runtime_error(ss.str());
                }
            }
        public:
            //------------------- ctors -------------------
            template<typename T>
            chameleon_type(T const & t) {
                set(t);
            }
            chameleon_type(): any(new any_der<double>(0)), type(tn<double>()) {
            }
            chameleon_type(chameleon_type const & c): any(c.any), type(c.type) {
            }
            template<typename T>
            void operator=(T const & t) {
                set(t);
            }
            //------------------- cast-ops -------------------
            #define create_cast_op(T)\
            operator T() const {\
                if(tn<T>() != type){\
                    throw std::runtime_error("wrong conversion to T");\
                }\
                return *(T const *)(any->get());\
            } 
            
            create_cast_op(std::string)
            create_cast_op(double)
            
            #undef create_cast_op
            //------------------- print & serialize-------------------
            void print(std::ostream & os = std::cout) const {
                #define print_op(T)\
                if(type == tn<T>())\
                    os << T(*this);\
                else
                
                print_op(std::string)
                print_op(double)
                    os<< "chameleon is empty";
                
                #undef print_op
            };
            template<typename Archive>
            void serialize(Archive & ar) {
                ar & type;
                
                #define serialize_op(T)\
                if(type == tn<T>()){ \
                    T t; \
                    if(Archive::type == archive_enum::output) \
                        t = T(*this); \
                    ar & t; \
                    set(t); \
                } \
                else
                
                serialize_op(std::string)
                serialize_op(double)
                    ; //last else
                
                #undef serialize_op
            }
        private:
            std::shared_ptr<any_base> any;
            type_type type;
        };
        std::ostream & operator<<(std::ostream & os, chameleon_type const & b) {
            b.print(os);
            return os;
        }
    }//end namespace detail
    ///  \brief handles argv and assigns the arguments automatically
    ///  
    ///  takes references (or not) and the name for the parameter. upon read it assigns the the incomming argv correctly.
    ///  there is an instance called parameter already created, so no need to construct one in your code
    class bash_parameter_class {
    public:
        typedef detail::chameleon_type store_type;
        typedef std::map<std::string, store_type> map_type;
        ///  \brief the only constructor
        ///  
        bash_parameter_class() {}
        ///  \brief internal variable
        ///  
        ///  @param name is the wanted name for bash. name = "energy" -> bash option = "-energy"
        ///  
        ///  internal variable that is changed on read, if option is given, otherwise 0
        void set(std::string name) {
            set(name, 0);
        }
        ///  \brief internal variable with default
        ///  
        ///  @param name is the wanted name for bash. name = "energy" -> bash option = "-energy"
        ///  @param def is the default value
        ///  
        ///  the internal variable is set to def and changed on read, if option is given
        void set(std::string name, store_type const & val) {
            dict[name] = val;
        }
        ///  \brief reads the argv
        ///  
        ///  looks for options that correspond to those set earier. if found, it changes the variable to the value
        ///  after the option
        template<typename T>
        void read(int argc, T argv) {
            //------------------- set path -------------------
            std::string dir(argv[0]);
            auto pos = dir.rfind("/");
            dir.erase(pos, dir.size() - pos);
            set("prog_dir", dir + "/");
            //------------------- if there is a bash_in.txt file -------------------
            std::ifstream ifs;
            ifs.open(dir + "/bash_in.txt");
            std::vector<std::string> argf;
            if(ifs.is_open()) {
                
                std::cout << GREEN << "bash_in.txt file found" << NONE << std::endl;
                while(ifs) {
                    argf.push_back("");
                    ifs >> argf.back();
                }
                ifs.close();
            }
            //------------------- join manual input -------------------
            for(int i = 0; i < argc; ++i) {
                argf.push_back(argv[i]);
            }
            
            for(unsigned i = 0; i < argf.size(); ++i) {
                if(argf[i][0] == '-') {
                    if(i + 1 < argf.size() and argf[i+1][0] != '-') {
                        #define convert_op(T)\
                        if(detail::convertable_to<T>(argf[i+1]))\
                            set(argf[i].erase(0, 1), detail::convert<T>(argf[i+1]));\
                        else
                        
                        convert_op(double)
                        convert_op(std::string)
                        std::cout << "couldn't convert to anything" << std::endl;
                        
                        #undef convert_op
                    }
                    else
                        set(argf[i].erase(0, 1), "not_set");
                }
            }
            
            
            
            
        }
        ///  \brief returns the dict
        map_type get() const {
            return dict;
        }
        ///  \brief returns the value to key
        ///  
        ///  just an abreviation for parameter.get()["key"]
        store_type & operator[](std::string const & key) {
            //~ if(dict.find(key) != dict.end())
                return dict[key];
            //~ else
                //~ throw std::runtime_error("key not found in bash_parameter");
        }
        bool contains(std::string const & key) const {
            if(dict.find(key) != dict.end())
                return true;
            return false;
        }
        template<typename Archive>
        void serialize(Archive & ar) {
            ar & dict;
        }
        void print() {
            std::for_each(dict.begin(), dict.end(), 
                [&](std::pair<std::string const, store_type const &> p) {
                    std::cout << p.first << ": " << p.second << std::endl;
                }
            );
        }
    private:
        map_type dict;    ///< the dictionary
    } parameter;
}//end namespace addon

#endif //__BASH_PARAMETER_HEADER_MSK
