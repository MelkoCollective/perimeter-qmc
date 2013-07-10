// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    18.04.2013 21:03:11 EDT
// File:    timer2_msk.hpp

#ifndef __TIMER2_MSK_HEADER
#define __TIMER2_MSK_HEADER

//~ #define __TIMER2_MSK_USE_BOOST_CLOCK
#define __TIMER2_MSK_USE_C_CLOCK

//~ #define __TIMER2_MSK_USE_C11_CLOCK
//~ #define __TIMER2_MSK_WHAT_C11_CLOCK std::chrono::system_clock
//~ #define __TIMER2_MSK_WHAT_C11_CLOCK std::chrono::high_resolution_clock

#ifdef __TIMER2_MSK_USE_BOOST_CLOCK
    #include <boost/timer/timer.hpp>
#elif defined(__TIMER2_MSK_USE_C11_CLOCK)
    #include <chrono>
#endif

#include <ctime>

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <typeinfo>
#include <stdexcept>

#include "color.hpp"

/* minimal code

#include <timer2_msk.hpp>
#include "timer2_msk.hpp"

//addon::timer_class<addon::normal> timer(workload, filename);
addon::timer_class<addon::normal> timer;
timer.set_names("constant-1", "constant-2");
timer.set_comment("test");
timer.print(1, 2);
timer.write(1, 2);
* 
*/

///  \brief usefull tools
///  
///  diverse tools like timer, random number generator and bash_parameter
namespace addon {
    // +---------------------------------------------------+
    // |                   constants                       |
    // +---------------------------------------------------+
    const std::string default_name = "benchmark.txt";
    const std::string default_time_format = "%d.%m.%y at %H:%M:%S   "; //check strftime of time.h

    ///  \brief prints date and time
    ///  
    ///  the format is "dd.mm.yy at hh:mm:ss"
    void timer_date(std::ostream & os) {
        time_t t;
        time(&t);
        char buffer[80];
        struct tm * timeinfo;
        timeinfo = localtime(&t);
        
        strftime(buffer, 80, default_time_format.c_str(), timeinfo);
        
        os << buffer;
    }
    ///  \brief prints date and time a little different
    ///  
    ///  the format is "dd.mm.yy---hh:mm"
    void timer_time(std::ostream & os) {
        time_t t;
        time(&t);
        char buffer[80];
        struct tm * timeinfo;
        timeinfo = localtime(&t);
        
        strftime(buffer, 80, "%d.%m.%y---%H:%M", timeinfo);
        
        os << buffer;
    }
    ///  \brief default template parameter
    ///  
    ///  if the user doesn't specify all template parameters, the remaining will have the type
    ///  not_used, what can be checked via a typeid compare
    struct not_used {
        not_used() {}
    };
    ///  \brief stream for not_used
    ///  since the default template argument is assumed to be streamable
    std::ostream & operator<<(std::ostream & os, not_used nu) {
        return os;
    }
    ///  \brief option for timer_class
    ///  
    ///  if choosen, the output will be in a work friendly format, but not so well readable
    struct data {
    };
    ///  \brief option for timer_class
    ///  
    ///  if choosen, the output will be readable but not so convenient to work with.
    ///  it will print the same output in the file that is shown with print(...)
    struct normal {
    };
    namespace dest {
        enum dest_enum {
              console = 0
            , file
        };
    }//end namespace dest
    
    namespace detail {
        unsigned active_timer_ = 0;
    }//end namespace detail
    ///  \brief the improved timer class
    ///  
    ///  it uses boost cpu_timer s.t. it only measures the usertime, not the walltime.
    ///  for now one has to link libboost_timer.so and libboost_system.so to the binary
    template<typename T>
    class timer_class {
    public:
        ///  \brief the only constructor
        ///  
        ///  @param workload takes a uint64_t that specifies, how much work that has to be done
        ///  @param name defines the output-filename. the default is "benchmark.txt"
        ///  the timer starts as soon as the class is constructed
        timer_class<T>(uint64_t workload = 1, std::string const & name = default_name): 
                                                              name_(name)
                                                            , work_(workload)
                                                            #ifdef __TIMER2_MSK_USE_C_CLOCK
                                                              , start_t(std::clock())
                                                            #elif defined(__TIMER2_MSK_USE_C11_CLOCK)
                                                              , start_t(__TIMER2_MSK_WHAT_C11_CLOCK::now())
                                                            #endif
                                                            , written_(0)
                                                            , comment_("")
                                                            , last_print_(0)
                                                            , last_i_(0)
                                                            , first_i_(-1)
                                                            , mod_(1)
                                                            , loop_time_(0)
                                                            , timer_id_(detail::active_timer_)
                                                            {
            ++detail::active_timer_;
        }
        ///  \brief the destructor
        ///  
        ///  it prints the finish, depending if it is in data or normal mode
        ~timer_class() {
            --detail::active_timer_;
            of_.open(name_.c_str(), std::ios_base::app);
            if(written_ > 1 or (written_ > 0 and typeid(T) == typeid(data)))
            {
                of_ << "--------------------------------------------------------";
                timer_time(of_);
            }
                
            if(written_ != 0)
                of_ << std::endl;
            of_.close();
        }
        ///  \brief returns the elapsed usertime in [s]
        operator double() const {
            return elapsed();
        }
        ///  \brief returns the last mesured loop-time in [us]
        double loop_time() {
            if(loop_time_ != 0)
                return loop_time_;
            else
                return (elapsed() * 1000000) / work_;
        }
        ///  \brief name the data
        ///  
        ///  the names cannot contain spaces. a runtime_error will be thrown is so.
        ///  in the normal-mode these names will preceed the data.
        ///  in the data-mode, this will be the title of the rows.
        ///  if the file already exists, it will check, if the names are the same, i.o. to avoid printing different
        ///  data-rows into the same file. a runtime_error will be thrown upon conflict
        void set_names(   std::string n0 = "empty", std::string n1 = "empty", std::string n2 = "empty", std::string n3 = "empty", std::string n4 = "empty"
                        , std::string n5 = "empty", std::string n6 = "empty", std::string n7 = "empty", std::string n8 = "empty", std::string n9 = "empty"
                        ) {
            if(    n0.find(' ') != std::string::npos or n1.find(' ') != std::string::npos or n2.find(' ') != std::string::npos or n3.find(' ') != std::string::npos or n4.find(' ') != std::string::npos 
                or n5.find(' ') != std::string::npos or n6.find(' ') != std::string::npos or n7.find(' ') != std::string::npos or n8.find(' ') != std::string::npos or n9.find(' ') != std::string::npos) {
                throw std::runtime_error("There is a space in a name in the method timer.set_names. Please remove spaces from names.");
            }
            if(written_ == 0) //after a first write it's impossible to change names
            {
                names_[0] = n0;
                names_[1] = n1;
                names_[2] = n2;
                names_[3] = n3;
                names_[4] = n4;
                names_[5] = n5;
                names_[6] = n6;
                names_[7] = n7;
                names_[8] = n8;
                names_[9] = n9;
                
                if(typeid(T) == typeid(data))
                {
                    std::stringstream o;
                    o << "Time[s]";
                    if(names_[0] != "empty") o << " " << names_[0];
                    if(names_[1] != "empty") o << " " << names_[1];
                    if(n2 != "empty") o << " " << names_[2];
                    if(n3 != "empty") o << " " << names_[3];
                    if(n4 != "empty") o << " " << names_[4];
                    if(n5 != "empty") o << " " << names_[5];
                    if(n6 != "empty") o << " " << names_[6];
                    if(n7 != "empty") o << " " << names_[7];
                    if(n8 != "empty") o << " " << names_[8];
                    if(n9 != "empty") o << " " << names_[9];
                    
                    check_descriptor_line(o.str());
                }
            }
        }
        double p(uint64_t const & i) const {
            return double(i) / work_;
        }
        std::string time_predict(uint64_t const & i) const {
            std::stringstream ss;
            double e = elapsed();
            double p2(double(i - first_i_) / (work_ - first_i_));
            unsigned time_pred((1-p2)/p2*e);
            ss << std::setfill('0') << std::setw(2) << time_pred/3600 << ":"
                                    << std::setw(2) << (time_pred/60)%60 << ":"
                                    << std::setw(2) << time_pred%60 << std::setfill(' ');
            return ss.str();
        }
        void write_state(uint64_t const & i) {
            int pos = name_.rfind("/");
            std::string state_f = name_;
            state_f = state_f.erase(pos, name_.size() - pos) + "/state.txt";
            std::ofstream of_;
            of_.open(state_f.c_str());
            of_ << "p " << p(i) << std::endl;
            of_ << "t " << time_predict(i) << std::endl;
            of_ << "l " << loop_time() << std::endl;
            of_.close();
        }
        ///  \brief autotuned low-overhead progress display
        ///  
        ///  @param i is the loop, variable that is assumed to be in the intervall [0-workload)
        ///  call this function in the loop that corresponds to the parameter i. It will print every 1-2s
        ///  the progress in % and an expected remaining duration (constant workload assumed). It achives this
        ///  behavior using an power 2 modulo operation (binary &) at the start of the function. The modulo-factor
        ///  is changed dynamically, s.t. the print only happens every 1-2s. The overhead is minimal since it's only
        ///  a binary & operation most of the time. (takes 20% of the time of a single mersenne-rng() call)
        void progress(uint64_t const & i, unsigned const & mode = dest::console) {
            if((i&(mod_-1)) == 0)
            {
                if(elapsed() - last_print_ > 1) {
                    double e = elapsed();
                    if(e - last_print_ > 3 and mod_ > 1)
                        mod_ >>= 1;

                    last_print_ = e;
                    last_i_ = i;
                    loop_time_ = (e * 1000000) / (i - first_i_);
                    
                    if(mode == dest::console) {
                        if(last_print_ != 0)
                            std::cout << "\x1B[1A" << "\x1B[2K" << "\x1B[80D"; //move 1 up // kill all("2K") /to left("1K") /to right("K") line // move 80 left
                        if(first_i_ == -1)
                            first_i_ = i;
                        
                        
                        for(unsigned i = 0; i < timer_id_; ++i) {
                            std::cout << "        ";
                        }
                        
                        std::cout   << "progress: " << REDB << std::setprecision(4) << std::setw(3) << int(100 * p(i)) << "% " << NONE
                                    << "  loop-time: " << std::setw(7) << std::setprecision(4) << loop_time_ << " us"
                                    << "  mod: " << std::setw(10) << mod_ 
                                    << "  left: " << GREENB << time_predict(i) << NONE
                                    << std::endl;
                    }
                    else if(mode == dest::file) {
                        write_state(i);
                    }
                }
                else {
                    if(first_i_ == -1)
                        first_i_ = i;
                    
                    if(i != first_i_)
                        mod_ <<= 1;
                    
                }
            }
        }
        template<typename F>
        void progress_trigger(uint64_t const & i, F fct) {
            if((i&((mod_<<1)-1)) == 0) {
                fct();
            }
        }
        ///  \brief burns time
        ///  
        ///  @param dt is the amount of time in seconds one wants to "wait"
        void burn(double const & dt) {
            double end = elapsed() + dt;
            while(elapsed() < end) ;
        }
        
        ///  \brief free comment
        ///  
        ///  is only used in normal-mode
        void set_comment( std::string com) {
            comment_ = com;
        }
        ///  \brief retruns the elapsed usertime in [s]
        ///  
        ///  works the same as the double cast
        double elapsed() const {
            #ifdef __TIMER2_MSK_USE_BOOST_CLOCK
                return timer_.elapsed().user / 1000000000.0;
            #elif defined(__TIMER2_MSK_USE_C_CLOCK)
                return (std::clock() - start_t) / double(CLOCKS_PER_SEC);
            #elif defined(__TIMER2_MSK_USE_C11_CLOCK)
                return std::chrono::duration_cast<std::chrono::duration<double>>(__TIMER2_MSK_WHAT_C11_CLOCK::now() - start_t).count();
            #else
                return 0;
            #endif
        }
        ///  \brief the file-writer
        ///  
        ///  writes to the file depending on the mode. the parameter can be anything that is printable.
        ///  the order of the parameters should be the same as the one in set_names
        template< typename T0 = not_used, typename T1 = not_used, typename T2 = not_used, typename T3 = not_used, typename T4 = not_used
                , typename T5 = not_used, typename T6 = not_used, typename T7 = not_used, typename T8 = not_used, typename T9 = not_used
                >
        void write(   T0 t0 = not_used(), T1 t1 = not_used(), T2 t2 = not_used(), T3 t3 = not_used(), T4 t4 = not_used()
                    , T5 t5 = not_used(), T6 t6 = not_used(), T7 t7 = not_used(), T8 t8 = not_used(), T9 t9 = not_used()
                  ) {
            if(typeid(T) == typeid(data)) {
                data_writer(t0, t1, t2, t3, t4, t5, t6, t7, t8, t9);
            }
            else {
                of_.open(name_.c_str(), std::ios_base::app);
                normal_writer(of_, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9);
                of_.close();
            }
        }
        ///  \brief prints to console
        ///  
        ///  always prints in normal-mode.
        ///  the order of the parameters should be the same as the one in set_names
        template< typename T0 = not_used, typename T1 = not_used, typename T2 = not_used, typename T3 = not_used, typename T4 = not_used
                , typename T5 = not_used, typename T6 = not_used, typename T7 = not_used, typename T8 = not_used, typename T9 = not_used
                >
        void print(   T0 t0 = not_used(), T1 t1 = not_used(), T2 t2 = not_used(), T3 t3 = not_used(), T4 t4 = not_used()
                    , T5 t5 = not_used(), T6 t6 = not_used(), T7 t7 = not_used(), T8 t8 = not_used(), T9 t9 = not_used()
                  ) {
            normal_writer(std::cout, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9);
        }
    private:
        ///  \brief data-mode only
        ///  
        ///  if the file exists, this checks if the names of the rows are identical.
        ///  runttime_error otherwise
        void check_descriptor_line(std::string const compare) {
            std::ifstream i(name_.c_str(), std::ios_base::in);
            if(i.is_open()) //file exists
            {
                std::string original = "";
                
                while(original == "" and !i.eof())
                    std::getline(i, original);
                
                if(original == compare) //all is ok, descriptor matches
                    ;
                else
                {
                    if(original != "") //the file exists, but isn't empty
                        throw std::runtime_error("Descriptor line in write_data is different. Please be consistent.");
                    else
                    { //the file exists, but is empty
                        of_.open(name_.c_str(), std::ios_base::trunc | std::ios_base::out);
                        of_ << compare;
                        of_ << std::endl;
                        of_.close();
                    }
                }
            }
            else
            {
                of_.open(name_.c_str(), std::ios_base::trunc | std::ios_base::out);
                of_ << compare;
                of_ << std::endl;
                of_.close();
            }
            
            i.close();
        }
        ///  \brief normal-mode
        ///  
        ///  prints the data to os. this is used for print and write (normal-mode)
        template< typename T0 = not_used, typename T1 = not_used, typename T2 = not_used, typename T3 = not_used, typename T4 = not_used
                , typename T5 = not_used, typename T6 = not_used, typename T7 = not_used, typename T8 = not_used, typename T9 = not_used
                >
        void normal_writer(std::ostream & os = std::cout
                    , T0 t0 = not_used(), T1 t1 = not_used(), T2 t2 = not_used(), T3 t3 = not_used(), T4 t4 = not_used()
                    , T5 t5 = not_used(), T6 t6 = not_used(), T7 t7 = not_used(), T8 t8 = not_used(), T9 t9 = not_used()
                  ) {
            timer_date(os);
            int s0 = 20;
            int s1 = 10;
            std::string endl;
            if(typeid(t4) != typeid(not_used))
                endl = "\n";
            else
            {
                endl = "  ";
                s0 = 0;
                s1 = 0;
            }
            
            os << comment_ << endl;
            
            os << std::setw(s0) << "Time[s]" << ": " << std::setw(s1) << elapsed() << endl;
            
            if(typeid(t0) != typeid(not_used))
                os << std::setw(s0) << names_[0] << ": " << std::setw(s1) << t0 << endl;
            if(typeid(t1) != typeid(not_used))
                os << std::setw(s0) << names_[1] << ": " << std::setw(s1) << t1 << endl;
            if(typeid(t2) != typeid(not_used))
                os << std::setw(s0) << names_[2] << ": " << std::setw(s1) << t2 << endl;
            if(typeid(t3) != typeid(not_used))
                os << std::setw(s0) << names_[3] << ": " << std::setw(s1) << t3 << endl;
            if(typeid(t4) != typeid(not_used))
                os << std::setw(s0) << names_[4] << ": " << std::setw(s1) << t4 << endl;
            if(typeid(t5) != typeid(not_used))
                os << std::setw(s0) << names_[5] << ": " << std::setw(s1) << t5 << endl;
            if(typeid(t6) != typeid(not_used))
                os << std::setw(s0) << names_[6] << ": " << std::setw(s1) << t6 << endl;
            if(typeid(t7) != typeid(not_used))
                os << std::setw(s0) << names_[7] << ": " << std::setw(s1) << t7 << endl;
            if(typeid(t8) != typeid(not_used))
                os << std::setw(s0) << names_[8] << ": " << std::setw(s1) << t8 << endl;
            if(typeid(t9) != typeid(not_used))
                os << std::setw(s0) << names_[9] << ": " << std::setw(s1) << t9 << endl;
            if(endl != "\n")
                os << std::endl;
            if(os != std::cout)
                ++written_;
        }
        ///  \brief data-mode
        ///  
        ///  writes the data into the file. only used by write (data-mode).
        ///  the order of the parameters should be the same as the one in set_names
        template< typename T0 = not_used, typename T1 = not_used, typename T2 = not_used, typename T3 = not_used, typename T4 = not_used
                , typename T5 = not_used, typename T6 = not_used, typename T7 = not_used, typename T8 = not_used, typename T9 = not_used
                >
        void data_writer(T0 t0 = not_used(), T1 t1 = not_used(), T2 t2 = not_used(), T3 t3 = not_used(), T4 t4 = not_used()
                     , T5 t5 = not_used(), T6 t6 = not_used(), T7 t7 = not_used(), T8 t8 = not_used(), T9 t9 = not_used()
                  ) {
            of_.open(name_.c_str(), std::ios_base::app);
            of_ << elapsed();
            if(typeid(t0) != typeid(not_used)) of_ << " " << t0;
            if(typeid(t1) != typeid(not_used)) of_ << " " << t1;
            if(typeid(t2) != typeid(not_used)) of_ << " " << t2;
            if(typeid(t3) != typeid(not_used)) of_ << " " << t3;
            if(typeid(t4) != typeid(not_used)) of_ << " " << t4;
            if(typeid(t5) != typeid(not_used)) of_ << " " << t5;
            if(typeid(t6) != typeid(not_used)) of_ << " " << t6;
            if(typeid(t7) != typeid(not_used)) of_ << " " << t7;
            if(typeid(t8) != typeid(not_used)) of_ << " " << t8;
            if(typeid(t9) != typeid(not_used)) of_ << " " << t9;
            of_ << std::endl;
            
            ++written_;
            of_.close();
        }
        
        std::string const name_;
        std::ofstream of_;
        uint64_t work_;
        
        #ifdef __TIMER2_MSK_USE_BOOST_CLOCK
            boost::timer::cpu_timer timer_;
        #elif defined(__TIMER2_MSK_USE_C_CLOCK)
            std::clock_t start_t;
        #elif defined(__TIMER2_MSK_USE_C11_CLOCK)
            __TIMER2_MSK_WHAT_C11_CLOCK::time_point start_t;
        #endif
        
        int written_;
        std::string names_[10];
        std::string comment_;
        double last_print_;
        double last_i_;
        double first_i_; //in case the simulation starts at i = 50 / 100
        uint64_t mod_;
        
        double loop_time_;
        unsigned timer_id_;
    };
}//end namespace addon
#endif //__TIMER2_MSK_HEADER
