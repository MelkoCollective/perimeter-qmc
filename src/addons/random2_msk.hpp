// Author:  Mario S. Könz <mskoenz@gmx.net>
// Date:    18.04.2013 20:36:17 EDT
// File:    random2_msk.hpp

#ifndef __RANDOM2_MSK_HEADER
#define __RANDOM2_MSK_HEADER

/* minimal code

#include <random2_msk.hpp>
addon::random_class<int, mersenne> rng(0, 6);
rng.seed(0);
rng();

*/

#include <boost/random.hpp>
#include <time.h>

//timer2_msk.hpp documents addon
namespace addon
{
    //+-------------------------------------------------+
    //|                 flags for choose                |
    //+-------------------------------------------------+
    ///  \brief option for random_class
    ///  
    ///  uses the mersenne twister boost::mt19937
    struct mersenne {
        typedef boost::mt19937 type;
    };
    ///  \brief option for random_class
    ///  
    ///  uses the mersenne twister boost::mt11213b (fastes according to boost)
    struct mersenne_fast {
        typedef boost::mt11213b type;
    };
    ///  \brief option for random_class
    ///  
    ///  uses the lagged fibonacci boost::lagged_fibonacci44497
    struct fibonacci {
        typedef boost::lagged_fibonacci44497 type;
    };
    
    ///  \brief handy rng class
    ///  
    ///  @tparam T is int or double, dependig on what rng you need
    ///  @tparam RNG is the option for the rng. it has to be a struct that has a typename type that is a valid rng.
    ///  
    ///  This class is basically a compiletime wrapper. it allows for an easy change of the rng. it also remaps
    ///  the native rng output (fibonacci: double [0-1), mersenne uint64_t [whole range]) to the desired range.
    ///  for optimality it uses a const int called shift that is defined via the constructor. with this variable
    ///  it's possible to avoid rescaling if it's not needed.
    template<typename T, typename RNG>
    class random_class {
        public:
            ///  \brief range [0, 1)
            ///  
            ///  doesn't make much sense for int-rng, since it will always be 0
            random_class(): shift(0) {
                init();
            }
            ///  \brief range [0-end)
            ///  
            ///  remaps the native rng output to [0, end). requires an additional multiplication.
            random_class(T const & end): scale(end), offset(0), shift(1) {
                init();
            }
            ///  \brief range [start, end)
            ///  
            ///  remaps the native rng output to [start, end). requires an additional multiplication and addition.
            random_class(T const & start, T const & end): scale(end-start), offset(start), shift(2) {
                init();
            }
            ///  \brief returns a random number
            ///  
            ///  here the const shift variable makes sure, that the minimal amount of operations needed is executed
            inline T operator()() {
                if (shift == 0) {
                    return impl_picker<RNG, int>()(rng());
                }
                if (shift == 1) {
                    return scale*impl_picker<RNG, int>()(rng());
                }
                return offset + scale*impl_picker<RNG, int>()(rng());
            }
            ///  \brief rescale
            ///  
            ///  one can change the scale during runtime. so the range [0, end) changes to [0, scale)
            inline void set_scale(T const & scaleIn) {
                scale = scaleIn;
            }
            ///  \brief set seed
            ///  
            ///  one can set an own seed.
            void seed(long int const & sd) {
                rng.seed(sd);
                seed_ = sd;
            }
            ///  \brief get seed
            ///  
            ///  you can get the automatic seed and store it, in case something strange/nice happens
            inline unsigned long int seed() {
                return seed_;
            }
            
        private:
            ///  \brief initializes rng
            ///  
            ///  first it sets the seed to time(NULL) + prevent_correlations. the second static term is increased, 
            ///  every time init is called. this is neccesary b.c. two rngs of the same type (mersenne) created in
            ///  the same second will have the same seed, thus identical.
            ///  it also burns the first 100 numbers. since it should be a chaotic system the a small difference in the seed
            ///  doesn't matter anymore after 100 calls
            inline void init() {
                seed_ = time(NULL) + prevent_correlations;
                ++prevent_correlations;
                rng.seed(seed_);
                for(uint i = 0; i < 100; ++i) {
                    rng(); //warming up the rng
                    //fibbonaci has for several seconds the same first number even though the seed changes every second
                }
            }
            ///  \brief rescale native
            ///  
            ///  the different rngs have different outputs depending on the type.
            ///  this impl_picker needs to be specialized for every rng that doesn't output a double in [0,1)
            template<typename RNG2, typename U>
            struct impl_picker {
                inline double operator()(double const r) {
                    return r;
                }
            };
            ///  \brief spec fo mersenne
            ///  
            ///  remaps the native random uint32_t to a double in [0, 1)
            template<typename U>
            struct impl_picker<mersenne, U> {
                inline double operator()(double const r) {
                    return r/uint32_t(-1);
                }
            };
            ///  \brief spec fo fast_mersenne
            ///  
            ///  remaps the native random uint32_t to a double in [0, 1)
            template<typename U>
            struct impl_picker<mersenne_fast, U> {
                inline double operator()(double const r) {
                    return r/uint32_t(-1);
                }
            };
            
            typename RNG::type rng; ///< the native rng
            T scale;    ///< length of the range
            T offset;   ///< start of the range
            const short int shift;  ///< shows, what operations are needed (for speedup. Tested!)
            uint64_t seed_; ///< the used seed
            static uint64_t prevent_correlations;   ///< needed for multiple same-type-rngs
    };
    
    template<typename T, typename RNG>
    uint64_t random_class<T, RNG>::prevent_correlations = 0;
    
}//end namespace addon

#endif //__RANDOM2_MSK_HEADER