#include "papilo/core/ProbingView.hpp"

#include <cassert>
#include <random>

using namespace papilo;

template <typename REAL>
class FarkasStrategy : public RoundingStrategy<REAL>
{

   const Num<REAL> num;

   typedef std::mt19937 MyRNG;
   uint32_t seed;

   MyRNG random_generator;

 public:
   FarkasStrategy( uint32_t seed_, Num<REAL> num_ ) : seed( seed_ ), num( num_ )
   {
      random_generator.seed( seed );
   }

   Fixing<REAL>
   select_diving_variable( const Vec<REAL>& cont_solution,
                           const ProbingView<REAL>& view ) override
   {
      // this is currently fractional diving
      REAL value = -1;
      int variable = -1;
      REAL score = -1;
      std::uniform_int_distribution<int> dist_rounding( 0, 1e5 );
      const Vec<REAL>& obj = view.get_obj();

      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isEq( view.getProbingUpperBounds()[i],
                       view.getProbingLowerBounds()[i] ) )
            continue;

         REAL current_score = abs( obj[i] ) + dist_rounding( random_generator );

         if( current_score > score )
         {
            if( num.isLT( obj[i], 0 ) )
            {
               variable = i;
               value = ceil( cont_solution[i] );
            }
            else if( num.isGT( obj[i], 0 ) )
            {
               variable = i;
               value = floor( cont_solution[i] );
            }
            else
            {
               REAL frac = cont_solution[i] - floor( cont_solution[i] );
               variable = i;
               if( frac > 0.5 )
                  value = ceil( cont_solution[i] );
               else
                  value = floor( cont_solution[i] );
            }
         }
      }
      return { variable, value };
   }
};
