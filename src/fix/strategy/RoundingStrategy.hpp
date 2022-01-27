#include "papilo/core/ProbingView.hpp"

#include <cassert>
#include <random>

using namespace papilo;

template <typename REAL>
class RoundingStrategy
{
 public:
   virtual Fixing<REAL>
   select_diving_variable( const Vec<REAL>& cont_solution,
                           const ProbingView<REAL>& view ) = 0;
};

template <typename REAL>
class RandomRoundingStrategy : public RoundingStrategy<REAL>
{

   const Num<REAL> num;

   typedef std::mt19937 MyRNG;
   uint32_t seed;

   MyRNG random_generator;

 public:
   RandomRoundingStrategy( uint32_t seed_, Num<REAL> num_ )
       : seed( seed_ ), num( num_ )
   {
      random_generator.seed( seed );
   }

   Fixing<REAL>
   select_diving_variable( const Vec<REAL>& cont_solution,
                           const ProbingView<REAL>& view ) override
   {
      // TODO: this does not work since fixed variable could be obtained

      Vec<int> remaining_unfixed_cols{};
      for( int i = 0; i < cont_solution.size(); i++ )
      {
         if( num.isEq( view.getProbingUpperBounds()[i],
                       view.getProbingLowerBounds()[i] ) )
            continue;
         remaining_unfixed_cols.push_back( i );
      }
      if( remaining_unfixed_cols.empty() )
         return { -1, -1 };

      std::uniform_int_distribution<uint32_t> dist_variable(
          0, remaining_unfixed_cols.size() - 1 );
      std::uniform_int_distribution<uint32_t> dist_rounding( 0, 1 );
      int variable = remaining_unfixed_cols[dist_variable( random_generator )];
      REAL value = -1;
      if( dist_rounding( random_generator ) )
         value = round( cont_solution[variable] + 0.5 );
      else
         value = round( cont_solution[variable] - 0.5 );

      return { variable, value };
   }
};

template <typename REAL>
class FractionalRoundingStrategy : public RoundingStrategy<REAL>
{

   const Num<REAL> num;

 public:
   FractionalRoundingStrategy( Num<REAL> num_ ) : num( num_ ) {}

   Fixing<REAL>
   select_diving_variable( const Vec<REAL>& cont_solution,
                           const ProbingView<REAL>& view ) override
   {

      // this is currently fractional diving
      REAL value = -1;
      int variable = -1;
      REAL score = -1;

      for( int i = 0; i < cont_solution.size(); i++ )
      {
         REAL frac = cont_solution[i] - floor( cont_solution[i] );
         if( frac == 0 || num.isEq( view.getProbingUpperBounds()[i],
                                    view.getProbingLowerBounds()[i] ) )
            continue;
         else if( frac > 0.5 )
         {
            if( variable == -1 || 1 - frac > score )
            {
               score = 1 - frac;
               variable = i;
               value = ceil( cont_solution[i] );
            }
         }
         else
         {
            if( variable == -1 || frac > score )
            {
               score = frac;
               variable = i;
               value = floor( cont_solution[i] );
            }
         }
      }
      return { variable, value };
   }
};