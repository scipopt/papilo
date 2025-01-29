/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software: you can redistribute it and/or modify      */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation, either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program.  If not, see <https://www.gnu.org/licenses/>.    */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _PAPILO_CORE_PRESOLVE_METHOD_HPP_
#define _PAPILO_CORE_PRESOLVE_METHOD_HPP_

#include "papilo/Config.hpp"
#include "papilo/core/ConstraintMatrix.hpp"
#include "papilo/core/PresolveOptions.hpp"
#include "papilo/core/Reductions.hpp"
#include "papilo/core/RowFlags.hpp"
#include "papilo/core/VariableDomains.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Timer.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/fmt.hpp"
#ifdef PAPILO_TBB
#include "papilo/misc/tbb.hpp"
#else
#include <chrono>
#endif

#include "papilo/verification/ArgumentType.hpp"
#include <bitset>


namespace papilo
{

// forward declaration of problem and reduction
template <typename REAL>
class Presolve;

template <typename REAL>
class Problem;

template <typename REAL>
class ProblemUpdate;

enum class PresolveStatus : int
{
   kUnchanged = 0,

   kReduced = 1,

   kUnbndOrInfeas = 2,

   kUnbounded = 3,

   kInfeasible = 4,
};

enum class PresolverTiming : int
{
   kFast = 0,
   kMedium = 1,
   kExhaustive = 2,
};

enum class PresolverType
{
   kAllCols,
   kIntegralCols,
   kContinuousCols,
   kMixedCols,
};

template <typename REAL>
class PresolveMethod
{
 public:
   PresolveMethod()
   {
      ncalls = 0;
      nsuccessCall = 0;
      name = "unnamed";
      argument = ArgumentType::kPrimal;
      type = PresolverType::kAllCols;
      timing = PresolverTiming::kExhaustive;
      execTime = 0.0;
      enabled = true;
      delayed = false;
      symmetries_active = true;
      skip = 0;
      nconsecutiveUnsuccessCall = 0;
   }

   virtual ~PresolveMethod() = default;

   virtual void
   compress( const Vec<int>& rowmap, const Vec<int>& colmap )
   {
   }

   virtual bool
   initialize( const Problem<REAL>& problem,
               const PresolveOptions& presolveOptions )
   {
      return false;
   }

   virtual void
   addPresolverParams( ParameterSet& paramSet )
   {
   }

   void
   addParameters( ParameterSet& paramSet )
   {
      paramSet.addParameter(
          fmt::format( "{}.enabled", this->name ).c_str(),
          fmt::format( "is presolver {} enabled", this->name ).c_str(),
          this->enabled );

      addPresolverParams( paramSet );
   }

   PresolveStatus
   run( const Problem<REAL>& problem, const ProblemUpdate<REAL>& problemUpdate,
        const Num<REAL>& num, Reductions<REAL>& reductions, const Timer& timer, int& cause )
   {
      if( !enabled || delayed )
         return PresolveStatus::kUnchanged;

      if( skip != 0 )
      {
         --skip;
         return PresolveStatus::kUnchanged;
      }

      if( problem.getNumIntegralCols() == 0 &&
          ( type == PresolverType::kIntegralCols ||
            type == PresolverType::kMixedCols ) )
         return PresolveStatus::kUnchanged;

      if( problem.getNumContinuousCols() == 0 &&
          ( type == PresolverType::kContinuousCols ||
            type == PresolverType::kMixedCols ) )
         return PresolveStatus::kUnchanged;

      ++ncalls;

#ifdef PAPILO_TBB
      auto start = tbb::tick_count::now();
#else
      auto start = std::chrono::steady_clock::now();
#endif
      PresolveStatus result =
          execute( problem, problemUpdate, num, reductions, timer, cause );
#ifdef PAPILO_TBB
      auto end = tbb::tick_count::now();
      auto duration = end - start;
      execTime = execTime + duration.seconds();
#else
      auto end = std::chrono::steady_clock::now();
      execTime = execTime + std::chrono::duration_cast<std::chrono::milliseconds>(
                                end- start ).count()/1000;
#endif


      switch( result )
      {
      case PresolveStatus::kUnbounded:
      case PresolveStatus::kUnbndOrInfeas:
      case PresolveStatus::kInfeasible:
         Message::debug( &problemUpdate,
                         "[{}:{}] {} detected unboundedness or infeasibility\n",
                         __FILE__, __LINE__, this->name );
         break;
      case PresolveStatus::kReduced:
         ++nsuccessCall;
         nconsecutiveUnsuccessCall = 0;
         break;
      case PresolveStatus::kUnchanged:
         ++nconsecutiveUnsuccessCall;
         if( timing != PresolverTiming::kFast )
            skip += nconsecutiveUnsuccessCall;
         break;
      }

      return result;
   }

   void
   printStats( const Message& message, std::pair<int, int> stats )
   {
      double success =
          ncalls == 0 ? 0.0
                      : ( double( nsuccessCall ) / double( ncalls ) ) * 100.0;
      double applied =
          stats.first == 0
              ? 0.0
              : ( double( stats.second ) / double( stats.first ) ) * 100.0;
      message.info( " {:>18} {:>12} {:>18.1f} {:>18} {:>18.1f} {:>18.3f}\n",
                    name, ncalls, success, stats.first, applied, execTime );
   }

   PresolverTiming
   getTiming() const
   {
      return this->timing;
   }

   bool
   isEnabled() const
   {
      return this->enabled;
   }

   bool
   isDelayed() const
   {
      return this->delayed;
   }

   const std::string&
   getName() const
   {
      return this->name;
   }

   ArgumentType
   getArgument() const
   {
      return this->argument;
   }

   unsigned int
   getNCalls() const
   {
      return ncalls;
   }

   void
   setDelayed( bool value )
   {
      this->delayed = value;
   }

   void
   setEnabled( bool value )
   {
      this->enabled = value;
   }

   void
   set_symmetries_enabled( bool value )
   {
      this->symmetries_active = value;
   }

   PresolveStatus
   run_symmetries( const Problem<REAL>& problem, const ProblemUpdate<REAL>& problemUpdate,
                    const Num<REAL>& num, Reductions<REAL>& reductions, const Timer& timer )
   {
      if( !enabled || !symmetries_active )
         return PresolveStatus::kUnchanged;

#ifdef PAPILO_TBB
      auto start = tbb::tick_count::now();
#else
      auto start = std::chrono::steady_clock::now();
#endif
      PresolveStatus result =
          execute_symmetries( problem, problemUpdate, num, reductions, timer );
#ifdef PAPILO_TBB
      auto end = tbb::tick_count::now();
      auto duration = end - start;
      execTime = execTime + duration.seconds();
#else
      auto end = std::chrono::steady_clock::now();
      execTime = execTime + std::chrono::duration_cast<std::chrono::milliseconds>(
                                end- start ).count()/1000;
#endif

      return result;
   }



 protected:
   /// execute member function for a presolve method gets the constant problem
   /// and can communicate reductions via the given reductions object
   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate,
            const Num<REAL>& num, Reductions<REAL>& reductions,
            const Timer& timer, int& reason_of_infeasibility) = 0;

   virtual
   PresolveStatus
   execute_symmetries( const Problem<REAL>& problem,
                       const ProblemUpdate<REAL>& problemUpdate,
                       const Num<REAL>& num, Reductions<REAL>& reductions,
                       const Timer& timer )
   {
      return PresolveStatus::kUnchanged;
   }

   void
   setName( const std::string& value )
   {
      this->name = value;
   }

   void
   setArgument( ArgumentType value )
   {
      this->argument = value;
   }


   void
   setTiming( PresolverTiming value )
   {
      this->timing = value;
   }

   void
   setType( PresolverType value )
   {
      this->type = value;
   }

   static bool
   is_time_exceeded( const Timer& timer,  double tlim)
   {
      return tlim != std::numeric_limits<double>::max() &&
             timer.getTime() >= tlim;
   }

   bool
   check_if_substitution_generates_huge_or_small_coefficients( const Num<REAL>& num,
                                       const ConstraintMatrix<REAL>& constMatrix,
                                       int equalityrow, int col )
   {
      assert(col < constMatrix.getNCols());
      assert(equalityrow < constMatrix.getNRows());
      assert(constMatrix.getRowFlags()[equalityrow].test(RowFlag::kEquation));

      const auto colvec = constMatrix.getColumnCoefficients( col );
      const int* colindices = colvec.getIndices();

      const auto rowvec = constMatrix.getRowCoefficients( equalityrow );
      const int row_length = rowvec.getLength();
      const int* row_indices = rowvec.getIndices();
      REAL coeff = 0;
      REAL min_col_vector = -1;
      REAL max_col_vector = 0;
      for( int i = 0; i < colvec.getLength(); i++ )
      {
         REAL val = colvec.getValues()[i];
         if( colindices[i] == equalityrow )
            coeff = val;
         if( min_col_vector == -1 || abs( min_col_vector ) > abs( val ) )
            min_col_vector = abs( val );
         if( abs( max_col_vector ) < abs( val ) )
            max_col_vector = abs( val );
      }
      assert( coeff != 0 );

      if( num.isHugeVal( coeff ) )
         return false;

      REAL min_row_vector = -1;
      REAL max_row_vector = 0;
      for( int j = 0; j < row_length; j++ )
      {
         REAL val = rowvec.getValues()[j];
         if( row_indices[j] == col )
            continue;
         if( min_row_vector == -1 || abs( min_row_vector ) > abs( val ) )
            min_row_vector = abs( val );
         if( abs( max_row_vector ) < abs( val ) )
            max_row_vector = abs( val );
      }
      if( num.isHugeVal( max_col_vector ) || num.isHugeVal( max_row_vector ) ||
          num.isHugeVal( max_col_vector / coeff * max_row_vector ) ||
          num.isZero( min_col_vector / coeff * min_row_vector ) )
      {
         return false;
      }
      return true;
   }

   void
   skipRounds( unsigned int nrounds )
   {
      this->skip += nrounds;
   }

   template <typename LOOP>
   void
   loop( int start, int end, LOOP&& loop_instruction )
   {
#ifdef PAPILO_TBB
      tbb::parallel_for( tbb::blocked_range<int>( start, end ),
                         [&]( const tbb::blocked_range<int>& r )
                         {
                            for( int i = r.begin(); i != r.end(); ++i )
                               loop_instruction( i );
                         } );
#else
      for( int i = 0; i < end; i++ )
      {
         loop_instruction( i );
      }
#endif
   }

 private:
   std::string name;
   ArgumentType argument;
   double execTime;
   bool enabled;
   bool delayed;
   bool symmetries_active;
   PresolverTiming timing;
   PresolverType type;
   unsigned int ncalls;
   // number of times execute returns REDUCED
   unsigned int nsuccessCall;
   unsigned int nconsecutiveUnsuccessCall;
   unsigned int skip;
   };

} // namespace papilo

// for MSVS 2022, the handling of operator<< for papilo::PresolveStatus fails
// instead, one gets some cryptic compiling error
// this adds an specific implementation as a workaround
#ifdef _MSC_VER
inline
std::ostream& operator<<(std::ostream& os, const papilo::PresolveStatus& p)
{
   os << (int)p;
   return os;
}
#endif

#endif
