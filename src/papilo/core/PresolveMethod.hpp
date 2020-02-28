/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020  Konrad-Zuse-Zentrum                                   */
/*                     fuer Informationstechnik Berlin                       */
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

#include "papilo/core/PresolveOptions.hpp"
#include "papilo/core/Reductions.hpp"
#include "papilo/core/RowFlags.hpp"
#include "papilo/core/VariableDomains.hpp"
#include "papilo/io/Message.hpp"
#include "papilo/misc/Num.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/fmt.hpp"
#include "papilo/misc/tbb.hpp"
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

/// result codes of a presolving routine
enum class PresolveStatus : int
{
   /// problem was not changed
   kUnchanged = 0,

   /// problem was reduced
   kReduced = 1,

   /// problem was detected to be unbounded or infeasible
   kUnbndOrInfeas = 2,

   /// problem was detected to be unbounded
   kUnbounded = 3,

   /// problem was detected to be infeasible
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
      type = PresolverType::kAllCols;
      timing = PresolverTiming::kExhaustive;
      delayed = false;
      execTime = 0.0;
      enabled = true;
      skip = 0;
      nconsecutiveUnsuccessCall = 0;
   }

   virtual ~PresolveMethod() {}

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
        const Num<REAL>& num, Reductions<REAL>& reductions )
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

      auto start = tbb::tick_count::now();
      PresolveStatus result =
          execute( problem, problemUpdate, num, reductions );
      auto end = tbb::tick_count::now();
      auto duration = end - start;
      execTime = execTime + duration.seconds();

      switch( result )
      {
      case PresolveStatus::kUnbounded:
      case PresolveStatus::kUnbndOrInfeas:
      case PresolveStatus::kInfeasible:
         Message::debug( &problemUpdate,
                         "[{}:{}] {} detected unboundedness or infeasibility\n",
                         __FILE__, __LINE__, this->name );
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

   PresolverType
   getType() const
   {
      return this->type;
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

   bool
   runInRound( int roundCounter )
   {
      assert( roundCounter < 4 );

      if( ( roundCounter == 0 && timing == PresolverTiming::kFast ) ||
          ( roundCounter == 1 && timing == PresolverTiming::kMedium ) ||
          ( roundCounter == 2 && timing == PresolverTiming::kExhaustive ) )
         return true;

      // always finish with a fast round
      if( roundCounter == 3 && timing == PresolverTiming::kFast )
         return true;

      return false;
   }

   /// todo interface for certificate function and flag to indicate whther
   /// presolver writes certificates
   virtual void
   getCertificate()
   {
   }

   unsigned int
   getNCalls() const
   {
      return ncalls;
   }

   void
   setDelayed( bool delayed )
   {
      this->delayed = delayed;
   }

   void
   setEnabled( bool enabled )
   {
      this->enabled = enabled;
   }

 protected:
   /// execute member function for a presolve method gets the constant problem
   /// and can communicate reductions via the given reductions object
   virtual PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions ) = 0;

   void
   setName( const std::string& name )
   {
      this->name = name;
   }

   void
   setTiming( PresolverTiming timing )
   {
      this->timing = timing;
   }

   void
   setType( PresolverType type )
   {
      this->type = type;
   }

   void
   skipRounds( unsigned int nrounds )
   {
      this->skip += nrounds;
   }

 private:
   std::string name;
   double execTime;
   bool enabled;
   bool delayed;
   PresolverTiming timing;
   PresolverType type;
   unsigned int ncalls;
   // number of times execute returns REDUCED
   unsigned int nsuccessCall;
   unsigned int nconsecutiveUnsuccessCall;
   unsigned int skip;
};

} // namespace papilo

#endif
