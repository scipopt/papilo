/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2026 Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/* You should have received a copy of the Apache-2.0 license                 */
/* along with PaPILO; see the file LICENSE. If not visit scipopt.org.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _PAPILO_MISC_VERSION_LOGGER_HPP_
#define _PAPILO_MISC_VERSION_LOGGER_HPP_

#include "papilo/Config.hpp"
#include "papilo/core/Presolve.hpp"
#include "papilo/core/postsolve/Postsolve.hpp"
#include "papilo/io/MpsParser.hpp"
#include "papilo/io/MpsWriter.hpp"
#include "papilo/io/SolParser.hpp"
#include "papilo/io/SolWriter.hpp"
#include "papilo/misc/NumericalStatistics.hpp"
#include "papilo/misc/OptionsParser.hpp"
#if PAPILO_HAS_TBB
#include "papilo/misc/tbb.hpp"
#endif

#if PAPILO_HAS_SOPLEX
#include "soplex.h"
#include "soplex/spxgithash.h"
#endif /* PAPILO_HAS_SOPLEX */

#if PAPILO_HAS_HIGHS
#include "Highs.h"
#endif /* PAPILO_HAS_HIGHS */

#if PAPILO_HAS_SCIP
#include "scip/scipgithash.h"
#endif /* PAPILO_HAS_SCIP */

#if PAPILO_HAS_GUROBI
#include "gurobi_c++.h"
#endif /* PAPILO_HAS_GUROBI */

#if PAPILO_HAS_GLOP
#include "ortools/base/version.h"
#endif /* PAPILO_HAS_GLOP */

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <string>
#include <utility>

namespace papilo
{

void
join( const Vec<std::string>& v, char c, std::string& s )
{

   s.clear();

   for( auto p = v.begin(); p != v.end(); ++p )
   {
      s += *p;
      if( p != v.end() - 1 )
         s += c;
   }
}

void
print_header()
{
   std::string mode = "optimized";
#ifndef NDEBUG
   mode = "debug";
#endif

   Vec<std::string> list_of_solvers;
   std::string solvers = "";
#if PAPILO_HAS_HIGHS
   list_of_solvers.push_back( "HiGHS" );
#endif /* PAPILO_HAS_HIGHS */
#if PAPILO_HAS_SCIP
   list_of_solvers.push_back( "SCIP" );
#endif /* PAPILO_HAS_SCIP */
#if PAPILO_HAS_GUROBI
   list_of_solvers.push_back( "Gurobi" );
#endif /* PAPILO_HAS_GUROBI */
#if PAPILO_HAS_GLOP
   list_of_solvers.push_back( "Ortools" );
#endif /* PAPILO_HAS_GLOP */
#if PAPILO_HAS_SOPLEX
   list_of_solvers.push_back( "SoPlex" );
#endif /* PAPILO_HAS_SOPLEX */
#if PAPILO_HAS_ROUNDINGSAT
   list_of_solvers.push_back( "ROUNDINGSAT" );
#endif /* PAPILO_HAS_ROUNDINGSAT */
   if( list_of_solvers.empty() )
      solvers = "none";
   else
      join( list_of_solvers, ',', solvers );
#if PAPILO_HAS_GITHASH
   fmt::print( "PaPILO version {}.{}.{} [mode: {}][Solvers: {}][GitHash: {}]\n",
               PAPILO_VERSION_MAJOR, PAPILO_VERSION_MINOR, PAPILO_VERSION_PATCH,
               mode, solvers, PAPILO_GITHASH );
#else
   fmt::print( "PaPILO version {}.{}.{} [mode: {}][Solvers: {}][GitHash: ]\n",
               PAPILO_VERSION_MAJOR, PAPILO_VERSION_MINOR, PAPILO_VERSION_PATCH,
               mode, solvers );
#endif
   fmt::print( "Copyright (C) 2020-2026 Zuse Institute Berlin (ZIB)\n" );
   fmt::print( "\n" );

   fmt::print( "External libraries: \n" );

   fmt::print( "  Boost    {}.{}.{} \t (https://www.boost.org/)\n",
               BOOST_VERSION_NUMBER_MINOR( BOOST_VERSION ),
               BOOST_VERSION_NUMBER_PATCH( BOOST_VERSION ) / 100,
               BOOST_VERSION_NUMBER_MAJOR( BOOST_VERSION ) );

#if PAPILO_HAS_TBB
   fmt::print( "  TBB      {}.{}  \t Thread building block https://github.com/oneapi-src/oneTBB developed by Intel\n", TBB_VERSION_MAJOR, TBB_VERSION_MINOR);
#endif
#if PAPILO_HAS_GMP
   fmt::print( "  GMP      {}  \t GNU Multiple Precision Arithmetic Library "
               "developed by T. Granlund (gmplib.org)\n",
               PAPILO_GMP_VERSION );
#endif
#if PAPILO_HAS_HIGHS
   fmt::print( "  HiGHS    {}.{}.{} \t high performance software "
               "for linear optimization (https://www.maths.ed.ac.uk/hall/HiGHS/) [GitHash: {}]\n" , HIGHS_VERSION_MAJOR, HIGHS_VERSION_MINOR, HIGHS_VERSION_PATCH,
               HIGHS_GITHASH);
#endif /* PAPILO_HAS_HIGHS */
#if PAPILO_HAS_ROUNDINGSAT
    fmt::print( "  RoundingSAT    \t SAT Solver "
                "(https://gitlab.com/MIAOresearch/software/roundingsat)\n");
#endif /* PAPILO_HAS_ROUNDINGSAT */
    //TODO
#if PAPILO_HAS_GLOP
    fmt::print( "  ORTOOLS  {}.{}   \t fast and portable software for combinatorial optimization developed by Google.\n" ,
                operations_research::OrToolsMajorVersion(), operations_research::OrToolsMinorVersion() );
#endif /* PAPILO_HAS_GLOP */

#if PAPILO_HAS_SCIP
   fmt::print( "  SCIP     {}.{}.{} \t Mixed Integer Programming Solver "
               "developed at Zuse "
               "Institute Berlin (scip.zib.de) [GitHash: {}]\n",
               SCIP_VERSION_MAJOR, SCIP_VERSION_MINOR, SCIP_VERSION_PATCH,
               SCIPgetGitHash() );
#endif /* PAPILO_HAS_SCIP */
#if PAPILO_HAS_GUROBI
   fmt::print( "  Gurobi   {}.{}.{} \t Gurobi Optimizer \n",
               GRB_VERSION_MAJOR, GRB_VERSION_MINOR, GRB_VERSION_TECHNICAL );
#endif /* PAPILO_HAS_GUROBI */
#if PAPILO_HAS_SOPLEX
   fmt::print(
       "  SoPlex   {}.{}.{} \t Linear Programming Solver developed at Zuse "
       "Institute Berlin (soplex.zib.de) [GitHash: {}]\n",
       SOPLEX_VERSION / 100, ( SOPLEX_VERSION % 100 ) / 10, SOPLEX_VERSION % 10,
       soplex::getGitHash() );
#endif /* PAPILO_HAS_SOPLEX */


   fmt::print( "\n" );
}



} // namespace papilo

#endif
