/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                       */
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

#ifndef _PAPILO_CORE_STATISTICS_HPP_
#define _PAPILO_CORE_STATISTICS_HPP_

namespace papilo
{

struct Statistics
{
   double presolvetime;
   int ntsxapplied;
   int ntsxconflicts;
   int nboundchgs;
   int nsidechgs;
   int ncoefchgs;
   int nrounds;
   int ndeletedcols;
   int ndeletedrows;
   int consecutive_rounds_of_only_boundchanges;
   // variable substitutions and constraint deletions are excluded
   int single_matrix_coefficient_changes;

   Statistics( double _presolvetime, int _ntsxapplied, int _ntsxconflicts,
               int _nboundchgs, int _nsidechgs, int _ncoefchgs, int _nrounds,
               int _ndeletedcols, int _ndeletedrows, int _consecutive_rounds_of_only_boundchanges,
               int _single_matrix_coefficient_changes = 0
                   )
       : presolvetime( _presolvetime ), ntsxapplied( _ntsxapplied ),
         ntsxconflicts( _ntsxconflicts ), nboundchgs( _nboundchgs ),
         nsidechgs( _nsidechgs ), ncoefchgs( _ncoefchgs ), nrounds( _nrounds ),
         ndeletedcols( _ndeletedcols ), ndeletedrows( _ndeletedrows ),
         consecutive_rounds_of_only_boundchanges(_consecutive_rounds_of_only_boundchanges),
         single_matrix_coefficient_changes( _single_matrix_coefficient_changes)
   {
   }

   Statistics()
       : presolvetime( 0.0 ), ntsxapplied( 0 ), ntsxconflicts( 0 ),
         nboundchgs( 0 ), nsidechgs( 0 ), ncoefchgs( 0 ), nrounds( 0 ),
         ndeletedcols( 0 ), ndeletedrows( 0 ),
         consecutive_rounds_of_only_boundchanges( 0 ),
         single_matrix_coefficient_changes( 0 )
   {
   }
};

inline Statistics
operator-( const Statistics& a, const Statistics& b )
{
   return {
       0.0, a.ntsxapplied - b.ntsxapplied, a.ntsxconflicts - b.ntsxconflicts,
       a.nboundchgs - b.nboundchgs, a.nsidechgs - b.nsidechgs,
       a.ncoefchgs - b.ncoefchgs, a.nrounds - b.nrounds,
       a.ndeletedcols - b.ndeletedcols, a.ndeletedrows - b.ndeletedrows,
       a.consecutive_rounds_of_only_boundchanges - b.consecutive_rounds_of_only_boundchanges,
       a.single_matrix_coefficient_changes - b.single_matrix_coefficient_changes
   };
}

} // namespace papilo

#endif
