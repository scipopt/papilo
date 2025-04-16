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

#ifndef _PAPILO_CORE_REDUCTION_TYPE_HPP_
#define _PAPILO_CORE_REDUCTION_TYPE_HPP_


/// possible types of post solving
enum class ReductionType : int
{
   kFixedCol = 0,
   kFixedInfCol = 5,
   kParallelCol = 2,
   kSubstitutedColWithDual = 3,
   kSubstitutedCol = 1,
   kVarBoundChange = 4,

   kRedundantRow = 7,
   kRowBoundChange = 8,
   kReasonForRowBoundChangeForcedByRow = 9,
   kRowBoundChangeForcedByRow = 10,

   kSaveRow = 11,

   kReducedBoundsCost = 12,
   kColumnDualValue = 13,
   kRowDualValue = 14,
   kCoefficientChange = 15,
};

#endif
