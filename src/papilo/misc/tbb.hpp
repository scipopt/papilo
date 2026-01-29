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

#ifndef _PAPILO_MISC_TBB_HPP_
#define _PAPILO_MISC_TBB_HPP_

/* if those macros are not defined and tbb includes windows.h
 * then many macros are defined that can interfere with standard C++ code
 */
#ifndef NOMINMAX
#define NOMINMAX
#define PAPILO_DEFINED_NOMINMAX
#endif

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#define PAPILO_DEFINED_WIN32_LEAN_AND_MEAN
#endif

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#define PAPILO_DEFINED_WIN32_LEAN_AND_MEAN
#endif

#ifndef NOGDI
#define NOGDI
#define PAPILO_DEFINED_NOGDI
#endif

#ifdef _MSC_VER
#pragma push_macro( "__TBB_NO_IMPLICIT_LINKAGE" )
#define __TBB_NO_IMPLICIT_LINKAGE 1
#endif

#include "tbb/blocked_range.h"
#include "tbb/combinable.h"
#include "tbb/concurrent_vector.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"
#include "tbb/partitioner.h"
#include "tbb/task_arena.h"
#include "tbb/tick_count.h"

#ifdef _MSC_VER
#pragma pop_macro( "__TBB_NO_IMPLICIT_LINKAGE" )
#endif

#ifdef PAPILO_DEFINED_NOGDI
#undef NOGDI
#undef PAPILO_DEFINED_NOGDI
#endif

#ifdef PAPILO_DEFINED_NOMINMAX
#undef NOMINMAX
#undef PAPILO_DEFINED_NOMINMAX
#endif

#ifdef PAPILO_DEFINED_WIN32_LEAN_AND_MEAN
#undef WIN32_LEAN_AND_MEAN
#undef PAPILO_DEFINED_WIN32_LEAN_AND_MEAN
#endif

#endif
