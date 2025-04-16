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

#ifndef _PAPILO_MISC_TIMER_HPP_
#define _PAPILO_MISC_TIMER_HPP_
#include "papilo/Config.hpp"

#ifdef PAPILO_TBB
#include "papilo/misc/tbb.hpp"
#else
#include <chrono>
#endif

namespace papilo
{

#ifdef PAPILO_TBB
class Timer
{
 public:
   Timer( double& time_ ) : time( time_ ) { start = tbb::tick_count::now(); }

   double
   getTime() const
   {
      return ( tbb::tick_count::now() - start ).seconds();
   }

   ~Timer() { time += ( tbb::tick_count::now() - start ).seconds(); }

 private:
   tbb::tick_count start;
   double& time;
};
#else
class Timer
{
 public:
   Timer( double& time_ ) : time( time_ ) { start = std::chrono::steady_clock::now(); }

   double
   getTime() const
   {
      return std::chrono::duration_cast<std::chrono::milliseconds>(
                 std::chrono::steady_clock::now() - start )
                 .count() /1000.0;
   }

   ~Timer() {
      time += std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - start )
                  .count() /1000.0;
   }

 private:
   std::chrono::steady_clock::time_point start;
   double& time;
};
#endif

} // namespace papilo

#endif
