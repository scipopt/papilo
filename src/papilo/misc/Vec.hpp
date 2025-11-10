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

#ifndef _PAPILO_MISC_VEC_HPP_
#define _PAPILO_MISC_VEC_HPP_

#include "papilo/misc/Alloc.hpp"

// GCC 12 is confused about array access in small_vector and raises stringop warnings
// Boost >= 1.85 suppresses this warnings, but for earlier Boost, we do this here
#include <boost/version.hpp>
#include <boost/config.hpp>  // for BOOST_GCC
#if BOOST_VERSION < 108500 && defined(BOOST_GCC) && ((BOOST_GCC >= 120000) && (BOOST_GCC < 130000))
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wstringop-overread"
#  pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

#include <boost/container/small_vector.hpp>

#if BOOST_VERSION < 108500 && defined(BOOST_GCC) && ((BOOST_GCC >= 120000) && (BOOST_GCC < 130000))
#pragma GCC diagnostic pop
#endif

#include <vector>
namespace papilo
{
template <typename T>
using Vec = std::vector<T, Allocator<T>>;

template <typename T, int N>
using SmallVec = boost::container::small_vector<T, N, Allocator<T>>;

} // namespace papilo

#endif
