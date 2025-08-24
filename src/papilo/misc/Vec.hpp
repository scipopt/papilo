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
