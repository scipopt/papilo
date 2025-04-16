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

#ifndef _PAPILO_MISC_MULTIPRECISION_HPP_
#define _PAPILO_MISC_MULTIPRECISION_HPP_

#include "papilo/Config.hpp"

// work around build failure with boost on Fedora 37
#include <memory>

#ifdef PAPILO_SERIALIZATION_AVAILABLE
#include <boost/serialization/split_free.hpp>
#endif

#ifdef PAPILO_HAVE_FLOAT128
#include <boost/multiprecision/float128.hpp>
namespace papilo
{
using Quad = boost::multiprecision::float128;
} // namespace papilo
#elif defined( PAPILO_HAVE_GMP )
#include <boost/multiprecision/gmp.hpp>

namespace papilo
{

using Quad =
    boost::multiprecision::number<boost::multiprecision::gmp_float<35>>;
} // namespace papilo

#ifdef PAPILO_SERIALIZATION_AVAILABLE
BOOST_SERIALIZATION_SPLIT_FREE( papilo::Quad )
#endif

#else
#include <boost/multiprecision/cpp_bin_float.hpp>
#ifdef PAPILO_SERIALIZATION_AVAILABLE
#include <boost/serialization/nvp.hpp>
#endif
namespace papilo
{
using Quad = boost::multiprecision::cpp_bin_float_quad;
} // namespace papilo
#endif

#ifdef PAPILO_HAVE_GMP
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/gmp.hpp>
#ifdef PAPILO_SERIALIZATION_AVAILABLE
#include <boost/serialization/nvp.hpp>
#endif

// unfortunately the multiprecision gmp types do not provide an overload for
// serialization
namespace papilo
{
using Integral = boost::multiprecision::mpz_int;
using Rational = boost::multiprecision::mpq_rational;
using Float100 = boost::multiprecision::mpf_float_100;
using Float500 = boost::multiprecision::mpf_float_500;
using Float1000 = boost::multiprecision::mpf_float_1000;
} // namespace papilo

#ifdef PAPILO_SERIALIZATION_AVAILABLE
BOOST_SERIALIZATION_SPLIT_FREE( papilo::Rational )
BOOST_SERIALIZATION_SPLIT_FREE( papilo::Float100 )
BOOST_SERIALIZATION_SPLIT_FREE( papilo::Float500 )
BOOST_SERIALIZATION_SPLIT_FREE( papilo::Float1000 )

namespace boost
{
namespace serialization
{

template <class Archive>
void
save( Archive& ar, const papilo::Rational& num, const unsigned int version )
{
   boost::multiprecision::cpp_rational t( num );
   ar& t;
}

template <class Archive>
void
load( Archive& ar, papilo::Rational& num, const unsigned int version )
{
   boost::multiprecision::cpp_rational t;
   ar& t;
   num = papilo::Rational( t );
}

template <class Archive, unsigned M>
void
save( Archive& ar,
      const boost::multiprecision::number<boost::multiprecision::gmp_float<M>>&
          num,
      const unsigned int version )
{
   boost::multiprecision::number<boost::multiprecision::cpp_bin_float<M>> t(
       num );
   ar& t;
}

template <class Archive, unsigned M>
void
load( Archive& ar,
      boost::multiprecision::number<boost::multiprecision::gmp_float<M>>& num,
      const unsigned int version )
{
   boost::multiprecision::number<boost::multiprecision::cpp_bin_float<M>> t;
   ar& t;
   num =
       boost::multiprecision::number<boost::multiprecision::gmp_float<M>>( t );
}

} // namespace serialization
} // namespace boost
#endif

#else
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#ifdef PAPILO_SERIALIZATION_AVAILABLE
#include <boost/serialization/nvp.hpp>
#endif
namespace papilo
{
using Integral = boost::multiprecision::cpp_int;
using Rational = boost::multiprecision::cpp_rational;
using Float100 =
    boost::multiprecision::number<boost::multiprecision::cpp_bin_float<100>>;
using Float500 =
    boost::multiprecision::number<boost::multiprecision::cpp_bin_float<500>>;
using Float1000 =
    boost::multiprecision::number<boost::multiprecision::cpp_bin_float<1000>>;
} // namespace papilo
#endif

#endif
