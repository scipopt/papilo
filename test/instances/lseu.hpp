/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2021  Konrad-Zuse-Zentrum                              */
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

#ifndef PAPILO_TEST_INSTANCES_LSEU
#define PAPILO_TEST_INSTANCES_LSEU

#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

namespace papilo
{
namespace instances
{

Problem<double>
lseu()
{
   /// PROBLEM BUILDER CODE
   Vec<double> coeffobj{
       7.0,   10.0,  179.0, 186.0, 179.0, 0.0,   0.0,   6.0,   0.0,   0.0,
       164.0, 164.0, 170.0, 164.0, 346.0, 346.0, 248.0, 253.0, 248.0, 346.0,
       346.0, 346.0, 346.0, 160.0, 161.0, 160.0, 160.0, 161.0, 160.0, 160.0,
       278.0, 278.0, 278.0, 86.0,  86.0,  86.0,  86.0,  86.0,  86.0,  188.0,
       188.0, 188.0, 85.0,  85.0,  85.0,  85.0,  85.0,  78.0,  78.0,  78.0,
       78.0,  78.0,  78.0,  78.0,  78.0,  78.0,  171.0, 171.0, 171.0, 171.0,
       163.0, 163.0, 163.0, 69.0,  69.0,  183.0, 183.0, 183.0, 183.0, 49.0,
       183.0, 258.0, 517.0, 250.0, 500.0, 250.0, 500.0, 159.0, 318.0, 159.0,
       318.0, 159.0, 318.0, 159.0, 318.0, 114.0, 228.0, 159.0, 318.0,
   };
   Vec<double> lbs{
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   };
   Vec<uint8_t> lbInf{
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   };
   Vec<double> ubs{
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   };
   Vec<uint8_t> ubInf{
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   };
   Vec<uint8_t> isIntegral{
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   };
   Vec<uint8_t> lhsIsInf{
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   };
   Vec<double> lhs{
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   };
   Vec<uint8_t> rhsIsInf{
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   };
   Vec<double> rhs{
       1.0,    1.0,     1.0,    1.0,     1.0,    1.0,     1.0,
       1.0,    1.0,     1.0,    1.0,     1.0,    1.0,     1.0,
       1.0,    1.0,     1.0,    -190.0,  2700.0, -2600.0, -630.0,
       -900.0, -1656.0, -335.0, -1026.0, -150.0, -500.0,  -270.0,
   };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 2, 1.0 },
       std::tuple<int, int, double>{ 0, 3, 1.0 },
       std::tuple<int, int, double>{ 0, 4, 1.0 },
       std::tuple<int, int, double>{ 1, 5, 1.0 },
       std::tuple<int, int, double>{ 1, 6, 1.0 },
       std::tuple<int, int, double>{ 1, 7, 1.0 },
       std::tuple<int, int, double>{ 1, 8, 1.0 },
       std::tuple<int, int, double>{ 1, 9, 1.0 },
       std::tuple<int, int, double>{ 2, 10, 1.0 },
       std::tuple<int, int, double>{ 2, 11, 1.0 },
       std::tuple<int, int, double>{ 2, 12, 1.0 },
       std::tuple<int, int, double>{ 2, 13, 1.0 },
       std::tuple<int, int, double>{ 3, 14, 1.0 },
       std::tuple<int, int, double>{ 3, 15, 1.0 },
       std::tuple<int, int, double>{ 4, 16, 1.0 },
       std::tuple<int, int, double>{ 4, 17, 1.0 },
       std::tuple<int, int, double>{ 4, 18, 1.0 },
       std::tuple<int, int, double>{ 5, 19, 1.0 },
       std::tuple<int, int, double>{ 5, 20, 1.0 },
       std::tuple<int, int, double>{ 5, 21, 1.0 },
       std::tuple<int, int, double>{ 5, 22, 1.0 },
       std::tuple<int, int, double>{ 6, 23, 1.0 },
       std::tuple<int, int, double>{ 6, 24, 1.0 },
       std::tuple<int, int, double>{ 6, 25, 1.0 },
       std::tuple<int, int, double>{ 6, 26, 1.0 },
       std::tuple<int, int, double>{ 6, 27, 1.0 },
       std::tuple<int, int, double>{ 6, 28, 1.0 },
       std::tuple<int, int, double>{ 6, 29, 1.0 },
       std::tuple<int, int, double>{ 7, 30, 1.0 },
       std::tuple<int, int, double>{ 7, 31, 1.0 },
       std::tuple<int, int, double>{ 7, 32, 1.0 },
       std::tuple<int, int, double>{ 8, 33, 1.0 },
       std::tuple<int, int, double>{ 8, 34, 1.0 },
       std::tuple<int, int, double>{ 8, 35, 1.0 },
       std::tuple<int, int, double>{ 8, 36, 1.0 },
       std::tuple<int, int, double>{ 8, 37, 1.0 },
       std::tuple<int, int, double>{ 9, 39, 1.0 },
       std::tuple<int, int, double>{ 9, 40, 1.0 },
       std::tuple<int, int, double>{ 9, 41, 1.0 },
       std::tuple<int, int, double>{ 10, 42, 1.0 },
       std::tuple<int, int, double>{ 10, 43, 1.0 },
       std::tuple<int, int, double>{ 10, 44, 1.0 },
       std::tuple<int, int, double>{ 10, 45, 1.0 },
       std::tuple<int, int, double>{ 10, 46, 1.0 },
       std::tuple<int, int, double>{ 11, 47, 1.0 },
       std::tuple<int, int, double>{ 11, 48, 1.0 },
       std::tuple<int, int, double>{ 11, 49, 1.0 },
       std::tuple<int, int, double>{ 11, 50, 1.0 },
       std::tuple<int, int, double>{ 12, 51, 1.0 },
       std::tuple<int, int, double>{ 12, 52, 1.0 },
       std::tuple<int, int, double>{ 12, 53, 1.0 },
       std::tuple<int, int, double>{ 12, 54, 1.0 },
       std::tuple<int, int, double>{ 12, 55, 1.0 },
       std::tuple<int, int, double>{ 13, 56, 1.0 },
       std::tuple<int, int, double>{ 13, 57, 1.0 },
       std::tuple<int, int, double>{ 13, 58, 1.0 },
       std::tuple<int, int, double>{ 13, 59, 1.0 },
       std::tuple<int, int, double>{ 14, 60, 1.0 },
       std::tuple<int, int, double>{ 14, 61, 1.0 },
       std::tuple<int, int, double>{ 14, 62, 1.0 },
       std::tuple<int, int, double>{ 15, 63, 1.0 },
       std::tuple<int, int, double>{ 15, 64, 1.0 },
       std::tuple<int, int, double>{ 16, 65, 1.0 },
       std::tuple<int, int, double>{ 16, 66, 1.0 },
       std::tuple<int, int, double>{ 16, 67, 1.0 },
       std::tuple<int, int, double>{ 16, 68, 1.0 },
       std::tuple<int, int, double>{ 16, 70, 1.0 },
       std::tuple<int, int, double>{ 17, 5, -450.0 },
       std::tuple<int, int, double>{ 17, 10, -435.0 },
       std::tuple<int, int, double>{ 17, 19, -435.0 },
       std::tuple<int, int, double>{ 17, 30, -350.0 },
       std::tuple<int, int, double>{ 17, 65, -230.0 },
       std::tuple<int, int, double>{ 17, 71, -200.0 },
       std::tuple<int, int, double>{ 17, 72, -400.0 },
       std::tuple<int, int, double>{ 18, 0, 525.0 },
       std::tuple<int, int, double>{ 18, 1, 500.0 },
       std::tuple<int, int, double>{ 18, 2, 475.0 },
       std::tuple<int, int, double>{ 18, 3, 475.0 },
       std::tuple<int, int, double>{ 18, 4, 475.0 },
       std::tuple<int, int, double>{ 18, 12, 435.0 },
       std::tuple<int, int, double>{ 18, 13, 435.0 },
       std::tuple<int, int, double>{ 18, 15, 435.0 },
       std::tuple<int, int, double>{ 18, 16, 435.0 },
       std::tuple<int, int, double>{ 18, 17, 435.0 },
       std::tuple<int, int, double>{ 18, 18, 435.0 },
       std::tuple<int, int, double>{ 18, 26, 425.0 },
       std::tuple<int, int, double>{ 18, 27, 425.0 },
       std::tuple<int, int, double>{ 18, 28, 425.0 },
       std::tuple<int, int, double>{ 18, 35, 330.0 },
       std::tuple<int, int, double>{ 18, 36, 330.0 },
       std::tuple<int, int, double>{ 18, 37, 330.0 },
       std::tuple<int, int, double>{ 18, 38, 330.0 },
       std::tuple<int, int, double>{ 18, 40, 330.0 },
       std::tuple<int, int, double>{ 18, 41, 330.0 },
       std::tuple<int, int, double>{ 18, 44, 325.0 },
       std::tuple<int, int, double>{ 18, 45, 325.0 },
       std::tuple<int, int, double>{ 18, 46, 325.0 },
       std::tuple<int, int, double>{ 18, 48, 300.0 },
       std::tuple<int, int, double>{ 18, 49, 300.0 },
       std::tuple<int, int, double>{ 18, 53, 300.0 },
       std::tuple<int, int, double>{ 18, 54, 300.0 },
       std::tuple<int, int, double>{ 18, 55, 300.0 },
       std::tuple<int, int, double>{ 18, 58, 300.0 },
       std::tuple<int, int, double>{ 18, 59, 300.0 },
       std::tuple<int, int, double>{ 18, 60, 285.0 },
       std::tuple<int, int, double>{ 18, 61, 285.0 },
       std::tuple<int, int, double>{ 18, 63, 265.0 },
       std::tuple<int, int, double>{ 18, 64, 265.0 },
       std::tuple<int, int, double>{ 18, 67, 230.0 },
       std::tuple<int, int, double>{ 18, 68, 230.0 },
       std::tuple<int, int, double>{ 18, 69, 190.0 },
       std::tuple<int, int, double>{ 18, 77, 200.0 },
       std::tuple<int, int, double>{ 18, 78, 400.0 },
       std::tuple<int, int, double>{ 18, 79, 200.0 },
       std::tuple<int, int, double>{ 18, 80, 400.0 },
       std::tuple<int, int, double>{ 18, 81, 200.0 },
       std::tuple<int, int, double>{ 18, 82, 400.0 },
       std::tuple<int, int, double>{ 18, 83, 200.0 },
       std::tuple<int, int, double>{ 18, 84, 400.0 },
       std::tuple<int, int, double>{ 18, 85, 200.0 },
       std::tuple<int, int, double>{ 18, 86, 400.0 },
       std::tuple<int, int, double>{ 19, 0, -525.0 },
       std::tuple<int, int, double>{ 19, 1, -500.0 },
       std::tuple<int, int, double>{ 19, 2, -475.0 },
       std::tuple<int, int, double>{ 19, 3, -475.0 },
       std::tuple<int, int, double>{ 19, 4, -475.0 },
       std::tuple<int, int, double>{ 19, 12, -435.0 },
       std::tuple<int, int, double>{ 19, 13, -435.0 },
       std::tuple<int, int, double>{ 19, 15, -435.0 },
       std::tuple<int, int, double>{ 19, 16, -435.0 },
       std::tuple<int, int, double>{ 19, 17, -435.0 },
       std::tuple<int, int, double>{ 19, 18, -435.0 },
       std::tuple<int, int, double>{ 19, 26, -425.0 },
       std::tuple<int, int, double>{ 19, 27, -425.0 },
       std::tuple<int, int, double>{ 19, 28, -425.0 },
       std::tuple<int, int, double>{ 19, 35, -330.0 },
       std::tuple<int, int, double>{ 19, 36, -330.0 },
       std::tuple<int, int, double>{ 19, 37, -330.0 },
       std::tuple<int, int, double>{ 19, 38, -330.0 },
       std::tuple<int, int, double>{ 19, 40, -330.0 },
       std::tuple<int, int, double>{ 19, 41, -330.0 },
       std::tuple<int, int, double>{ 19, 44, -325.0 },
       std::tuple<int, int, double>{ 19, 45, -325.0 },
       std::tuple<int, int, double>{ 19, 46, -325.0 },
       std::tuple<int, int, double>{ 19, 48, -300.0 },
       std::tuple<int, int, double>{ 19, 49, -300.0 },
       std::tuple<int, int, double>{ 19, 53, -300.0 },
       std::tuple<int, int, double>{ 19, 54, -300.0 },
       std::tuple<int, int, double>{ 19, 55, -300.0 },
       std::tuple<int, int, double>{ 19, 58, -300.0 },
       std::tuple<int, int, double>{ 19, 59, -300.0 },
       std::tuple<int, int, double>{ 19, 60, -285.0 },
       std::tuple<int, int, double>{ 19, 61, -285.0 },
       std::tuple<int, int, double>{ 19, 63, -265.0 },
       std::tuple<int, int, double>{ 19, 64, -265.0 },
       std::tuple<int, int, double>{ 19, 67, -230.0 },
       std::tuple<int, int, double>{ 19, 68, -230.0 },
       std::tuple<int, int, double>{ 19, 69, -190.0 },
       std::tuple<int, int, double>{ 19, 77, -200.0 },
       std::tuple<int, int, double>{ 19, 78, -400.0 },
       std::tuple<int, int, double>{ 19, 79, -200.0 },
       std::tuple<int, int, double>{ 19, 80, -400.0 },
       std::tuple<int, int, double>{ 19, 81, -200.0 },
       std::tuple<int, int, double>{ 19, 82, -400.0 },
       std::tuple<int, int, double>{ 19, 83, -200.0 },
       std::tuple<int, int, double>{ 19, 84, -400.0 },
       std::tuple<int, int, double>{ 19, 85, -200.0 },
       std::tuple<int, int, double>{ 19, 86, -400.0 },
       std::tuple<int, int, double>{ 20, 13, -435.0 },
       std::tuple<int, int, double>{ 20, 21, -400.0 },
       std::tuple<int, int, double>{ 20, 32, -350.0 },
       std::tuple<int, int, double>{ 20, 37, -330.0 },
       std::tuple<int, int, double>{ 20, 41, -330.0 },
       std::tuple<int, int, double>{ 20, 46, -325.0 },
       std::tuple<int, int, double>{ 20, 49, -300.0 },
       std::tuple<int, int, double>{ 20, 55, -300.0 },
       std::tuple<int, int, double>{ 20, 59, -300.0 },
       std::tuple<int, int, double>{ 20, 85, -200.0 },
       std::tuple<int, int, double>{ 20, 86, -400.0 },
       std::tuple<int, int, double>{ 21, 0, -525.0 },
       std::tuple<int, int, double>{ 21, 1, -500.0 },
       std::tuple<int, int, double>{ 21, 3, -475.0 },
       std::tuple<int, int, double>{ 21, 4, -190.0 },
       std::tuple<int, int, double>{ 21, 7, -450.0 },
       std::tuple<int, int, double>{ 21, 8, -165.0 },
       std::tuple<int, int, double>{ 21, 17, -435.0 },
       std::tuple<int, int, double>{ 21, 18, -300.0 },
       std::tuple<int, int, double>{ 21, 24, -400.0 },
       std::tuple<int, int, double>{ 21, 25, -115.0 },
       std::tuple<int, int, double>{ 21, 33, -330.0 },
       std::tuple<int, int, double>{ 21, 38, -330.0 },
       std::tuple<int, int, double>{ 21, 39, -330.0 },
       std::tuple<int, int, double>{ 21, 42, -325.0 },
       std::tuple<int, int, double>{ 21, 47, -300.0 },
       std::tuple<int, int, double>{ 21, 51, -300.0 },
       std::tuple<int, int, double>{ 21, 56, -300.0 },
       std::tuple<int, int, double>{ 21, 61, -285.0 },
       std::tuple<int, int, double>{ 21, 64, -265.0 },
       std::tuple<int, int, double>{ 21, 69, -190.0 },
       std::tuple<int, int, double>{ 21, 81, -200.0 },
       std::tuple<int, int, double>{ 21, 82, -400.0 },
       std::tuple<int, int, double>{ 22, 0, -525.0 },
       std::tuple<int, int, double>{ 22, 1, -500.0 },
       std::tuple<int, int, double>{ 22, 3, -475.0 },
       std::tuple<int, int, double>{ 22, 4, -190.0 },
       std::tuple<int, int, double>{ 22, 7, -450.0 },
       std::tuple<int, int, double>{ 22, 8, -165.0 },
       std::tuple<int, int, double>{ 22, 12, -435.0 },
       std::tuple<int, int, double>{ 22, 17, -435.0 },
       std::tuple<int, int, double>{ 22, 18, -300.0 },
       std::tuple<int, int, double>{ 22, 20, -400.0 },
       std::tuple<int, int, double>{ 22, 24, -400.0 },
       std::tuple<int, int, double>{ 22, 25, -115.0 },
       std::tuple<int, int, double>{ 22, 27, -425.0 },
       std::tuple<int, int, double>{ 22, 28, -140.0 },
       std::tuple<int, int, double>{ 22, 33, -330.0 },
       std::tuple<int, int, double>{ 22, 36, -330.0 },
       std::tuple<int, int, double>{ 22, 38, -330.0 },
       std::tuple<int, int, double>{ 22, 39, -330.0 },
       std::tuple<int, int, double>{ 22, 42, -325.0 },
       std::tuple<int, int, double>{ 22, 45, -325.0 },
       std::tuple<int, int, double>{ 22, 47, -300.0 },
       std::tuple<int, int, double>{ 22, 51, -300.0 },
       std::tuple<int, int, double>{ 22, 54, -300.0 },
       std::tuple<int, int, double>{ 22, 56, -300.0 },
       std::tuple<int, int, double>{ 22, 58, -300.0 },
       std::tuple<int, int, double>{ 22, 61, -285.0 },
       std::tuple<int, int, double>{ 22, 64, -265.0 },
       std::tuple<int, int, double>{ 22, 68, -230.0 },
       std::tuple<int, int, double>{ 22, 69, -190.0 },
       std::tuple<int, int, double>{ 22, 81, -200.0 },
       std::tuple<int, int, double>{ 22, 82, -400.0 },
       std::tuple<int, int, double>{ 22, 83, -200.0 },
       std::tuple<int, int, double>{ 22, 84, -400.0 },
       std::tuple<int, int, double>{ 23, 2, -475.0 },
       std::tuple<int, int, double>{ 23, 4, -285.0 },
       std::tuple<int, int, double>{ 23, 6, -450.0 },
       std::tuple<int, int, double>{ 23, 8, -285.0 },
       std::tuple<int, int, double>{ 23, 9, -150.0 },
       std::tuple<int, int, double>{ 23, 11, -435.0 },
       std::tuple<int, int, double>{ 23, 14, -435.0 },
       std::tuple<int, int, double>{ 23, 16, -435.0 },
       std::tuple<int, int, double>{ 23, 18, -135.0 },
       std::tuple<int, int, double>{ 23, 22, -100.0 },
       std::tuple<int, int, double>{ 23, 23, -400.0 },
       std::tuple<int, int, double>{ 23, 25, -285.0 },
       std::tuple<int, int, double>{ 23, 29, -100.0 },
       std::tuple<int, int, double>{ 23, 31, -350.0 },
       std::tuple<int, int, double>{ 23, 35, -330.0 },
       std::tuple<int, int, double>{ 23, 40, -330.0 },
       std::tuple<int, int, double>{ 23, 44, -325.0 },
       std::tuple<int, int, double>{ 23, 48, -300.0 },
       std::tuple<int, int, double>{ 23, 53, -300.0 },
       std::tuple<int, int, double>{ 23, 60, -285.0 },
       std::tuple<int, int, double>{ 23, 63, -265.0 },
       std::tuple<int, int, double>{ 23, 66, -230.0 },
       std::tuple<int, int, double>{ 23, 77, -200.0 },
       std::tuple<int, int, double>{ 23, 78, -400.0 },
       std::tuple<int, int, double>{ 24, 2, -475.0 },
       std::tuple<int, int, double>{ 24, 4, -285.0 },
       std::tuple<int, int, double>{ 24, 6, -450.0 },
       std::tuple<int, int, double>{ 24, 8, -285.0 },
       std::tuple<int, int, double>{ 24, 9, -150.0 },
       std::tuple<int, int, double>{ 24, 11, -435.0 },
       std::tuple<int, int, double>{ 24, 14, -435.0 },
       std::tuple<int, int, double>{ 24, 15, -435.0 },
       std::tuple<int, int, double>{ 24, 16, -435.0 },
       std::tuple<int, int, double>{ 24, 18, -135.0 },
       std::tuple<int, int, double>{ 24, 22, -100.0 },
       std::tuple<int, int, double>{ 24, 23, -400.0 },
       std::tuple<int, int, double>{ 24, 25, -285.0 },
       std::tuple<int, int, double>{ 24, 26, -425.0 },
       std::tuple<int, int, double>{ 24, 28, -285.0 },
       std::tuple<int, int, double>{ 24, 29, -100.0 },
       std::tuple<int, int, double>{ 24, 31, -350.0 },
       std::tuple<int, int, double>{ 24, 35, -330.0 },
       std::tuple<int, int, double>{ 24, 40, -330.0 },
       std::tuple<int, int, double>{ 24, 44, -325.0 },
       std::tuple<int, int, double>{ 24, 48, -300.0 },
       std::tuple<int, int, double>{ 24, 53, -300.0 },
       std::tuple<int, int, double>{ 24, 60, -285.0 },
       std::tuple<int, int, double>{ 24, 63, -265.0 },
       std::tuple<int, int, double>{ 24, 66, -230.0 },
       std::tuple<int, int, double>{ 24, 67, -230.0 },
       std::tuple<int, int, double>{ 24, 77, -200.0 },
       std::tuple<int, int, double>{ 24, 78, -400.0 },
       std::tuple<int, int, double>{ 24, 79, -200.0 },
       std::tuple<int, int, double>{ 24, 80, -400.0 },
       std::tuple<int, int, double>{ 25, 29, -300.0 },
       std::tuple<int, int, double>{ 25, 34, -330.0 },
       std::tuple<int, int, double>{ 25, 43, -325.0 },
       std::tuple<int, int, double>{ 25, 52, -300.0 },
       std::tuple<int, int, double>{ 25, 57, -300.0 },
       std::tuple<int, int, double>{ 25, 73, -200.0 },
       std::tuple<int, int, double>{ 25, 74, -400.0 },
       std::tuple<int, int, double>{ 26, 22, -300.0 },
       std::tuple<int, int, double>{ 26, 29, -300.0 },
       std::tuple<int, int, double>{ 26, 34, -330.0 },
       std::tuple<int, int, double>{ 26, 43, -325.0 },
       std::tuple<int, int, double>{ 26, 52, -300.0 },
       std::tuple<int, int, double>{ 26, 57, -300.0 },
       std::tuple<int, int, double>{ 26, 73, -200.0 },
       std::tuple<int, int, double>{ 26, 74, -400.0 },
       std::tuple<int, int, double>{ 26, 75, -200.0 },
       std::tuple<int, int, double>{ 26, 76, -400.0 },
       std::tuple<int, int, double>{ 27, 50, -300.0 },
       std::tuple<int, int, double>{ 27, 62, -285.0 },
       std::tuple<int, int, double>{ 27, 87, -200.0 },
       std::tuple<int, int, double>{ 27, 88, -400.0 },
   };
   Vec<std::string> rnames{
       "R101", "R102", "R103", "R104", "R105", "R106", "R107",
       "R108", "R109", "R110", "R111", "R112", "R113", "R114",
       "R115", "R116", "R117", "R118", "R119", "R120", "R121",
       "R122", "R123", "R124", "R125", "R126", "R127", "R128",
   };
   Vec<std::string> cnames{
       "C101", "C102", "C103", "C104", "C105", "C106", "C107", "C108", "C109",
       "C110", "C111", "C112", "C113", "C114", "C115", "C116", "C117", "C118",
       "C119", "C120", "C121", "C122", "C123", "C124", "C125", "C126", "C127",
       "C128", "C129", "C130", "C131", "C132", "C133", "C134", "C135", "C136",
       "C137", "C138", "C139", "C140", "C141", "C142", "C143", "C144", "C145",
       "C146", "C147", "C148", "C149", "C150", "C151", "C152", "C153", "C154",
       "C155", "C156", "C157", "C158", "C159", "C160", "C161", "C162", "C163",
       "C164", "C165", "C166", "C167", "C168", "C169", "C170", "C171", "C172",
       "C173", "C174", "C175", "C176", "C177", "C178", "C179", "C180", "C181",
       "C182", "C183", "C184", "C185", "C186", "C187", "C188", "C189",
   };
   int nCols = 89;
   int nRows = 28;
   ProblemBuilder<double> pb;
   pb.reserve( 309, 28, 89 );
   pb.setNumRows( nRows );
   pb.setNumCols( nCols );
   pb.setObjAll( coeffobj );
   pb.setObjOffset( 0.0 );
   pb.setColLbAll( lbs );
   pb.setColLbInfAll( lbInf );
   pb.setColUbAll( ubs );
   pb.setColUbInfAll( ubInf );
   pb.setColIntegralAll( isIntegral );
   pb.setRowLhsInfAll( lhsIsInf );
   pb.setRowRhsInfAll( rhsIsInf );
   pb.setRowLhsAll( lhs );
   pb.setRowRhsAll( rhs );
   pb.setRowNameAll( rnames );
   pb.addEntryAll( entries );
   pb.setColNameAll( cnames );
   pb.setProblemName( "lseu.hpp" );
   Problem<double> problem = pb.build();
   /// PROBLEM BUILDER CODE END

   return problem;
}

} // namespace instances
} // namespace papilo

#endif
