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

#ifndef PAPILO_TEST_INSTANCES_BELL5
#define PAPILO_TEST_INSTANCES_BELL5

#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

namespace papilo
{
namespace instances
{

Problem<double>
bell5()
{
   /// PROBLEM BUILDER CODE
   Vec<double> coeffobj{
       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
       43000.0, 43000.0, 43000.0, 43000.0, 43000.0, 43000.0, 43000.0, 43000.0,
       43000.0, 43000.0, 43000.0, 43000.0, 43000.0, 43000.0, 58000.0, 58000.0,
       58000.0, 59000.0, 59000.0, 59000.0, 59000.0, 60000.0, 59000.0, 59000.0,
       59000.0, 58000.0, 58000.0, 58000.0, 10000.0, 10000.0, 10000.0, 10000.0,
       10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0,
       10000.0, 10000.0, 24.5645, 20.3962, 14.1693, 50.2605, 58.0423, 36.6095,
       39.201,  48.034,  29.4336, 36.0182, 18.7245, 30.3169, 5.3655,  25.55,
       20.7977, 1.825,   2.45645, 2.03962, 1.41693, 5.02605, 5.80423, 3.66095,
       3.9201,  4.8034,  2.94336, 3.60182, 1.87245, 3.03169, 0.53655, 2.555,
       2.07977, 0.1825,  0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
   };
   Vec<double> lbs{
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   };
   Vec<uint8_t> lbInf{
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   };
   Vec<double> ubs{
       1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,
       1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,
       1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,
       1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     10000.0, 10000.0,
       10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 1000.0,  1000.0,
       1000.0,  1000.0,  1000.0,  1000.0,  1000.0,  1000.0,  1000.0,  1000.0,
       1000.0,  1000.0,  1000.0,  1000.0,  100.0,   100.0,   100.0,   100.0,
       100.0,   100.0,   0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
   };
   Vec<uint8_t> ubInf{
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   };
   Vec<uint8_t> isIntegral{
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   };
   Vec<uint8_t> lhsIsInf{
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   };
   Vec<double> lhs{
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   };
   Vec<uint8_t> rhsIsInf{
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   };
   Vec<double> rhs{
       0.0,   0.0,     0.0,     0.0,     0.0,     0.0,    0.0,     0.0,
       0.0,   0.0,     0.0,     0.0,     0.0,     0.0,    0.0,     0.0,
       0.0,   0.0,     0.0,     0.0,     0.0,     0.0,    0.0,     0.0,
       0.0,   0.0,     0.0,     0.0,     0.0,     0.0,    0.0,     0.0,
       0.0,   0.0,     0.0,     0.0,     0.0,     0.0,    0.0,     0.0,
       0.0,   0.0,     0.0,     0.0,     0.0,     0.0,    0.0,     0.0,
       0.0,   0.0,     0.0,     0.0,     0.0,     0.0,    0.0,     0.0,
       0.0,   0.0,     0.0,     -800.0,  -3100.0, -50.0,  -40.0,   -900.0,
       -40.0, -7150.0, -4800.0, -3000.0, -250.0,  -400.0, -1000.0, -4000.0,
       -90.0, -700.0,  -24.0,   8.0,     8.0,     8.0,    8.0,     1.0,
       1.0,   13.0,    13.0,    1.0,     1.0,     13.0,   13.0,    8.0,
       8.0,   2.0,     2.0,
   };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, -1.0 },
       std::tuple<int, int, double>{ 0, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 1, -1.0 },
       std::tuple<int, int, double>{ 1, 2, 1.0 },
       std::tuple<int, int, double>{ 2, 2, -1.0 },
       std::tuple<int, int, double>{ 2, 3, 1.0 },
       std::tuple<int, int, double>{ 3, 3, -1.0 },
       std::tuple<int, int, double>{ 3, 4, 1.0 },
       std::tuple<int, int, double>{ 4, 4, -1.0 },
       std::tuple<int, int, double>{ 4, 5, 1.0 },
       std::tuple<int, int, double>{ 5, 4, -1.0 },
       std::tuple<int, int, double>{ 5, 6, 1.0 },
       std::tuple<int, int, double>{ 6, 6, -1.0 },
       std::tuple<int, int, double>{ 6, 7, 1.0 },
       std::tuple<int, int, double>{ 7, 7, -1.0 },
       std::tuple<int, int, double>{ 7, 8, 1.0 },
       std::tuple<int, int, double>{ 8, 3, -1.0 },
       std::tuple<int, int, double>{ 8, 9, 1.0 },
       std::tuple<int, int, double>{ 9, 9, -1.0 },
       std::tuple<int, int, double>{ 9, 10, 1.0 },
       std::tuple<int, int, double>{ 10, 10, -1.0 },
       std::tuple<int, int, double>{ 10, 11, 1.0 },
       std::tuple<int, int, double>{ 11, 1, -1.0 },
       std::tuple<int, int, double>{ 11, 12, 1.0 },
       std::tuple<int, int, double>{ 12, 12, -1.0 },
       std::tuple<int, int, double>{ 12, 13, 1.0 },
       std::tuple<int, int, double>{ 13, 13, -1.0 },
       std::tuple<int, int, double>{ 13, 14, 1.0 },
       std::tuple<int, int, double>{ 14, 0, -1.0 },
       std::tuple<int, int, double>{ 14, 15, 1.0 },
       std::tuple<int, int, double>{ 15, 0, -20.0 },
       std::tuple<int, int, double>{ 15, 16, 1.0 },
       std::tuple<int, int, double>{ 15, 30, 1.0 },
       std::tuple<int, int, double>{ 16, 1, -20.0 },
       std::tuple<int, int, double>{ 16, 17, 1.0 },
       std::tuple<int, int, double>{ 16, 31, 1.0 },
       std::tuple<int, int, double>{ 17, 2, -20.0 },
       std::tuple<int, int, double>{ 17, 18, 1.0 },
       std::tuple<int, int, double>{ 17, 32, 1.0 },
       std::tuple<int, int, double>{ 18, 3, -20.0 },
       std::tuple<int, int, double>{ 18, 19, 1.0 },
       std::tuple<int, int, double>{ 18, 33, 1.0 },
       std::tuple<int, int, double>{ 19, 4, -20.0 },
       std::tuple<int, int, double>{ 19, 20, 1.0 },
       std::tuple<int, int, double>{ 19, 34, 1.0 },
       std::tuple<int, int, double>{ 20, 5, -20.0 },
       std::tuple<int, int, double>{ 20, 21, 1.0 },
       std::tuple<int, int, double>{ 20, 35, 1.0 },
       std::tuple<int, int, double>{ 21, 6, -20.0 },
       std::tuple<int, int, double>{ 21, 22, 1.0 },
       std::tuple<int, int, double>{ 21, 36, 1.0 },
       std::tuple<int, int, double>{ 22, 8, -20.0 },
       std::tuple<int, int, double>{ 22, 23, 1.0 },
       std::tuple<int, int, double>{ 22, 37, 1.0 },
       std::tuple<int, int, double>{ 23, 9, -20.0 },
       std::tuple<int, int, double>{ 23, 24, 1.0 },
       std::tuple<int, int, double>{ 23, 38, 1.0 },
       std::tuple<int, int, double>{ 24, 10, -20.0 },
       std::tuple<int, int, double>{ 24, 25, 1.0 },
       std::tuple<int, int, double>{ 24, 39, 1.0 },
       std::tuple<int, int, double>{ 25, 11, -20.0 },
       std::tuple<int, int, double>{ 25, 26, 1.0 },
       std::tuple<int, int, double>{ 25, 40, 1.0 },
       std::tuple<int, int, double>{ 26, 12, -20.0 },
       std::tuple<int, int, double>{ 26, 27, 1.0 },
       std::tuple<int, int, double>{ 26, 41, 1.0 },
       std::tuple<int, int, double>{ 27, 13, -20.0 },
       std::tuple<int, int, double>{ 27, 28, 1.0 },
       std::tuple<int, int, double>{ 27, 42, 1.0 },
       std::tuple<int, int, double>{ 28, 15, -20.0 },
       std::tuple<int, int, double>{ 28, 29, 1.0 },
       std::tuple<int, int, double>{ 28, 43, 1.0 },
       std::tuple<int, int, double>{ 29, 16, -672.0 },
       std::tuple<int, int, double>{ 29, 30, -1344.0 },
       std::tuple<int, int, double>{ 29, 74, -1.0 },
       std::tuple<int, int, double>{ 29, 75, 1.0 },
       std::tuple<int, int, double>{ 29, 89, 1.0 },
       std::tuple<int, int, double>{ 29, 90, 1.0 },
       std::tuple<int, int, double>{ 30, 17, -672.0 },
       std::tuple<int, int, double>{ 30, 31, -1344.0 },
       std::tuple<int, int, double>{ 30, 75, -1.0 },
       std::tuple<int, int, double>{ 30, 76, 1.0 },
       std::tuple<int, int, double>{ 30, 86, 1.0 },
       std::tuple<int, int, double>{ 30, 91, 1.0 },
       std::tuple<int, int, double>{ 31, 18, -672.0 },
       std::tuple<int, int, double>{ 31, 32, -1344.0 },
       std::tuple<int, int, double>{ 31, 76, -1.0 },
       std::tuple<int, int, double>{ 31, 77, 1.0 },
       std::tuple<int, int, double>{ 31, 92, 1.0 },
       std::tuple<int, int, double>{ 32, 19, -672.0 },
       std::tuple<int, int, double>{ 32, 33, -1344.0 },
       std::tuple<int, int, double>{ 32, 77, -1.0 },
       std::tuple<int, int, double>{ 32, 78, 1.0 },
       std::tuple<int, int, double>{ 32, 83, 1.0 },
       std::tuple<int, int, double>{ 32, 93, 1.0 },
       std::tuple<int, int, double>{ 33, 20, -672.0 },
       std::tuple<int, int, double>{ 33, 34, -1344.0 },
       std::tuple<int, int, double>{ 33, 78, -1.0 },
       std::tuple<int, int, double>{ 33, 79, 1.0 },
       std::tuple<int, int, double>{ 33, 80, 1.0 },
       std::tuple<int, int, double>{ 33, 94, 1.0 },
       std::tuple<int, int, double>{ 34, 21, -672.0 },
       std::tuple<int, int, double>{ 34, 35, -1344.0 },
       std::tuple<int, int, double>{ 34, 79, -1.0 },
       std::tuple<int, int, double>{ 34, 95, 1.0 },
       std::tuple<int, int, double>{ 35, 22, -672.0 },
       std::tuple<int, int, double>{ 35, 36, -1344.0 },
       std::tuple<int, int, double>{ 35, 80, -1.0 },
       std::tuple<int, int, double>{ 35, 81, 1.0 },
       std::tuple<int, int, double>{ 35, 96, 1.0 },
       std::tuple<int, int, double>{ 36, 81, -1.0 },
       std::tuple<int, int, double>{ 36, 82, 1.0 },
       std::tuple<int, int, double>{ 37, 23, -672.0 },
       std::tuple<int, int, double>{ 37, 37, -1344.0 },
       std::tuple<int, int, double>{ 37, 82, -1.0 },
       std::tuple<int, int, double>{ 37, 97, 1.0 },
       std::tuple<int, int, double>{ 38, 24, -672.0 },
       std::tuple<int, int, double>{ 38, 38, -1344.0 },
       std::tuple<int, int, double>{ 38, 83, -1.0 },
       std::tuple<int, int, double>{ 38, 84, 1.0 },
       std::tuple<int, int, double>{ 38, 98, 1.0 },
       std::tuple<int, int, double>{ 39, 25, -672.0 },
       std::tuple<int, int, double>{ 39, 39, -1344.0 },
       std::tuple<int, int, double>{ 39, 84, -1.0 },
       std::tuple<int, int, double>{ 39, 85, 1.0 },
       std::tuple<int, int, double>{ 39, 99, 1.0 },
       std::tuple<int, int, double>{ 40, 26, -672.0 },
       std::tuple<int, int, double>{ 40, 40, -1344.0 },
       std::tuple<int, int, double>{ 40, 85, -1.0 },
       std::tuple<int, int, double>{ 40, 100, 1.0 },
       std::tuple<int, int, double>{ 41, 27, -672.0 },
       std::tuple<int, int, double>{ 41, 41, -1344.0 },
       std::tuple<int, int, double>{ 41, 86, -1.0 },
       std::tuple<int, int, double>{ 41, 87, 1.0 },
       std::tuple<int, int, double>{ 41, 101, 1.0 },
       std::tuple<int, int, double>{ 42, 28, -672.0 },
       std::tuple<int, int, double>{ 42, 42, -1344.0 },
       std::tuple<int, int, double>{ 42, 87, -1.0 },
       std::tuple<int, int, double>{ 42, 88, 1.0 },
       std::tuple<int, int, double>{ 42, 102, 1.0 },
       std::tuple<int, int, double>{ 43, 88, -1.0 },
       std::tuple<int, int, double>{ 44, 29, -672.0 },
       std::tuple<int, int, double>{ 44, 43, -1344.0 },
       std::tuple<int, int, double>{ 44, 89, -1.0 },
       std::tuple<int, int, double>{ 44, 103, 1.0 },
       std::tuple<int, int, double>{ 45, 44, -24.0 },
       std::tuple<int, int, double>{ 45, 90, 1.0 },
       std::tuple<int, int, double>{ 46, 45, -24.0 },
       std::tuple<int, int, double>{ 46, 91, 1.0 },
       std::tuple<int, int, double>{ 47, 46, -24.0 },
       std::tuple<int, int, double>{ 47, 92, 1.0 },
       std::tuple<int, int, double>{ 48, 47, -24.0 },
       std::tuple<int, int, double>{ 48, 93, 1.0 },
       std::tuple<int, int, double>{ 49, 48, -24.0 },
       std::tuple<int, int, double>{ 49, 94, 1.0 },
       std::tuple<int, int, double>{ 50, 49, -24.0 },
       std::tuple<int, int, double>{ 50, 95, 1.0 },
       std::tuple<int, int, double>{ 51, 50, -24.0 },
       std::tuple<int, int, double>{ 51, 96, 1.0 },
       std::tuple<int, int, double>{ 52, 51, -24.0 },
       std::tuple<int, int, double>{ 52, 97, 1.0 },
       std::tuple<int, int, double>{ 53, 52, -24.0 },
       std::tuple<int, int, double>{ 53, 98, 1.0 },
       std::tuple<int, int, double>{ 54, 53, -24.0 },
       std::tuple<int, int, double>{ 54, 99, 1.0 },
       std::tuple<int, int, double>{ 55, 54, -24.0 },
       std::tuple<int, int, double>{ 55, 100, 1.0 },
       std::tuple<int, int, double>{ 56, 55, -24.0 },
       std::tuple<int, int, double>{ 56, 101, 1.0 },
       std::tuple<int, int, double>{ 57, 56, -24.0 },
       std::tuple<int, int, double>{ 57, 102, 1.0 },
       std::tuple<int, int, double>{ 58, 57, -24.0 },
       std::tuple<int, int, double>{ 58, 103, 1.0 },
       std::tuple<int, int, double>{ 59, 58, -1.0 },
       std::tuple<int, int, double>{ 59, 59, 1.0 },
       std::tuple<int, int, double>{ 59, 73, 1.0 },
       std::tuple<int, int, double>{ 59, 90, -1.0 },
       std::tuple<int, int, double>{ 60, 59, -1.0 },
       std::tuple<int, int, double>{ 60, 60, 1.0 },
       std::tuple<int, int, double>{ 60, 70, 1.0 },
       std::tuple<int, int, double>{ 60, 91, -1.0 },
       std::tuple<int, int, double>{ 61, 60, -1.0 },
       std::tuple<int, int, double>{ 61, 61, 1.0 },
       std::tuple<int, int, double>{ 61, 92, -1.0 },
       std::tuple<int, int, double>{ 62, 61, -1.0 },
       std::tuple<int, int, double>{ 62, 62, 1.0 },
       std::tuple<int, int, double>{ 62, 67, 1.0 },
       std::tuple<int, int, double>{ 62, 93, -1.0 },
       std::tuple<int, int, double>{ 63, 62, -1.0 },
       std::tuple<int, int, double>{ 63, 63, 1.0 },
       std::tuple<int, int, double>{ 63, 64, 1.0 },
       std::tuple<int, int, double>{ 63, 94, -1.0 },
       std::tuple<int, int, double>{ 64, 63, -1.0 },
       std::tuple<int, int, double>{ 64, 95, -1.0 },
       std::tuple<int, int, double>{ 65, 64, -1.0 },
       std::tuple<int, int, double>{ 65, 65, 1.0 },
       std::tuple<int, int, double>{ 65, 96, -1.0 },
       std::tuple<int, int, double>{ 66, 65, -1.0 },
       std::tuple<int, int, double>{ 66, 66, 1.0 },
       std::tuple<int, int, double>{ 67, 66, -1.0 },
       std::tuple<int, int, double>{ 67, 97, -1.0 },
       std::tuple<int, int, double>{ 68, 67, -1.0 },
       std::tuple<int, int, double>{ 68, 68, 1.0 },
       std::tuple<int, int, double>{ 68, 98, -1.0 },
       std::tuple<int, int, double>{ 69, 68, -1.0 },
       std::tuple<int, int, double>{ 69, 69, 1.0 },
       std::tuple<int, int, double>{ 69, 99, -1.0 },
       std::tuple<int, int, double>{ 70, 69, -1.0 },
       std::tuple<int, int, double>{ 70, 100, -1.0 },
       std::tuple<int, int, double>{ 71, 70, -1.0 },
       std::tuple<int, int, double>{ 71, 71, 1.0 },
       std::tuple<int, int, double>{ 71, 101, -1.0 },
       std::tuple<int, int, double>{ 72, 71, -1.0 },
       std::tuple<int, int, double>{ 72, 72, 1.0 },
       std::tuple<int, int, double>{ 72, 102, -1.0 },
       std::tuple<int, int, double>{ 73, 72, -1.0 },
       std::tuple<int, int, double>{ 74, 73, -1.0 },
       std::tuple<int, int, double>{ 74, 103, -1.0 },
       std::tuple<int, int, double>{ 75, 0, 1.0 },
       std::tuple<int, int, double>{ 75, 58, 0.000833 },
       std::tuple<int, int, double>{ 75, 74, 8.3e-05 },
       std::tuple<int, int, double>{ 76, 1, 1.0 },
       std::tuple<int, int, double>{ 76, 59, 0.000833 },
       std::tuple<int, int, double>{ 76, 75, 8.3e-05 },
       std::tuple<int, int, double>{ 77, 2, 1.0 },
       std::tuple<int, int, double>{ 77, 60, 0.000833 },
       std::tuple<int, int, double>{ 77, 76, 8.3e-05 },
       std::tuple<int, int, double>{ 78, 3, 1.0 },
       std::tuple<int, int, double>{ 78, 61, 0.000833 },
       std::tuple<int, int, double>{ 78, 77, 8.3e-05 },
       std::tuple<int, int, double>{ 79, 4, 1.0 },
       std::tuple<int, int, double>{ 79, 62, 0.000833 },
       std::tuple<int, int, double>{ 79, 78, 8.3e-05 },
       std::tuple<int, int, double>{ 80, 5, 1.0 },
       std::tuple<int, int, double>{ 80, 63, 0.000833 },
       std::tuple<int, int, double>{ 80, 79, 8.3e-05 },
       std::tuple<int, int, double>{ 81, 6, 1.0 },
       std::tuple<int, int, double>{ 81, 64, 0.000833 },
       std::tuple<int, int, double>{ 81, 80, 8.3e-05 },
       std::tuple<int, int, double>{ 82, 7, 1.0 },
       std::tuple<int, int, double>{ 82, 65, 0.000833 },
       std::tuple<int, int, double>{ 82, 81, 8.3e-05 },
       std::tuple<int, int, double>{ 83, 8, 1.0 },
       std::tuple<int, int, double>{ 83, 66, 0.000833 },
       std::tuple<int, int, double>{ 83, 82, 8.3e-05 },
       std::tuple<int, int, double>{ 84, 9, 1.0 },
       std::tuple<int, int, double>{ 84, 67, 0.000833 },
       std::tuple<int, int, double>{ 84, 83, 8.3e-05 },
       std::tuple<int, int, double>{ 85, 10, 1.0 },
       std::tuple<int, int, double>{ 85, 68, 0.000833 },
       std::tuple<int, int, double>{ 85, 84, 8.3e-05 },
       std::tuple<int, int, double>{ 86, 11, 1.0 },
       std::tuple<int, int, double>{ 86, 69, 0.000833 },
       std::tuple<int, int, double>{ 86, 85, 8.3e-05 },
       std::tuple<int, int, double>{ 87, 12, 1.0 },
       std::tuple<int, int, double>{ 87, 70, 0.000833 },
       std::tuple<int, int, double>{ 87, 86, 8.3e-05 },
       std::tuple<int, int, double>{ 88, 13, 1.0 },
       std::tuple<int, int, double>{ 88, 71, 0.000833 },
       std::tuple<int, int, double>{ 88, 87, 8.3e-05 },
       std::tuple<int, int, double>{ 89, 14, 1.0 },
       std::tuple<int, int, double>{ 89, 72, 0.000833 },
       std::tuple<int, int, double>{ 89, 88, 8.3e-05 },
       std::tuple<int, int, double>{ 90, 15, 1.0 },
       std::tuple<int, int, double>{ 90, 73, 0.000833 },
       std::tuple<int, int, double>{ 90, 89, 8.3e-05 },
   };
   Vec<std::string> rnames{
       "A1",  "A2",  "A3",  "A4",  "A5",  "A6",  "A7",  "A8",  "A9",  "A10",
       "A11", "A12", "A13", "A14", "A15", "B1",  "B2",  "B3",  "B4",  "B5",
       "B6",  "B7",  "B9",  "B10", "B11", "B12", "B13", "B14", "B16", "C1",
       "C2",  "C3",  "C4",  "C5",  "C6",  "C7",  "C8",  "C9",  "C10", "C11",
       "C12", "C13", "C14", "C15", "C16", "D1",  "D2",  "D3",  "D4",  "D5",
       "D6",  "D7",  "D9",  "D10", "D11", "D12", "D13", "D14", "D16", "E1",
       "E2",  "E3",  "E4",  "E5",  "E6",  "E7",  "E8",  "E9",  "E10", "E11",
       "E12", "E13", "E14", "E15", "E16", "F0",  "F1",  "F2",  "F3",  "F4",
       "F5",  "F6",  "F7",  "F8",  "F9",  "F10", "F11", "F12", "F13", "F14",
       "F15",
   };
   Vec<std::string> cnames{
       "c1",  "c2",  "c3",  "c4",  "c5",  "c6",  "c7",  "c8",  "c9",  "c10",
       "c11", "c12", "c13", "c14", "c15", "c16", "d1",  "d2",  "d3",  "d4",
       "d5",  "d6",  "d7",  "d9",  "d10", "d11", "d12", "d13", "d14", "d16",
       "h1",  "h2",  "h3",  "h4",  "h5",  "h6",  "h7",  "h9",  "h10", "h11",
       "h12", "h13", "h14", "h16", "g1",  "g2",  "g3",  "g4",  "g5",  "g6",
       "g7",  "g9",  "g10", "g11", "g12", "g13", "g14", "g16", "a1",  "a2",
       "a3",  "a4",  "a5",  "a6",  "a7",  "a8",  "a9",  "a10", "a11", "a12",
       "a13", "a14", "a15", "a16", "b1",  "b2",  "b3",  "b4",  "b5",  "b6",
       "b7",  "b8",  "b9",  "b10", "b11", "b12", "b13", "b14", "b15", "b16",
       "f1",  "f2",  "f3",  "f4",  "f5",  "f6",  "f7",  "f9",  "f10", "f11",
       "f12", "f13", "f14", "f16",
   };
   int nCols = 104;
   int nRows = 91;
   ProblemBuilder<double> pb;
   pb.reserve( 266, 91, 104 );
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
   pb.setProblemName( "bell5.hpp" );
   Problem<double> problem = pb.build();
   /// PROBLEM BUILDER CODE END

   return problem;
}

} // namespace instances
} // namespace papilo

#endif
