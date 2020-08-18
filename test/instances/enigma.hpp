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

#ifndef PAPILO_TEST_INSTANCES_ENIGMA
#define PAPILO_TEST_INSTANCES_ENIGMA

#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemBuilder.hpp"

namespace papilo
{
namespace instances
{

Problem<double>
enigma()
{
   /// PROBLEM BUILDER CODE
   Vec<double> coeffobj{
       0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   };
   Vec<double> lbs{
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   };
   Vec<uint8_t> lbInf{
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   };
   Vec<double> ubs{
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   };
   Vec<uint8_t> ubInf{
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   };
   Vec<uint8_t> isIntegral{
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   };
   Vec<uint8_t> lhsIsInf{
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   };
   Vec<double> lhs{
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   };
   Vec<uint8_t> rhsIsInf{
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   };
   Vec<double> rhs{
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   };
   Vec<std::tuple<int, int, double>> entries{
       std::tuple<int, int, double>{ 0, 0, 1.0 },
       std::tuple<int, int, double>{ 0, 10, 1.0 },
       std::tuple<int, int, double>{ 0, 20, 1.0 },
       std::tuple<int, int, double>{ 0, 30, 1.0 },
       std::tuple<int, int, double>{ 0, 40, 1.0 },
       std::tuple<int, int, double>{ 0, 50, 1.0 },
       std::tuple<int, int, double>{ 0, 60, 1.0 },
       std::tuple<int, int, double>{ 0, 70, 1.0 },
       std::tuple<int, int, double>{ 0, 80, 1.0 },
       std::tuple<int, int, double>{ 0, 90, 1.0 },
       std::tuple<int, int, double>{ 1, 1, 1.0 },
       std::tuple<int, int, double>{ 1, 11, 1.0 },
       std::tuple<int, int, double>{ 1, 21, 1.0 },
       std::tuple<int, int, double>{ 1, 31, 1.0 },
       std::tuple<int, int, double>{ 1, 41, 1.0 },
       std::tuple<int, int, double>{ 1, 51, 1.0 },
       std::tuple<int, int, double>{ 1, 61, 1.0 },
       std::tuple<int, int, double>{ 1, 71, 1.0 },
       std::tuple<int, int, double>{ 1, 81, 1.0 },
       std::tuple<int, int, double>{ 1, 91, 1.0 },
       std::tuple<int, int, double>{ 2, 2, 1.0 },
       std::tuple<int, int, double>{ 2, 12, 1.0 },
       std::tuple<int, int, double>{ 2, 22, 1.0 },
       std::tuple<int, int, double>{ 2, 32, 1.0 },
       std::tuple<int, int, double>{ 2, 42, 1.0 },
       std::tuple<int, int, double>{ 2, 52, 1.0 },
       std::tuple<int, int, double>{ 2, 62, 1.0 },
       std::tuple<int, int, double>{ 2, 72, 1.0 },
       std::tuple<int, int, double>{ 2, 82, 1.0 },
       std::tuple<int, int, double>{ 2, 92, 1.0 },
       std::tuple<int, int, double>{ 3, 3, 1.0 },
       std::tuple<int, int, double>{ 3, 13, 1.0 },
       std::tuple<int, int, double>{ 3, 23, 1.0 },
       std::tuple<int, int, double>{ 3, 33, 1.0 },
       std::tuple<int, int, double>{ 3, 43, 1.0 },
       std::tuple<int, int, double>{ 3, 53, 1.0 },
       std::tuple<int, int, double>{ 3, 63, 1.0 },
       std::tuple<int, int, double>{ 3, 73, 1.0 },
       std::tuple<int, int, double>{ 3, 83, 1.0 },
       std::tuple<int, int, double>{ 3, 93, 1.0 },
       std::tuple<int, int, double>{ 4, 4, 1.0 },
       std::tuple<int, int, double>{ 4, 14, 1.0 },
       std::tuple<int, int, double>{ 4, 24, 1.0 },
       std::tuple<int, int, double>{ 4, 34, 1.0 },
       std::tuple<int, int, double>{ 4, 44, 1.0 },
       std::tuple<int, int, double>{ 4, 54, 1.0 },
       std::tuple<int, int, double>{ 4, 64, 1.0 },
       std::tuple<int, int, double>{ 4, 74, 1.0 },
       std::tuple<int, int, double>{ 4, 84, 1.0 },
       std::tuple<int, int, double>{ 4, 94, 1.0 },
       std::tuple<int, int, double>{ 5, 5, 1.0 },
       std::tuple<int, int, double>{ 5, 15, 1.0 },
       std::tuple<int, int, double>{ 5, 25, 1.0 },
       std::tuple<int, int, double>{ 5, 35, 1.0 },
       std::tuple<int, int, double>{ 5, 45, 1.0 },
       std::tuple<int, int, double>{ 5, 55, 1.0 },
       std::tuple<int, int, double>{ 5, 65, 1.0 },
       std::tuple<int, int, double>{ 5, 75, 1.0 },
       std::tuple<int, int, double>{ 5, 85, 1.0 },
       std::tuple<int, int, double>{ 5, 95, 1.0 },
       std::tuple<int, int, double>{ 6, 6, 1.0 },
       std::tuple<int, int, double>{ 6, 16, 1.0 },
       std::tuple<int, int, double>{ 6, 26, 1.0 },
       std::tuple<int, int, double>{ 6, 36, 1.0 },
       std::tuple<int, int, double>{ 6, 46, 1.0 },
       std::tuple<int, int, double>{ 6, 56, 1.0 },
       std::tuple<int, int, double>{ 6, 66, 1.0 },
       std::tuple<int, int, double>{ 6, 76, 1.0 },
       std::tuple<int, int, double>{ 6, 86, 1.0 },
       std::tuple<int, int, double>{ 6, 96, 1.0 },
       std::tuple<int, int, double>{ 7, 7, 1.0 },
       std::tuple<int, int, double>{ 7, 17, 1.0 },
       std::tuple<int, int, double>{ 7, 27, 1.0 },
       std::tuple<int, int, double>{ 7, 37, 1.0 },
       std::tuple<int, int, double>{ 7, 47, 1.0 },
       std::tuple<int, int, double>{ 7, 57, 1.0 },
       std::tuple<int, int, double>{ 7, 67, 1.0 },
       std::tuple<int, int, double>{ 7, 77, 1.0 },
       std::tuple<int, int, double>{ 7, 87, 1.0 },
       std::tuple<int, int, double>{ 7, 97, 1.0 },
       std::tuple<int, int, double>{ 8, 8, 1.0 },
       std::tuple<int, int, double>{ 8, 18, 1.0 },
       std::tuple<int, int, double>{ 8, 28, 1.0 },
       std::tuple<int, int, double>{ 8, 38, 1.0 },
       std::tuple<int, int, double>{ 8, 48, 1.0 },
       std::tuple<int, int, double>{ 8, 58, 1.0 },
       std::tuple<int, int, double>{ 8, 68, 1.0 },
       std::tuple<int, int, double>{ 8, 78, 1.0 },
       std::tuple<int, int, double>{ 8, 88, 1.0 },
       std::tuple<int, int, double>{ 8, 98, 1.0 },
       std::tuple<int, int, double>{ 9, 9, 1.0 },
       std::tuple<int, int, double>{ 9, 19, 1.0 },
       std::tuple<int, int, double>{ 9, 29, 1.0 },
       std::tuple<int, int, double>{ 9, 39, 1.0 },
       std::tuple<int, int, double>{ 9, 49, 1.0 },
       std::tuple<int, int, double>{ 9, 59, 1.0 },
       std::tuple<int, int, double>{ 9, 69, 1.0 },
       std::tuple<int, int, double>{ 9, 79, 1.0 },
       std::tuple<int, int, double>{ 9, 89, 1.0 },
       std::tuple<int, int, double>{ 9, 99, 1.0 },
       std::tuple<int, int, double>{ 10, 1, 202.0 },
       std::tuple<int, int, double>{ 10, 2, 404.0 },
       std::tuple<int, int, double>{ 10, 3, 606.0 },
       std::tuple<int, int, double>{ 10, 4, 808.0 },
       std::tuple<int, int, double>{ 10, 5, 1010.0 },
       std::tuple<int, int, double>{ 10, 6, 1212.0 },
       std::tuple<int, int, double>{ 10, 7, 1414.0 },
       std::tuple<int, int, double>{ 10, 8, 1616.0 },
       std::tuple<int, int, double>{ 10, 9, 1818.0 },
       std::tuple<int, int, double>{ 10, 11, -79.0 },
       std::tuple<int, int, double>{ 10, 12, -158.0 },
       std::tuple<int, int, double>{ 10, 13, -237.0 },
       std::tuple<int, int, double>{ 10, 14, -316.0 },
       std::tuple<int, int, double>{ 10, 15, -395.0 },
       std::tuple<int, int, double>{ 10, 16, -474.0 },
       std::tuple<int, int, double>{ 10, 17, -553.0 },
       std::tuple<int, int, double>{ 10, 18, -632.0 },
       std::tuple<int, int, double>{ 10, 19, -711.0 },
       std::tuple<int, int, double>{ 10, 21, 100023.0 },
       std::tuple<int, int, double>{ 10, 22, 200046.0 },
       std::tuple<int, int, double>{ 10, 23, 300069.0 },
       std::tuple<int, int, double>{ 10, 24, 400092.0 },
       std::tuple<int, int, double>{ 10, 25, 500115.0 },
       std::tuple<int, int, double>{ 10, 26, 600138.0 },
       std::tuple<int, int, double>{ 10, 27, 700161.0 },
       std::tuple<int, int, double>{ 10, 28, 800184.0 },
       std::tuple<int, int, double>{ 10, 29, 900207.0 },
       std::tuple<int, int, double>{ 10, 31, -89810.0 },
       std::tuple<int, int, double>{ 10, 32, -179620.0 },
       std::tuple<int, int, double>{ 10, 33, -269430.0 },
       std::tuple<int, int, double>{ 10, 34, -359240.0 },
       std::tuple<int, int, double>{ 10, 35, -449050.0 },
       std::tuple<int, int, double>{ 10, 36, -538860.0 },
       std::tuple<int, int, double>{ 10, 37, -628670.0 },
       std::tuple<int, int, double>{ 10, 38, -718480.0 },
       std::tuple<int, int, double>{ 10, 39, -808290.0 },
       std::tuple<int, int, double>{ 10, 41, -9980.0 },
       std::tuple<int, int, double>{ 10, 42, -19960.0 },
       std::tuple<int, int, double>{ 10, 43, -29940.0 },
       std::tuple<int, int, double>{ 10, 44, -39920.0 },
       std::tuple<int, int, double>{ 10, 45, -49900.0 },
       std::tuple<int, int, double>{ 10, 46, -59880.0 },
       std::tuple<int, int, double>{ 10, 47, -69860.0 },
       std::tuple<int, int, double>{ 10, 48, -79840.0 },
       std::tuple<int, int, double>{ 10, 49, -89820.0 },
       std::tuple<int, int, double>{ 10, 51, 1000.0 },
       std::tuple<int, int, double>{ 10, 52, 2000.0 },
       std::tuple<int, int, double>{ 10, 53, 3000.0 },
       std::tuple<int, int, double>{ 10, 54, 4000.0 },
       std::tuple<int, int, double>{ 10, 55, 5000.0 },
       std::tuple<int, int, double>{ 10, 56, 6000.0 },
       std::tuple<int, int, double>{ 10, 57, 7000.0 },
       std::tuple<int, int, double>{ 10, 58, 8000.0 },
       std::tuple<int, int, double>{ 10, 59, 9000.0 },
       std::tuple<int, int, double>{ 10, 61, 100.0 },
       std::tuple<int, int, double>{ 10, 62, 200.0 },
       std::tuple<int, int, double>{ 10, 63, 300.0 },
       std::tuple<int, int, double>{ 10, 64, 400.0 },
       std::tuple<int, int, double>{ 10, 65, 500.0 },
       std::tuple<int, int, double>{ 10, 66, 600.0 },
       std::tuple<int, int, double>{ 10, 67, 700.0 },
       std::tuple<int, int, double>{ 10, 68, 800.0 },
       std::tuple<int, int, double>{ 10, 69, 900.0 },
       std::tuple<int, int, double>{ 10, 71, 10000.0 },
       std::tuple<int, int, double>{ 10, 72, 20000.0 },
       std::tuple<int, int, double>{ 10, 73, 30000.0 },
       std::tuple<int, int, double>{ 10, 74, 40000.0 },
       std::tuple<int, int, double>{ 10, 75, 50000.0 },
       std::tuple<int, int, double>{ 10, 76, 60000.0 },
       std::tuple<int, int, double>{ 10, 77, 70000.0 },
       std::tuple<int, int, double>{ 10, 78, 80000.0 },
       std::tuple<int, int, double>{ 10, 79, 90000.0 },
       std::tuple<int, int, double>{ 10, 81, 100.0 },
       std::tuple<int, int, double>{ 10, 82, 200.0 },
       std::tuple<int, int, double>{ 10, 83, 300.0 },
       std::tuple<int, int, double>{ 10, 84, 400.0 },
       std::tuple<int, int, double>{ 10, 85, 500.0 },
       std::tuple<int, int, double>{ 10, 86, 600.0 },
       std::tuple<int, int, double>{ 10, 87, 700.0 },
       std::tuple<int, int, double>{ 10, 88, 800.0 },
       std::tuple<int, int, double>{ 10, 89, 900.0 },
       std::tuple<int, int, double>{ 10, 91, -1.0 },
       std::tuple<int, int, double>{ 10, 92, -2.0 },
       std::tuple<int, int, double>{ 10, 93, -3.0 },
       std::tuple<int, int, double>{ 10, 94, -4.0 },
       std::tuple<int, int, double>{ 10, 95, -5.0 },
       std::tuple<int, int, double>{ 10, 96, -6.0 },
       std::tuple<int, int, double>{ 10, 97, -7.0 },
       std::tuple<int, int, double>{ 10, 98, -8.0 },
       std::tuple<int, int, double>{ 10, 99, -9.0 },
       std::tuple<int, int, double>{ 11, 0, 1.0 },
       std::tuple<int, int, double>{ 11, 1, 1.0 },
       std::tuple<int, int, double>{ 11, 2, 1.0 },
       std::tuple<int, int, double>{ 11, 3, 1.0 },
       std::tuple<int, int, double>{ 11, 4, 1.0 },
       std::tuple<int, int, double>{ 11, 5, 1.0 },
       std::tuple<int, int, double>{ 11, 6, 1.0 },
       std::tuple<int, int, double>{ 11, 7, 1.0 },
       std::tuple<int, int, double>{ 11, 8, 1.0 },
       std::tuple<int, int, double>{ 12, 10, 1.0 },
       std::tuple<int, int, double>{ 12, 11, 1.0 },
       std::tuple<int, int, double>{ 12, 12, 1.0 },
       std::tuple<int, int, double>{ 12, 13, 1.0 },
       std::tuple<int, int, double>{ 12, 14, 1.0 },
       std::tuple<int, int, double>{ 12, 15, 1.0 },
       std::tuple<int, int, double>{ 12, 16, 1.0 },
       std::tuple<int, int, double>{ 12, 17, 1.0 },
       std::tuple<int, int, double>{ 12, 18, 1.0 },
       std::tuple<int, int, double>{ 12, 19, 1.0 },
       std::tuple<int, int, double>{ 13, 20, 1.0 },
       std::tuple<int, int, double>{ 13, 21, 1.0 },
       std::tuple<int, int, double>{ 13, 22, 1.0 },
       std::tuple<int, int, double>{ 13, 23, 1.0 },
       std::tuple<int, int, double>{ 13, 24, 1.0 },
       std::tuple<int, int, double>{ 13, 25, 1.0 },
       std::tuple<int, int, double>{ 13, 26, 1.0 },
       std::tuple<int, int, double>{ 13, 27, 1.0 },
       std::tuple<int, int, double>{ 13, 28, 1.0 },
       std::tuple<int, int, double>{ 13, 29, 1.0 },
       std::tuple<int, int, double>{ 14, 30, 1.0 },
       std::tuple<int, int, double>{ 14, 31, 1.0 },
       std::tuple<int, int, double>{ 14, 32, 1.0 },
       std::tuple<int, int, double>{ 14, 33, 1.0 },
       std::tuple<int, int, double>{ 14, 34, 1.0 },
       std::tuple<int, int, double>{ 14, 35, 1.0 },
       std::tuple<int, int, double>{ 14, 36, 1.0 },
       std::tuple<int, int, double>{ 14, 37, 1.0 },
       std::tuple<int, int, double>{ 14, 38, 1.0 },
       std::tuple<int, int, double>{ 14, 39, 1.0 },
       std::tuple<int, int, double>{ 15, 40, 1.0 },
       std::tuple<int, int, double>{ 15, 41, 1.0 },
       std::tuple<int, int, double>{ 15, 42, 1.0 },
       std::tuple<int, int, double>{ 15, 43, 1.0 },
       std::tuple<int, int, double>{ 15, 44, 1.0 },
       std::tuple<int, int, double>{ 15, 45, 1.0 },
       std::tuple<int, int, double>{ 15, 46, 1.0 },
       std::tuple<int, int, double>{ 15, 47, 1.0 },
       std::tuple<int, int, double>{ 15, 48, 1.0 },
       std::tuple<int, int, double>{ 15, 49, 1.0 },
       std::tuple<int, int, double>{ 16, 50, 1.0 },
       std::tuple<int, int, double>{ 16, 51, 1.0 },
       std::tuple<int, int, double>{ 16, 52, 1.0 },
       std::tuple<int, int, double>{ 16, 53, 1.0 },
       std::tuple<int, int, double>{ 16, 54, 1.0 },
       std::tuple<int, int, double>{ 16, 55, 1.0 },
       std::tuple<int, int, double>{ 16, 56, 1.0 },
       std::tuple<int, int, double>{ 16, 57, 1.0 },
       std::tuple<int, int, double>{ 16, 58, 1.0 },
       std::tuple<int, int, double>{ 16, 59, 1.0 },
       std::tuple<int, int, double>{ 17, 60, 1.0 },
       std::tuple<int, int, double>{ 17, 61, 1.0 },
       std::tuple<int, int, double>{ 17, 62, 1.0 },
       std::tuple<int, int, double>{ 17, 63, 1.0 },
       std::tuple<int, int, double>{ 17, 64, 1.0 },
       std::tuple<int, int, double>{ 17, 65, 1.0 },
       std::tuple<int, int, double>{ 17, 66, 1.0 },
       std::tuple<int, int, double>{ 17, 67, 1.0 },
       std::tuple<int, int, double>{ 17, 68, 1.0 },
       std::tuple<int, int, double>{ 17, 69, 1.0 },
       std::tuple<int, int, double>{ 18, 70, 1.0 },
       std::tuple<int, int, double>{ 18, 71, 1.0 },
       std::tuple<int, int, double>{ 18, 72, 1.0 },
       std::tuple<int, int, double>{ 18, 73, 1.0 },
       std::tuple<int, int, double>{ 18, 74, 1.0 },
       std::tuple<int, int, double>{ 18, 75, 1.0 },
       std::tuple<int, int, double>{ 18, 76, 1.0 },
       std::tuple<int, int, double>{ 18, 77, 1.0 },
       std::tuple<int, int, double>{ 18, 78, 1.0 },
       std::tuple<int, int, double>{ 18, 79, 1.0 },
       std::tuple<int, int, double>{ 19, 80, 1.0 },
       std::tuple<int, int, double>{ 19, 81, 1.0 },
       std::tuple<int, int, double>{ 19, 82, 1.0 },
       std::tuple<int, int, double>{ 19, 83, 1.0 },
       std::tuple<int, int, double>{ 19, 84, 1.0 },
       std::tuple<int, int, double>{ 19, 85, 1.0 },
       std::tuple<int, int, double>{ 19, 86, 1.0 },
       std::tuple<int, int, double>{ 19, 87, 1.0 },
       std::tuple<int, int, double>{ 19, 88, 1.0 },
       std::tuple<int, int, double>{ 19, 89, 1.0 },
       std::tuple<int, int, double>{ 20, 90, 1.0 },
       std::tuple<int, int, double>{ 20, 91, 1.0 },
       std::tuple<int, int, double>{ 20, 92, 1.0 },
       std::tuple<int, int, double>{ 20, 93, 1.0 },
       std::tuple<int, int, double>{ 20, 94, 1.0 },
       std::tuple<int, int, double>{ 20, 95, 1.0 },
       std::tuple<int, int, double>{ 20, 96, 1.0 },
       std::tuple<int, int, double>{ 20, 97, 1.0 },
       std::tuple<int, int, double>{ 20, 98, 1.0 },
       std::tuple<int, int, double>{ 20, 99, 1.0 },
   };
   Vec<std::string> rnames{
       "SOS0", "SOS1", "SOS2", "SOS3",     "SOS4", "SOS5", "SOS6",
       "SOS7", "SOS8", "SOS9", "BILANCIO", "SOSA", "SOSB", "SOSC",
       "SOSD", "SOSE", "SOSF", "SOSG",     "SOSH", "SOSI", "SOSL",
   };
   Vec<std::string> cnames{
       "A0",  "A1",  "A2",  "A3",  "A4",  "A5",  "A6",  "A7",  "A8",  "A9",
       "B0",  "B1",  "B2",  "B3",  "B4",  "B5",  "B6",  "B7",  "B8",  "B9",
       "C0",  "C1",  "C2",  "C3",  "C4",  "C5",  "C6",  "C7",  "C8",  "C9",
       "D0",  "D1",  "D2",  "D3",  "D4",  "D5",  "D6",  "D7",  "D8",  "D9",
       "E_0", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9",
       "F0",  "F1",  "F2",  "F3",  "F4",  "F5",  "F6",  "F7",  "F8",  "F9",
       "G0",  "G1",  "G2",  "G3",  "G4",  "G5",  "G6",  "G7",  "G8",  "G9",
       "H0",  "H1",  "H2",  "H3",  "H4",  "H5",  "H6",  "H7",  "H8",  "H9",
       "I0",  "I1",  "I2",  "I3",  "I4",  "I5",  "I6",  "I7",  "I8",  "I9",
       "L0",  "L1",  "L2",  "L3",  "L4",  "L5",  "L6",  "L7",  "L8",  "L9",
   };
   int nCols = 100;
   int nRows = 21;
   ProblemBuilder<double> pb;
   pb.reserve( 289, 21, 100 );
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
   pb.setProblemName( "enigma.hpp" );
   Problem<double> problem = pb.build();
   /// PROBLEM BUILDER CODE END

   return problem;
}

} // namespace instances
} // namespace papilo

#endif
