/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    PaPILO --- Parallel Presolve for Integer and Linear Optimization       */
/*                                                                           */
/* Copyright (C) 2020-2024 Zuse Institute Berlin (ZIB)                       */
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

#ifndef _PAPILO_PRESOLVERS_CliqueMerging_HPP_
#define _PAPILO_PRESOLVERS_CliqueMerging_HPP_

#include "papilo/core/PresolveMethod.hpp"
#include "papilo/core/Problem.hpp"
#include "papilo/core/ProblemUpdate.hpp"
#include "papilo/external/pdqsort/pdqsort.h"
#include <algorithm>
#include <map>
#include <set>
#include <vector>
#ifdef PAPILO_TBB
#include "papilo/misc/tbb.hpp"
#endif

namespace papilo
{

template <typename REAL>
class CliqueMerging : public PresolveMethod<REAL>
{
 public:
   CliqueMerging() : PresolveMethod<REAL>()
   {
      this->setName( "cliquemerging" );
      this->setTiming( PresolverTiming::kMedium );
      this->setType( PresolverType::kIntegralCols );
   }

   PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions, const Timer& timer,
            int& reason_of_infeasibility ) override;

   std::pair<std::set<int>, Vec<int>>
   greedyClique( const ConstraintMatrix<REAL>& matrix,
                 const std::set<std::pair<int, int>>& edges,
                 const std::map<int, std::set<int>>& Neighbourhoodlists,
                 int clique );

   bool
   isCovered( const ConstraintMatrix<REAL>& matrix, int row,
              const std::set<int>& newClique );
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class CliqueMerging<double>;
extern template class CliqueMerging<Quad>;
extern template class CliqueMerging<Rational>;
#endif

template <typename REAL>
std::pair<std::set<int>, Vec<int>>
CliqueMerging<REAL>::greedyClique(
    const ConstraintMatrix<REAL>& matrix,
    const std::set<std::pair<int, int>>& edges,
    const std::map<int, std::set<int>>& Neighbourhoodlists, int clique )
{
   std::set<int> newClique;
   Vec<int> newVertices;
   newClique.reserve(300);
   newVertices.reserve(200);
   assert( clique >= 0 );
   assert( clique < matrix.getNRows() );
   assert( clique >= 0 );
   assert( clique < matrix.getNRows() );
   const auto rowvec = matrix.getRowCoefficients( clique );
   const auto indices = rowvec.getIndices();

   for( int vertexIndex = 0; vertexIndex < rowvec.getLength(); ++vertexIndex )
   {
      newClique.emplace( indices[vertexIndex] );
   }
   std::set<int> Neighbourhood = Neighbourhoodlists.at( *newClique.begin() );
   for( int potentialNewVertex : Neighbourhood )
   {
      if( newClique.find( potentialNewVertex ) != newClique.end() )
         continue;
      bool isIn = true;
      for( int otherMember : newClique )
      {
         if( edges.find( std::pair<int, int>{ potentialNewVertex,
                                              otherMember } ) == edges.end() )
         {
            isIn = false;
            break;
         }
      }
      if( isIn )
      {
         newClique.emplace( potentialNewVertex );
         newVertices.emplace_back( potentialNewVertex );
      }
   }
   std::pair<std::set<int>, Vec<int>> output;
   output.first = newClique;
   output.second = newVertices;
   return output;
}

template <typename REAL>
bool
CliqueMerging<REAL>::isCovered( const ConstraintMatrix<REAL>& matrix,
                                int row, const std::set<int>& newClique )
{
   assert( row >= 0 );
   assert( row < matrix.getNRows() );
   const auto rowVec = matrix.getRowCoefficients( row );
   const auto rowInds = rowVec.getIndices();

   for( int vertexIndex = 0; vertexIndex < rowVec.getLength(); ++vertexIndex )
   {
      if( newClique.find( rowInds[vertexIndex] ) == newClique.end() )
         return false;
   }
   return true;
}

template <typename REAL>
PresolveStatus
CliqueMerging<REAL>::execute( const Problem<REAL>& problem,
                              const ProblemUpdate<REAL>& problemUpdate,
                              const Num<REAL>& num,
                              Reductions<REAL>& reductions, const Timer& timer,
                              int& reason_of_infeasibility )
{
   const auto& matrix = problem.getConstraintMatrix();

   const auto lb = problem.getLowerBounds();

   const auto ub = problem.getUpperBounds();

   std::vector<RowFlags> rowFlags = matrix.getRowFlags();

   const int nrows = matrix.getNRows();

   PresolveStatus result = PresolveStatus::kUnchanged;

   Vec<int> Cliques;
   Cliques.reserve(100000);

   std::set<std::pair<int, int>> edges;
   edges.reserve(1000000);

   std::map<int, std::set<int>> neighbourLists;
   Neighbourhoodlists.reserve(2000000);

   std::set<int> Vertices;
   Vertices.reserve(500000);

   for( int row = 0; row < nrows; ++row )
   {
      assert( row >= 0 );
      assert( row < matrix.getNRows() );
      const auto cliqueRow = matrix.getRowCoefficients( row );
      const auto cliqueIndices = cliqueRow.getIndices();

      std::tuple<bool, bool, bool> cliqueCheck =
          problem.is_clique_equation_or_sos1( matrix, row, num );
      if( !matrix.isRowRedundant( row ) &&
          std::get<0>( cliqueCheck ) & !std::get<1>( cliqueCheck ) &&
          !std::get<2>( cliqueCheck ) && cliqueRow.getLength() < 100 )
      {
         Cliques.push_back( row );
         rowFlags[row].set( RowFlag::kClique );
         for( int col = 0; col < cliqueRow.getLength(); ++col )
         {
            int vertex = cliqueIndices[col];
            Vertices.emplace( vertex );
            if( neighbourLists.find( vertex ) == neighbourLists.end() )
            {
               std::set<int> neighOfCol;
               neighbourLists.emplace( vertex, neighOfCol );
            }
            for( int otherCol = 0; otherCol < col; ++otherCol )
            {
               int otherVertex = cliqueIndices[otherCol];
               edges.emplace( std::pair<int, int>{ otherVertex, vertex } );
               edges.emplace( std::pair<int, int>{ vertex, otherVertex } );
               neighbourLists[vertex].emplace( otherVertex );
               neighbourLists[otherVertex].emplace( vertex );
            }
         }
      }
      else
      {
         rowFlags[row].unset( RowFlag::kClique );
      }
#ifdef PAPILO_TBB
      if( edges.size() > 1000000 )
#else
      if( edges.size() > 100000 )
#endif
         break;
   }

#ifdef PAPILO_TBB
   tbb::combinable<Vec<int>> completedCliquesComb;
#else   
   Vec<int> completedCliques;
#endif

#ifdef PAPILO_TBB
   tbb::combinable<Vec<std::pair<Vec<int>,Vec<int>>>> reductionResults;
   tbb::parallel_for(
       tbb::blocked_range<int>( 0, Cliques.end() - Cliques.begin() ),
       [&]( const tbb::blocked_range<int>& r )
       {
          for( int cliqueInd = r.begin(); cliqueInd != r.end(); ++cliqueInd )
#else
   for( int cliqueInd = 0; cliqueInd < Cliques.end() - Cliques.begin() ;
        ++cliqueInd )
#endif
      {
#ifdef PAPILO_TBB
      Vec<int> completedCliques = completedCliquesComb.combine( [](const Vec<int>& a, const Vec<int>& b) {
         std::vector<int> result = a;
         result.insert(result.end(), b.begin(), b.end() );
         return result;
      } );
#endif
      int clique = Cliques[cliqueInd];
      if( cliqueInd > 10000 )
         break;
      if( std::find( completedCliques.begin(), completedCliques.end(),
                     clique ) != completedCliques.end() )
         continue;

      const std::pair<std::set<int>, Vec<int>> resultOfSearch =
          greedyClique( matrix, edges, neighbourLists, clique );

      const std::set<int> newClique = resultOfSearch.first;

      const Vec<int> newVertices = resultOfSearch.second;

#ifdef PAPILO_TBB
      std::pair<Vec<int>,Vec<int>> reductionResult;
      reductionResult.first = newVertices ;
      reductionResult.first.push_back( clique );
#endif

      Vec<int> coveredCliques;
      for( int cl = 0; cl < Cliques.end() - Cliques.begin(); ++cl )
      {
         if( cl != cliqueInd && isCovered( matrix, Cliques[cl], newClique ) )
         {
            coveredCliques.push_back( cl );
            completedCliques.push_back( Cliques[cl] );
         }
      }

#ifdef PAPILO_TBB
      reductionResult.second = coveredCliques;
#endif


      if( coveredCliques.end() - coveredCliques.begin() > 0 )
      {

#ifdef PAPILO_TBB
         reductionResults.local().push_back( reductionResult );
#else

         result = PresolveStatus::kReduced;
         TransactionGuard<REAL> tg{ reductions };
         for( int covRow = 0;
              covRow < coveredCliques.end() - coveredCliques.begin(); ++covRow )
         {
            reductions.lockRow( Cliques[coveredCliques[covRow]] );
         }
         reductions.lockRow( clique );
         /*
         std::cout << "\nLocked Cols: "; */
         for( std::set<int>::iterator vertexIndex = newClique.begin();
              vertexIndex != newClique.end(); ++vertexIndex )
         {
            reductions.lockCol( *vertexIndex );
            reductions.lockColBounds( *vertexIndex );/*
            std::cout << *vertexIndex;
            std::cout << " ";*/
         }
         assert( clique >= 0 );
         assert( clique < matrix.getNRows() );
         auto rowVector = matrix.getRowCoefficients( clique );
         auto rowValues = rowVector.getValues();
         auto rowInds = rowVector.getIndices();
         auto val = rowValues[0] * ( ub[rowInds[0]] - abs( lb[rowInds[0]] ) );
         for( int vertexIndex = 0;
              vertexIndex < newVertices.end() - newVertices.begin();
              ++vertexIndex )
         {
            reductions.changeMatrixEntry(
                clique, newVertices[vertexIndex],
                val * ( ub[newVertices[vertexIndex]] -
                        lb[newVertices[vertexIndex]] ) );
            assert( !num.isEq( val * ( ub[newVertices[vertexIndex]] -
                                       lb[newVertices[vertexIndex]] ),
                               0.0 ) );
         }
         for( int row = 0; row < coveredCliques.end() - coveredCliques.begin();
              ++row )
         {
            reductions.markRowRedundant( Cliques[coveredCliques[row]] );
         }
#endif
      }
   }
#ifdef PAPILO_TBB
   });
   Vec<std::pair<Vec<int>,Vec<int>>> reductionResultsComb = reductionResults.combine( 
      [](const Vec<std::pair<Vec<int>,Vec<int>>>& a, const Vec<std::pair<Vec<int>,Vec<int>>>& b) {
      Vec<std::pair<Vec<int>,Vec<int>>> result = a;
      result.insert(result.end(), b.begin(), b.end() );
      return result;
   } );
   pdqsort(reductionResultsComb.begin(), reductionResultsComb.end());
   int totalreductionsize = reductionResultsComb.end() - reductionResultsComb.begin();
   for( int reductionsize = reductionResultsComb.end() - reductionResultsComb.begin(); reductionsize > 0; reductionsize /= 2 )
   {
      for( int transactionindex = 0; transactionindex * reductionsize < reductionResultsComb.end() - reductionResultsComb.begin(); ++transactionindex )
      {
         result = PresolveStatus::kReduced;
         TransactionGuard<REAL> tg{ reductions };
         for( int reductionindex = transactionindex * reductionsize; reductionindex < std::min((transactionindex+1)* reductionsize,totalreductionsize); ++reductionindex )
         {
            Vec<int> covCliques = reductionResultsComb[reductionindex].second;
            int clique = reductionResultsComb[reductionindex].first.back();
            Vec<int> newVertices = reductionResultsComb[reductionindex].first;
            auto cliqueRow = matrix.getRowCoefficients( clique );
            auto cliqueIndices = cliqueRow.getIndices();
            for( int cliqueIndex = 0; cliqueIndex < cliqueRow.getLength(); ++cliqueIndex )
            {
               reductions.lockCol( cliqueIndices[cliqueIndex] );
               reductions.lockColBounds( cliqueIndices[cliqueIndex] );
            }
            for( int vertexIndex = 0; vertexIndex < newVertices.end() - newVertices.begin() - 1; ++vertexIndex )
            {
               reductions.lockCol( newVertices[vertexIndex] );
               reductions.lockColBounds( newVertices[vertexIndex] );
            }
            for( int rowIndex = 0; rowIndex < covCliques.end() - covCliques.begin(); ++rowIndex )
            {
               reductions.lockRow( Cliques[covCliques[rowIndex]] );
            }
            reductions.lockRow( clique );
         }
         for( int reductionindex = transactionindex * reductionsize; reductionindex < std::min((transactionindex+1)* reductionsize,totalreductionsize); ++reductionindex )
         {
            int clique = reductionResultsComb[reductionindex].first.back();
            Vec<int> newVertices = reductionResultsComb[reductionindex].first;
            auto rowVector = matrix.getRowCoefficients( clique );
            auto rowValues = rowVector.getValues();
            auto rowInds = rowVector.getIndices();
            auto val = rowValues[0] * ( ub[rowInds[0]] - abs( lb[rowInds[0]] ) );
            for( int vertexIndex = 0;
                  vertexIndex < newVertices.end() - newVertices.begin() - 1;
                  ++vertexIndex )
               { 
               reductions.changeMatrixEntry(
                  clique, newVertices[vertexIndex],
                  val * ( ub[newVertices[vertexIndex]] -
                           lb[newVertices[vertexIndex]] ) );
               }
         }
         for( int reductionindex = transactionindex * reductionsize; reductionindex < std::min((transactionindex+1)* reductionsize,totalreductionsize); ++reductionindex )
         {
            Vec<int> covCliques = reductionResultsComb[reductionindex].second;
            for( int row = 0; row < covCliques.end() - covCliques.begin();
               ++row )
               {reductions.markRowRedundant( Cliques[covCliques[row]] );}
         }
      }
   }
#endif
   this->setEnabled( false );
   return result;
}

} // namespace papilo

#endif