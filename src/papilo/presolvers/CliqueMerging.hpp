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
#include <vector>
#include <algorithm>
#include <map>
#include <set>
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
      this->setTiming( PresolverTiming::kExhaustive );
      this->setType( PresolverType::kIntegralCols );
   }

   PresolveStatus
   execute( const Problem<REAL>& problem,
            const ProblemUpdate<REAL>& problemUpdate, const Num<REAL>& num,
            Reductions<REAL>& reductions, const Timer& timer,
            int& reason_of_infeasibility ) override;

   std::pair<std::set<int>,Vec<int>>
   greedyClique( const ConstraintMatrix<REAL>& matrix, const std::set<std::pair<int,int>> Edges,
               const std::map<int,std::set<int>> Neighbourhoodlists, const int clique );

   bool
   isCovered( const ConstraintMatrix<REAL>& matrix, 
            const int row, const std::set<int> newClique );

};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class CliqueMerging<double>;
extern template class CliqueMerging<Quad>;
extern template class CliqueMerging<Rational>;
#endif
            
template <typename REAL>
std::pair<std::set<int>,Vec<int>>
CliqueMerging<REAL>::greedyClique( const ConstraintMatrix<REAL>& matrix, const std::set<std::pair<int,int>> Edges,
               const std::map<int,std::set<int>> Neighbourhoodlists, const int clique )
{
   std::set<int> newClique;
   Vec<int> newVertices;

   const auto rowvec = matrix.getRowCoefficients( clique );
   const auto indices = rowvec.getIndices();

   for( int vertexIndex = 0; vertexIndex < rowvec.getLength(); ++vertexIndex)
   {
      newClique.emplace(indices[vertexIndex]);
   }
   std::set<int> Neighbourhood = Neighbourhoodlists.at(*newClique.begin());
   for( std::set<int>::iterator vertexIndex = Neighbourhood.begin(); vertexIndex != Neighbourhood.end(); ++vertexIndex )
   {
      int potentialNewVertex = *vertexIndex;
      if( newClique.find(potentialNewVertex) != newClique.end() )
         continue;
      bool isIn = true;
      for( std::set<int>::iterator otherMemberIndex = newClique.begin(); otherMemberIndex != newClique.end(); ++otherMemberIndex )
      {
         int otherMember = *otherMemberIndex;
         if( Edges.find( std::pair<int,int> {potentialNewVertex, otherMember}) == Edges.end() )
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
   std::pair<std::set<int>,Vec<int>> output;
   output.first = newClique;
   output.second = newVertices;
   return output;
}

template <typename REAL>
bool
CliqueMerging<REAL>::isCovered( const ConstraintMatrix<REAL>& matrix, 
         const int row, const std::set<int> newClique )
{
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
/*
     TransactionGuard<REAL> tg{ reductions };
   auto row = matrix.getRowCoefficients(100);
   std::cout << "\n";
   for( int i = 0; i < row.getLength(); ++i )
   {
      std::cout << row.getIndices()[i];
      std::cout << " ";
   }
   reductions.lockRow(100);
   reductions.lockCol(row.getIndices()[row.getLength()-1]+1);
   reductions.changeMatrixEntry(100, row.getIndices()[row.getLength()-1]-3, 1.0);
   PresolveStatus result1 = PresolveStatus::kReduced;
   return result1;*/

   const auto lb = problem.getLowerBounds();

   const auto ub = problem.getUpperBounds();

   std::vector<RowFlags> rowFlags = matrix.getRowFlags();

   const int nrows = matrix.getNRows();

   PresolveStatus result = PresolveStatus::kUnchanged;

   Vec<int> Cliques;

   std::set<std::pair<int,int>> Edges;

   std::map<int,std::set<int>> Neighbourlists;

   std::set<int> Vertices;

   for( int row = 0; row < nrows; ++row )
   {
      const auto cliqueRow = matrix.getRowCoefficients( row );
      const auto cliqueIndices = cliqueRow.getIndices();

      std::tuple<bool, bool, bool> cliqueCheck =
          problem.is_clique_equation_or_sos1( matrix, row, num );
      if( std::get<0>(cliqueCheck) & !std::get<1>(cliqueCheck) && !std::get<2>(cliqueCheck) && cliqueRow.getLength() < 100 )
      { /*
         std::cout << "\nClique: ";
         std::cout << matrix.getLeftHandSides()[row];
         std::cout << " ";
         std::cout << matrix.getRightHandSides()[row];
         std::cout << "\n";*/
         Cliques.push_back( row );
         rowFlags[row].set( RowFlag::kClique );
         for( int col = 0; col < cliqueRow.getLength(); ++col )
         {/*
            std::cout << cliqueIndices[col];
            std::cout << " ";
            std::cout << cliqueRow.getValues()[col] * (ub[cliqueIndices[col]]-abs(lb[cliqueIndices[col]]));
            std::cout << " ";*/
            int vertex = cliqueIndices[col];
            Vertices.emplace(vertex);
            if( Neighbourlists.find(vertex) == Neighbourlists.end() )
            {
               std::set<int> neighOfCol;
               Neighbourlists.emplace(vertex, neighOfCol);
            }
            for( int otherCol = 0; otherCol < col; ++otherCol )
            {
               int otherVertex = cliqueIndices[otherCol];
               Edges.emplace(std::pair<int,int> {otherVertex,vertex});
               Edges.emplace(std::pair<int,int> {vertex,otherVertex});
               Neighbourlists[vertex].emplace(otherVertex);
               Neighbourlists[otherVertex].emplace(vertex);
            }
         }
      }
      else
      {
         rowFlags[row].unset( RowFlag::kClique );
      }
   }

   Vec<bool> completedCliques(Cliques.end() - Cliques.begin(), false);

   for( int cliqueInd = 0; cliqueInd < Cliques.end() - Cliques.begin(); ++cliqueInd )
   {
      int clique = Cliques[cliqueInd];

      if( completedCliques[cliqueInd])
         continue;
      
      std::pair<std::set<int>,Vec<int>> resultOfSearch = greedyClique( matrix, Edges, Neighbourlists, clique );

      std::set<int> newClique = resultOfSearch.first;

      Vec<int> newVertices = resultOfSearch.second;

      Vec<int> coveredCliques;
      for( int cl = 0; cl < Cliques.end() - Cliques.begin(); ++cl )
      {
         if( cl != cliqueInd && isCovered( matrix, Cliques[cl], newClique) )
            coveredCliques.push_back( cl );
      }
      if( coveredCliques.end() - coveredCliques.begin() > 0 && newVertices.end() - newVertices.begin() > 0 )
      {
         result = PresolveStatus::kReduced;
         TransactionGuard<REAL> tg{ reductions };
         for( int cl = 0; cl < Cliques.end() - Cliques.begin(); ++cl )
         {
            reductions.lockRow( Cliques[cl] );
         }/*
         for( std::set<int>::iterator vertexIndex = newClique.begin(); vertexIndex != newClique.end(); ++vertexIndex )
         {
            reductions.lockCol( *vertexIndex );
            reductions.lockColBounds( *vertexIndex );
         }*/
         for( auto i = Vertices.begin(); i != Vertices.end(); ++i )
         {
            reductions.lockCol( *i );
            reductions.lockColBounds( *i );
         }/*
         for( int i = 0; i < newVertices.end() - newVertices.begin(); ++i )
         {
            reductions.lockCol( newVertices[i] );
            reductions.lockColBounds( newVertices[i] );
         }*/
        std::cout << "\nCovered Cliques:\n";
         for( int row = 0; row < coveredCliques.end() - coveredCliques.begin(); ++row )
         {
            std::cout << coveredCliques[row];
            reductions.markRowRedundant(Cliques[coveredCliques[row]]);
            completedCliques[coveredCliques[row]] = true;
         }
         auto rowVector = matrix.getRowCoefficients( clique );
         auto rowValues = rowVector.getValues();
         auto rowInds = rowVector.getIndices();
         auto val = rowValues[0] * ( ub[rowInds[0]] - abs(lb[rowInds[0]]) );
         /*
         for( auto i = newClique.begin(); i != newClique.end(); ++i )
         {
            for( auto j = newClique.begin(); j != i; ++j)
            {
               std::pair<int,int> edge;
               edge.first = *i;
               edge.second = *j;
               assert(Edges.find( edge) != Edges.end());
               edge.first = *j;
               edge.second = *i;
               assert(Edges.find( edge) != Edges.end());
               bool shareClique = false;
               for( int r = 0; r < Cliques.end() - Cliques.begin(); ++r)
               {
                  auto c = Cliques[r];
                  auto v = matrix.getRowCoefficients(c);
                  auto s1 = v.getIndices();
                  Vec<int> s;
                  for( int i1 = 0; i1< v.getLength(); ++i1)
                  {
                     s.push_back(s1[i1]);
                  }
                  auto va = v.getValues();
                  if( std::find( s.begin(), s.end(), *i ) != s.end() && std::find( s.begin(), s.end(), *j ) != s.end() )
                  {
                     int index1 = 0;
                     int index2 = 0;
                     while( s[index1] != *i)
                     {
                        index1 += 1;
                     }
                     while( s[index2] != *j)
                     {
                        index2 += 1;
                     }
                     if( num.isLE(matrix.getLeftHandSides()[c],va[index1]) &&
                     num.isLE(va[index1],matrix.getRightHandSides()[c]) &&
                     num.isLE(matrix.getLeftHandSides()[c],va[index2]) &&
                     num.isLE(va[index2],matrix.getRightHandSides()[c]) &&
                     ( (num.isLT(va[index2]+va[index1],matrix.getLeftHandSides()[c]) &&
                     num.isLE(0.0,matrix.getRightHandSides()[c])) ||
                     (num.isLT(matrix.getRightHandSides()[c],va[index2]+va[index1])
                     && num.isLE(matrix.getLeftHandSides()[c],0.0))))
                     {
                        shareClique = true;
                        break;
                     }/*
                    assert(s[index1] == *i && s[index2] == *j && va[index1] == 1.0 && va[index2] == 1.0 && matrix.getRightHandSides()[c] == 1.0);
                    assert(ub[*i] == 1.0 && ub[*j] == 1.0 && lb[*i] == 0.0  );
                    assert(lb[*j] == 0.0 );
                    assert(matrix.getLeftHandSides()[c] <= 0.0);
                    //assert(rowFlags[c].test( RowFlag::kLhsInf));
                     if( s[index1] == *i && s[index2] == *j && va[index1] == 1.0 && va[index2] == 1.0 && matrix.getRightHandSides()[c] == 1.0 &&
                     ub[*i] == 1.0 && ub[*j] == 1.0 && lb[*i] == 0.0 && lb[*j] == 0.0  )
                     {
                        shareClique = true;
                        break;
                     }*//*
                  }
               }
               assert( shareClique );
            }
         }
         std::set<int> Verifier;
         for( int i = 0; i < rowVector.getLength(); ++i )
         {
            Verifier.emplace( rowInds[i] );
         }
         for( int i = 0; i < newVertices.end() - newVertices.begin(); ++i )
         {
            Verifier.emplace(newVertices[i]);
         }
         assert( Verifier == newClique );
         std::cout << "\n new Clique Vertices for the following clique:\n";
         std::cout << cliqueInd;
         std::cout << " : ";*/
         for( int vertexIndex = 0; vertexIndex < newVertices.end() - newVertices.begin(); ++vertexIndex )
         { /*
            std::cout << newVertices[vertexIndex];
            std::cout << " ";
            std::cout << val * ( ub[newVertices[vertexIndex]] - lb[newVertices[vertexIndex]] );
            std::cout << " ";
            std::cout << ub[newVertices[vertexIndex]];
            std::cout << " ";*/
            reductions.changeMatrixEntry( clique, newVertices[vertexIndex], val * ( ub[newVertices[vertexIndex]] - lb[newVertices[vertexIndex]] ) );/*
            assert( !num.isEq(val * ( ub[newVertices[vertexIndex]] - abs(lb[newVertices[vertexIndex]]) ),0.0) );
            for( int i = 0; i < rowVector.getLength(); ++i )
            {
               assert( rowInds[i] != newVertices[vertexIndex]);
            }
            for( int i = 0; i < matrix.getRowCoefficients( clique ).getLength(); ++i )
            {
               assert(matrix.getRowCoefficients( clique ).getIndices()[i] != newVertices[vertexIndex]);
            }*/
         }
      }
   }
   return result;
}

} // namespace papilo

#endif