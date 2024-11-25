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
#include "papilo/core/SingleRow.hpp"
#include "papilo/misc/Array.hpp"
#include "papilo/misc/Hash.hpp"
#include "papilo/misc/Vec.hpp"
#include "papilo/misc/compress_vector.hpp"
#include "papilo/misc/fmt.hpp"
#include <vector>
#include <algorithm>
#ifdef PAPILO_TBB
#include "papilo/misc/tbb.hpp"
#endif
#include <atomic>
#include <boost/functional/hash.hpp>

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
            const ProblemUpdate<REAL>& problemUpdate,
            const Num<REAL>& num, Reductions<REAL>& reductions,
            const Timer& timer, int& reason_of_infeasibility) override;
            
    std::vector<std::vector<int>>
    maxClique( std::vector<std::vector<int>> Neighbourhoods, std::vector<std::vector<int>> maximalCLiques, 
    std::vector<int> OldClique, std::vector<int> Neighbourhood,
    int miniteration, int size, int n );
};

#ifdef PAPILO_USE_EXTERN_TEMPLATES
extern template class CliqueMerging<double>;
extern template class CliqueMerging<Quad>;
extern template class CliqueMerging<Rational>;
#endif

template <typename REAL>
std::vector<std::vector<int>>
CliqueMerging<REAL>::maxClique( std::vector<std::vector<int>> Neighbourhoods, std::vector<std::vector<int>> maximalCLiques, 
    std::vector<int> OldClique, std::vector<int> Neighbourhood,
    int miniteration, int size, int n ){
    if( Neighbourhood.size() == 0)
    {
        maximalCLiques.push_back( OldClique );
    }
    else if( miniteration < n ) // && Neighbourhood.size() + OldClique.size() >= size )
    {
        std::vector<std::vector<int>> maximalCliques1 = maxClique( Neighbourhoods, maximalCLiques, OldClique, Neighbourhood, miniteration + 1, size, n );
        std::vector<std::vector<int>> maximalCliques2 = {};
        if( std::find( OldClique.begin(), OldClique.end(), miniteration ) == OldClique.end() 
         && std::find( Neighbourhood.begin(), Neighbourhood.end(), miniteration ) != Neighbourhood.end() )
        {
            OldClique.push_back( miniteration );
            for( int vertex = 0; vertex < Neighbourhood.size(); ++vertex ){
                if( std::find( Neighbourhoods[miniteration].begin, Neighbourhoods[miniteration].end(), Neighbourhood[vertex] ) 
                == Neighbourhoods[miniteration].end())
                {
                    Neighbourhood.erase(Neighbourhood.begin() + vertex);
                    vertex -= 1;
                }
            }
            maximalCliques2 = maxClique( Neighbourhoods, maximalCliques, OldClique, Neighbourhood, miniteration + 1, size, n );
        }
        for( int maxClique = 0; maxClique < maximalCliques2.size(); ++ maxClique ){
            if( std::find( maximalCliques1.begin(), maximalCliques1.end(), maximalCliques2[maxClique]) == maximalCliques1.end() )
            {
                maximalCliques1.push_back(maximalCliques2[maxClique]);
            }
        }
        maximalCliques = maximalCliques1;
    }
    return(maximalCliques)
}

template <typename REAL>
PresolveStatus
CliqueMerging<REAL>::execute( const Problem<REAL>& problem,
                              const ProblemUpdate<REAL>& problemUpdate,
                              const Num<REAL>& num, Reductions<REAL>& reductions,
                              const Timer& timer, int& reason_of_infeasibility){

    const auto& matrix = problem.getConstraintMatrix();

    std::vector<RowFlags> rowFlags = matrix.getRowFlags();

    const int nrows = consMatrix.getNRows();

    PresolveStatus result = PresolveStatus::kUnchanged;

    std::vector<int> Cliques;

    for( int row = 0; row < nrows; ++row )
    {
        if( rowFlags[row].test( RowFlag::kClique ) && problem.is_clique_or_sos1(
                    matrix, row, num ) == {true, false})
        {
            Cliques.push_back(row)
        }
    }

    std::vector<int> Vertices;
    std::vector<std::vector<int>> Neighbourhoods;

    int nCliques = Cliques.size();

    for( int clique = 0; clique < nCliques; ++clique)
    {
        auto rowvec = matrix.getRowCoefficients(clique);
        auto indices = rowvec.getIndices();
        for( int v = 0; v < rowvec.getLength(); ++v )
        {
            int posOfv =  std::distance(Vertices.begin(), std::find(Vertices.begin(), Vertices.end(), indices[v]));
            if( posOfv == std::distance(Vertices.begin(), Vertices.end()) )
            {
                Vertices.push_back(indices[v]);
                Neighbourhoods.push_back({});
            }
            for( int w = 0; w < v; ++w )
            {
                int posOfw = std::distance(Vertices.begin(), std::find(Vertices.begin(), Vertices.end(), indices[w]));
                if( std::find(Neighbourhoods[posOfv].begin(), Neighbourhoods[posOfv].end(), posOfw) == Vertices.end() )
                {
                    Neighbourhoods[posOfv].push_back(posOfw);
                    Neighbourhoods[posOfw].push_back(posOfv);
                }
            }
        }
    }

    bool completedCliques[nCliques];
    completedCliques.fill(false);

    for( int clique = 0; clique < nCliques; ++clique )
    {
        if( completedCliques[clique] )
            continue;
        
        std::vector<std::vector<int>> maximalCliques = {};
        std::vector<int> OldClique;
        auto rowvec = matrix.getRowCoefficients(Cliques[clique]);
        auto indices = rowvec.getIndices();
        OldClique.push_back( std::distance(Vertices.begin(), std::find( Vertices.begin(), Vertices.end(), indices[0]) ));
        std::vector<int> Neighbourhood = Neighbourhoods[std::distance(Vertices.begin(), std::find( Vertices.begin(), Vertices.end(), indices[0]))];

        for( int i = 0; i < Neighbourhood.size(); ++i )
        {
            if( std::find(indices.begin(), indices.end(), Vertices[Neighbourhood[i]]) != indices.end() )
            {
                Neighbourhood.erase(Neighbourhood.begin() + i);
                i -= 1;
            }
        }

        for(int v = 1; v < rowvec.getLength(); ++v)
        {
            int vertex = std::distance(Vertices.begin(), std::find( Vertices.begin(), Vertices.end(), indices[v]));
            OldClique.push_back(vertex);
            for( int w = 0; w < Neighbourhood.size(); ++w )
            {
                if( std::find(Neighbourhoods[vertex].begin(), Neighbourhoods[vertex].end(), Neighbourhood[w]) == Neighbourhoods[vertex].end() )
                {
                    Neighbourhood.erase(Neighbourhood.begin() + w);
                    w -= 1;
                }
            }
        }

        int miniteration = 0;
        int size = OldClique.size();
        int n = Vertices.size();
        maximalCliques = maxClique( Neighbourhoods, maximalCliques, OldClique, Neighbourhood,
        miniteration, size, n );

        std::vector<int> maxCoveredCliques = {clique};
        auto rowvector = matrix.getRowCoefficients(Cliques[clique]);
        auto indices = rowvector.getIndices();
        std::vector<int> bestNewClique = indices;
        for( int newClique = 0; newClique < maximalCliques.size(); ++newClique )
        {
            std::vector<int> coveredCliques;
            for( int potCovClique = 0; potCovClique < nCliques; ++potCovClique )
            {
                if( completedCliques[potCovClique] )
                    continue;
                bool covered = true;
                auto rowvec = matrix.getRowCoefficients(Cliques[potCovClique]);
                auto indices = rowvec.getIndices();
                for( int col = 0; col < rowvec.getLength(); ++col )
                {
                    int vertex = std::distance(Vertices.begin(), std::find(Vertices.begin(), Vertices.end(), indices[col]));
                    if( std::find(maximalCliques[newClique].begin(), maximalCliques[newClique].end(), vertex) == maximalCliques[newClique].end() )
                    {
                        covered = false;
                        break;
                    }
                }
                if( covered )
                    coveredCliques.push_back(potCovClique);
            }
            if( coveredCliques.size() > maxCoveredCliques.size() )
            {
                maxCoveredCliques = coveredCliques;
                bestNewClique = maximalCliques[newClique];
            }
        }
        if( maxCoveredCliques.size() > 1 )
        {
            TransactionGuard<REAL> tg{ reductions };
            result = PresolveStatus::kReduced;
            reductions.lockRow( Cliques[maxCoveredCliques[0]] );
            reductions.changeRowRHS( Cliques[maxCoveredCliques[0]], 1 );
            reductions.changeRowLHS( Cliques[maxCoveredCliques[0]], -std::numeric_limits<REAL>::infinity() );
            rowFlags[Cliques[maxCoveredCliques[0]]].set( RowFlag::kLhsInf);
            rowFlags[Cliques[maxCoveredCliques[0]]].unset( RowFlag::kRhsInf);
            for( int entryOfNewClique = 0; entryOfNewClique < bestNewClique.size(); ++entryOfNewClique )
            {
                reductions.lockColBounds( bestNewClique[entryOfNewClique] );
                if( problem.getLowerBounds()[bestNewClique[entryOfNewClique]] == 0 )
                {
                    reductions.changeMatrixEntry( Cliques[maxCoveredCliques[0]], bestNewClique[entryOfNewClique], 1 );
                }
                else
                {
                    reductions.changeMatrixEntry( Cliques[maxCoveredCliques[0]], bestNewClique[entryOfNewClique], -1 );
                }
            }
            for( int replacedClique = 1; replacedClique < maxCoveredCliques.size(); ++replacedClique )
            {
                completedCliques[maxCoveredCliques[replacedClique]] = true;
                reductions.lockRow( Cliques[maxCoveredCliques[replacedClique]] );
                reductions.markRowRedundant( Cliques[maxCoveredCliques[replacedClique]] );
            }
        }
    }
    return result;
}

} // namespace papilo

#endif