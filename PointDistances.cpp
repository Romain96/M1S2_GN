/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <map>
#include <vector>
#include <iostream>

#include "Vertex.h"
#include "PointDistances.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

/**
 * @brief PointDistances::PointDistances
 */
PointDistances::PointDistances()
{
    // nothing
}

/**
 * @brief PointDistances::PointDistances builds the distance table
 * @param vertexList list of all vertices (point cloud)
 */
PointDistances::PointDistances(std::vector<Vertex *> &vertexList)
{
    std::cout << "building distance table" << std::endl;

    // building the table
    constructDistanceTable(vertexList);

    std::cout << "distance table built" << std::endl;
}


//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------

/**
 * @brief PointDistances::setVertices
 * @param v list of all vertices (point cloud)
 */
void PointDistances::setVertices(std::vector<Vertex *> &v)
{
    std::cout << "building distance table" << std::endl;

    constructDistanceTable(v);

    std::cout << "distance table built" << std::endl;
}

/**
 * @brief PointDistances::constructDistanceTable
 * @param vertexList
 */
void PointDistances::constructDistanceTable(std::vector<Vertex *> &vertexList)
{
    // Table is a triangle (upper triangle = lower triangle and diagonal = 0)
    // distance(i, j) = distance(j, i)
    // distance(i, i) = 0
    // we will store only the lower triangle

    int vertices = vertexList.size();
    int n = 1;
    Vertex *vRef = nullptr;     // reference vector (we compare each other to him)
    Vertex *vComp = nullptr;    // comparaison vector

    // for each vertex in vertexList (except 0)
    for (int i = 1; i < vertices - 1; i++)
    {
        // retrieving the reference vector
        vRef = vertexList[i];

        // computing n distances from (current_vertex vertex_0) to (current_vertex vertex_n)
        for (int j = 0; j < n; j++)
        {
            // vertexList is supposed sorted by increasing ID
            // TODO sort vertexList before

            // retrieving vComp
            vComp = vertexList[j];

            // retrieving vector
            std::vector<float> v = *(&_dist[i]);

            // computing the Euclidian distance
            float dist = Vertex::distance3(vRef->getPosition(), vComp->getPosition());

            v.push_back(dist);
        }

        // lower line has n+1 element
        n++;
    }
}

/**
 * @brief PointDistances::distance
 * @param i id of the first Vertex
 * @param j id of the second Vertex
 * @return the Euclidian distance between Vertex of id i and Vertex of id j
 */
float PointDistances::distance(int i, int j)
{
    if(i < 0 || j < 0)
    {
        std::cerr << "PointDistances::distance error i or j < 0" << std::endl;
        return 0.f;
    }
    else if (i >= _dist.size() || j >= _dist.size())
    {
        std::cerr << "PointDistances::distance error i or j > nb vertices" << std::endl;
        return 0.f;
    }
    else
    {
        // Dist(i,j) = (_dist[i])[j] if j < (i-1) or dis(j,i) otherwise
        if (j < (i - 1))
        {
            std::map<int, std::vector<float>>::iterator mapIt = _dist.begin();
            std::advance(mapIt, i);
            std::vector<float> f = (*mapIt).second;
            std::vector<float>::iterator fIt = f.begin();
            std::advance(fIt, j);
            return (*fIt);
        }
        else
        {
            return distance(j, i);
        }
    }
}
