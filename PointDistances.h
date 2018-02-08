#ifndef POINTDISTANCES_H
#define POINTDISTANCES_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <map>
#include <vector>

#include "Vertex.h"

/**
 * @brief The PointDistances class
 */
class PointDistances
{
protected:
    std::map<int, std::pair<int, std::vector<Vertex *>>> _dist;

public:
    // constructor(s)
    PointDistances(std::vector<Vertex *>& vertexList);

    // getter(s)

    // setter(s)

    // method(s)
    void constructDistanceTable(std::vector<Vertex *>& vertexList);
    void distance(int i, int j);
}

#endif // POINTDISTANCES_H
