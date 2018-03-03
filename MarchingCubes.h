#ifndef MARCHINGCUBES_H
#define MARCHINGCUBES_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// STD
#include <cmath>

#include "Vertex.h"

// triangle structure used to represent a triangle (by representing 3 Vertices)
typedef struct s_triangle
{
    Vertex p[3];
} TRIANGLE;

// structure used to represent a gridcell (8 vertices + 8 associated values)
typedef struct s_gridcell
{
    Vertex p[8];
    double val[8];
} GRIDCELL;

/**
 * @brief The MarchingCubes class
 */
class MarchingCubes
{
protected:
    float _subdivisionFactor;
    double _isolevel;

public:
    // constructor(s)
    MarchingCubes();
    MarchingCubes(float factor, double isolevel);

    // getter(s)
    float getSubdivisionFactor();
    double getIsolevel();

    // settet(s)
    void setSubdivisionFactor(float factor);
    void setIsolevel(double isolevel);

    // method(s)
    int polygonise(GRIDCELL grid, TRIANGLE *triangles);
    Vertex vertexInterpolate(Vertex p1, Vertex p2, double valp1, double valp2);
};

#endif // MARCHINGCUBES_H
