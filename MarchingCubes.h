#ifndef MARCHINGCUBES_H
#define MARCHINGCUBES_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

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
    unsigned float _subdivisionFactor;
    double _isolevel;

public:
    // constructor(s)

    // getter(s)
    unsigned float getSubdivisionFactor();

    // settet(s)
    void setSubdivisionFactor(unsigned float factor);

    // method(s)
    int polygonise(GRIDCELL grid, double isolevel, TRIANGLE *triangles);
    Vertex vertexInterpolate(double isolevel, Vertex p1, Vertex p2, double valp1, double valp2);
};

#endif // MARCHINGCUBES_H
