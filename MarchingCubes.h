#ifndef MARCHINGCUBES_H
#define MARCHINGCUBES_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

class MarchingCubes
{
protected:
    unsigned float _subdivisionFactor;

public:
    // constructor(s)

    // getter(s)
    unsigned float getSubdivisionFactor();

    // settet(s)
    void setSubdivisionFactor(unsigned float factor);

    // method(s)
};

#endif // MARCHINGCUBES_H
