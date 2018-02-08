#ifndef MESHRECONSTRUCTOR_H
#define MESHRECONSTRUCTOR_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include "glm/glm.hpp"
#include "glm/vec3.hpp"

#include <vector>

#include "Vertex.h"
#include "Mesh.h"
#include "PointDistances.h"

/**
 * @brief The MeshReconstructor class
 */
class MeshReconstructor
{
private:
    static int _k;

protected:
    Mesh _mesh;
    PointDistances _distances;
    std::vector<glm::vec3> _centroids;
    std::vector<glm::vec3> _normals;
    // TODO tangentPlanes

public:
    // constructor(s)
    MeshReconstructor(Mesh m);

    // getter(s)

    // setter(s)

    // method(s)
    void computeCentroids();
    void computeNormalsPCA();
    void computeTangentPlanes();
};

#endif // MESHRECONSTRUCTOR_H
