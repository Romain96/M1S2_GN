#ifndef MESHRECONSTRUCTOR_H
#define MESHRECONSTRUCTOR_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// glm
#include "glm/glm.hpp"
#include "glm/vec3.hpp"

#include <vector>

// Eigen
#include "Eigen/Dense"

#include "Vertex.h"
#include "Mesh.h"
#include "Octree.h"

using namespace Eigen;

/**
 * @brief The MeshReconstructor class
 */
class MeshReconstructor
{
private:
    int _k;

protected:
    Mesh _mesh;
    Octree *_tree;
    std::vector<glm::vec3> _centroids;
    std::vector<glm::vec3> _normals;
    // TODO tangentPlanes

public:
    // constructor(s)
    MeshReconstructor(Mesh m);

    // getter(s)
    int getK();

    // setter(s)
    void setK(int k);

    // method(s)
    void computePlanes();
};

#endif // MESHRECONSTRUCTOR_H
