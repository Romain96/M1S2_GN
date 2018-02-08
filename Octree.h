#ifndef OCTREE_H
#define OCTREE_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// glm
#include "glm/glm.hpp"
#include "glm/vec3.hpp"

// stl
#include <iostream>
#include <vector>

#include "Vertex.h"

/**
 * @brief The Octree class
 */
class Octree
{
protected:
    // the 8 most external points constrain the current node 3D space
    // they formed a cube
    glm::vec3 _borderLowerNW;
    glm::vec3 _borderLowerNE;
    glm::vec3 _borderLowerSW;
    glm::vec3 _borderLowerSE;
    glm::vec3 _borderUpperNW;
    glm::vec3 _borderUpperNE;
    glm::vec3 _borderUpperSW;
    glm::vec3 _borderUpperSE;

    // halting point (equivalent to "each child is null")
    bool _isLeaf;

    // each node has 8 children
    // the space is subdivided in 8 subregions
    Octree *_lowerNW;
    Octree *_lowerNE;
    Octree *_lowerSW;
    Octree *_lowerSE;
    Octree *_upperNW;
    Octree *_upperNE;
    Octree *_upperSW;
    Octree *_upperSE;

public:
    // constructor(s)
    Octree();

    // getter(s)
    glm::vec3& getBorderLowerNW();
    glm::vec3& getBorderLowerNE();
    glm::vec3& getBorderLowerSW();
    glm::vec3& getBorderLowerSE();
    glm::vec3& getBorderUpperNW();
    glm::vec3& getBorderUpperNE();
    glm::vec3& getBorderUpperSW();
    glm::vec3& getBorderUpperSE();

    Octree *getLowerNW();
    Octree *getLowerNE();
    Octree *getLowerSW();
    Octree *getLowerSE();
    Octree *getUpperNW();
    Octree *getUpperNE();
    Octree *getUpperSW();
    Octree *getUpperSE();

    // setter(s)
    void getBorderLowerNW(glm::vec3& v);
    void getBorderLowerNE(glm::vec3& v);
    void getBorderLowerSW(glm::vec3& v);
    void getBorderLowerSE(glm::vec3& v);
    void getBorderUpperNW(glm::vec3& v);
    void getBorderUpperNE(glm::vec3& v);
    void getBorderUpperSW(glm::vec3& v);
    void getBorderUpperSE(glm::vec3& v);

    void getLowerNW(Octree *t);
    void getLowerNE(Octree *t);
    void getLowerSW(Octree *t);
    void getLowerSE(Octree *t);
    void getUpperNW(Octree *t);
    void getUpperNE(Octree *t);
    void getUpperSW(Octree *t);
    void getUpperSE(Octree *t);

    // method(s)
    bool leaf();
    void findSpaceBorders();
    void constructWithMinSize(float size);
    void constructWithIterations(int k);
    std::vector<Vertex *>& findKNeartestNeighbours(Vertex *ref, int k);
};

#endif // OCTREE_H
