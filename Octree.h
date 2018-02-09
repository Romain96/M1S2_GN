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

    // max number of iterations when creating the Octree (with iterations method)
    static int _depth;
    static bool _depthTest;
    // min size of regions when creating the Octree (with min size method)
    static float _minSize;
    static bool _sizeTest;

    // points included in a leaf
    std::vector<Vertex *>_points;

public:
    // constructor(s)
    Octree();   // to create the root / initialize the class
    // better suited to instanciate children (whose boundaries are known through the father's ones)
    Octree(glm::vec3& lowerNW, glm::vec3& lowerNE, glm::vec3& lowerSW, glm::vec3& lowerSE,
           glm::vec3& upperNW, glm::vec3& upperNE, glm::vec3& upperSW, glm::vec3& upperSE);

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
    void setBorderLowerNW(glm::vec3& v);
    void setBorderLowerNE(glm::vec3& v);
    void setBorderLowerSW(glm::vec3& v);
    void setBorderLowerSE(glm::vec3& v);
    void setBorderUpperNW(glm::vec3& v);
    void setBorderUpperNE(glm::vec3& v);
    void setBorderUpperSW(glm::vec3& v);
    void setBorderUpperSE(glm::vec3& v);

    void setLowerNW(Octree *t);
    void setLowerNE(Octree *t);
    void setLowerSW(Octree *t);
    void setLowerSE(Octree *t);
    void setUpperNW(Octree *t);
    void setUpperNE(Octree *t);
    void setUpperSW(Octree *t);
    void setUpperSE(Octree *t);

    // method(s)
    bool leaf();
    void setLeaf();
    void setNotLeaf();
    static void deleteFromNode(Octree *t);
    void findSpaceBorders(std::vector<Vertex *>& vertices);
    void constructWithMinSize(float size, std::vector<Vertex *>& vertices);
    void constructWithIterations(int k, std::vector<Vertex *>& vertices);
    std::vector<Vertex *>& findKNeartestNeighbours(Vertex *ref, int k);

private:
    // internal method(s)
    void addPoint(Vertex *v);
    static void __buildOctreeNode(Octree *t, int depth, std::vector<Vertex *>& vertices);
    static void __findPointsInRegion(Octree *t, std::vector<Vertex *>& vertives);
};

#endif // OCTREE_H
