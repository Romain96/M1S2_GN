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
#include "Plane.h"
#include "Node.h"
#include "Edge.h"
#include "Graph.h"

using namespace Eigen;

/**
 * @brief The MeshReconstructor class
 */
class MeshReconstructor
{
private:
    int _k;
    int _iterations;
    float _size;

protected:
    // list of all vertices (point cloud)
    std::vector<Vertex *> _points;
    // Octree on the point cloud used to compute centroids
    Octree *_pointTree;
    // Octree on the centroids used to determine planes orientation
    Octree *_centroidTree;
    // list of all computed centroids
    std::vector<Vertex *> _centroids;
    // list of all computed tangent planes with their associated normal
    std::vector<Plane *> _planes;
    // Graph used to reorientate the tangent planes
    Graph _graph;
    // Minimum Spanning Tree (MST) derived from _graph and used to reorientate tangent planes
    Graph *_mst;
    // computed Mesh
    Mesh _result;

public:
    // constructor(s)
    MeshReconstructor(std::vector<Vertex *>& vertices);

    // getter(s)
    int getK();
    int getIterations();
    float getSize();
    Octree *getPointTree();
    Octree *getCentroidTree();
    std::vector<Vertex *>& getCentroids();
    std::vector<Plane *>& getPlanes();
    Mesh& getComputedMesh();

    // setter(s)
    void setK(int k);
    void setIterations(int it);
    void setSize(float size);
    void setPointTree(Octree *t);
    void setCentroidTree(Octree *t);
    void setCentroids(std::vector<Vertex *>& centroids);
    void setPlanes(std::vector<Plane *>& planes);
    void setComputedMesh(Mesh& m);

    // method(s)
    void buildPointTreeWithIterations();
    void buildPointTreeWithSize();
    void computeCentroidsAndTangentPlanes();
    void buildCentroidTreeWithIterations();
    void buildCentroidTreeWithSize();
    void buildGraph();
    void reorientateTangentPlanes();

private:
    void __buildCentroidOctree(std::vector<glm::vec3>& centroids);
};

#endif // MESHRECONSTRUCTOR_H
