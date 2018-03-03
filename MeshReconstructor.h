#ifndef MESHRECONSTRUCTOR_H
#define MESHRECONSTRUCTOR_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// STD
#include <vector>
#include <cmath>

// GLM
#include "glm/glm.hpp"
#include "glm/vec3.hpp"

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
 * @brief The MeshReconstructor class
 */
class MeshReconstructor
{
private:
    // k neighbourhood
    int _k;
    // number of iterations to create the Octree
    int _iterations;
    // minimum size to reach when creating the Octree
    float _size;
    // surface is supposed to be p-dense
    float _p;
    // surface is supposed to be d-noisy
    float _d;
    // subdivision factor (aka grid cell size for marching cubes algorithm)
    float _subdivisionFactor;
    // isolevel for marching cubes (aka level from which a point is considered inside/outside the isosurface)
    double _isolevel;

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
    Graph *_graph;
    // computed Mesh
    Mesh _result;

public:
    // constructor(s)
    MeshReconstructor(std::vector<Vertex *>& vertices);

    // getter(s)
    int getK();
    int getIterations();
    float getSize();
    float getDense();
    float getNoisy();
    float getSubdivisionFactor();
    float getIsolevel();
    Octree *getPointTree();
    Octree *getCentroidTree();
    std::vector<Vertex *>& getCentroids();
    std::vector<Plane *>& getPlanes();
    Mesh& getComputedMesh();

    // setter(s)
    void setK(int k);
    void setIterations(int it);
    void setSize(float size);
    void setDense(float p);
    void setNoisy(float d);
    void setSubdivisionFactor(float factor);
    void setIsolevel(float level);
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
    void reorientateTangentPlanes();
    // marching cubes
    int polygonise(GRIDCELL grid, TRIANGLE *triangles);
    Vertex vertexInterpolate(Vertex p1, Vertex p2, double valp1, double valp2);
    void createIsosurface();


    Vertex *__findNearestTangentPlaneAsCentroid(Vertex *p);
    float __signedDistanceToClosestTangentPlane(Vertex *p);
};

#endif // MESHRECONSTRUCTOR_H
