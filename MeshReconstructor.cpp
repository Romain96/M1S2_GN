/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <set>
#include <iostream>
#include <limits>

#include "MeshReconstructor.h"

#include "Octree.h"
#include "Plane.h"
#include "Graph.h"

// Eigen
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
using namespace Eigen;

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief MeshReconstructor::MeshReconstructor
 * @param vertices point cloud
 */
MeshReconstructor::MeshReconstructor(std::vector<Vertex *> &vertices) :
    _k(10),
    _iterations(2),
    _size(0.5f),
    _points(vertices),
    _p(0.f),
    _d(0.f),
    _subdivisionFactor(0.5f),
    _isolevel(0.f)
{
    // nothing
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief MeshReconstructor::getK
 * @return k the number of neighbours considered (user parameter)
 */
int MeshReconstructor::getK()
{
    return _k;
}

/**
 * @brief MeshReconstructor::getIterations
 * @return the number of iterations used to build the Octrees
 */
int MeshReconstructor::getIterations()
{
    return _iterations;
}

/**
 * @brief MeshReconstructor::getSize
 * @return the min size used to build the Octrees
 */
float MeshReconstructor::getSize()
{
    return _size;
}

/**
 * @brief MeshReconstructor::getDense
 * @return the coefficient of density
 */
float MeshReconstructor::getDense()
{
    return _d;
}

/**
 * @brief MeshReconstructor::getNoisy
 * @return the coefficient of noise
 */
float MeshReconstructor::getNoisy()
{
    return _d;
}

/**
 * @brief MeshReconstructor::getSubdivisionFactor
 * @return the subdivision factor (aka grid cell size)
 */
float MeshReconstructor::getSubdivisionFactor()
{
    return _subdivisionFactor;
}

/**
 * @brief MeshReconstructor::getIsolevel
 * @return the isolevel (aka level from which a point is considered inside/outside the isosurface)
 */
float MeshReconstructor::getIsolevel()
{
    return _isolevel;
}

/**
 * @brief MeshReconstructor::getPointTree
 * @return a pointer to the Octree based on the point cloud
 */
Octree *MeshReconstructor::getPointTree()
{
    return _pointTree;
}

/**
 * @brief MeshReconstructor::getCentroidTree
 * @return a pointer to the Octree based on the centroids
 */
Octree *MeshReconstructor::getCentroidTree()
{
    return _centroidTree;
}

/**
 * @brief MeshReconstructor::getCentroids
 * @return a vector containing all computed centroids
 */
std::vector<Vertex *>& MeshReconstructor::getCentroids()
{
    return _centroids;
}

/**
 * @brief MeshReconstructor::getPlanes
 * @return a vector containing all computed tangent planes
 */
std::vector<Plane *>& MeshReconstructor::getPlanes()
{
    return _planes;
}

/**
 * @brief MeshReconstructor::getComputedMesh
 * @return the computed Mesh
 */
Mesh& MeshReconstructor::getComputedMesh()
{
    return _result;
}

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief MeshReconstructor::setK
 * @param k number of neighbours to consider (k nearest neighbours search and centroid computation)
 */
void MeshReconstructor::setK(int k)
{
    _k = k;
}

/**
 * @brief MeshReconstructor::setIterations
 * @param it number of iterations too build the Octrees
 */
void MeshReconstructor::setIterations(int it)
{
    _iterations = it;
}

/**
 * @brief MeshReconstructor::setSize
 * @param size min size to reach when building the Octrees
 */
void MeshReconstructor::setSize(float size)
{
    _size = size;
}

/**
 * @brief MeshReconstructor::setDense
 * @param p new codfficient of density
 */
void MeshReconstructor::setDense(float p)
{
    _p = p;
}

/**
 * @brief MeshReconstructor::setNoisy
 * @param d new coefficient of noise
 */
void MeshReconstructor::setNoisy(float d)
{
    _d = d;
}

/**
 * @brief MeshReconstructor::setSubdivisionFactor
 * @param factor new subdivision faction
 */
void MeshReconstructor::setSubdivisionFactor(float factor)
{
    _subdivisionFactor = factor;
}

/**
 * @brief MeshReconstructor::setIsolevel
 * @param level new isolevel
 */
void MeshReconstructor::setIsolevel(float level)
{
    _isolevel = level;
}

/**
 * @brief MeshReconstructor::setPointTree
 * @param t new Octree based on the point cloud
 */
void MeshReconstructor::setPointTree(Octree *t)
{
    _pointTree = t;
}

/**
 * @brief MeshReconstructor::setCentroidTree
 * @param t new Octree based on the centroids
 */
void MeshReconstructor::setCentroidTree(Octree *t)
{
    _centroidTree = t;
}

/**
 * @brief MeshReconstructor::setCentroids
 * @param centroids new list of computed centroids
 */
void MeshReconstructor::setCentroids(std::vector<Vertex *> &centroids)
{
    _centroids = centroids;
}

/**
 * @brief MeshReconstructor::setPlanes
 * @param planes new list of computed tangent planes
 */
void MeshReconstructor::setPlanes(std::vector<Plane *> &planes)
{
    _planes = planes;
}

/**
 * @brief MeshReconstructor::setComputedMesh
 * @param m new computed Mesh
 */
void MeshReconstructor::setComputedMesh(Mesh &m)
{
    _result = m;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------

/**
 * @brief MeshReconstructor::buildPointTree
 */
void MeshReconstructor::buildPointTreeWithIterations()
{
    _pointTree = new Octree();
    _pointTree->findSpaceBorders(_points);
    _pointTree->constructWithIterations(_iterations, _points);
}

/**
 * @brief MeshReconstructor::buildPointTreeWithSize
 */
void MeshReconstructor::buildPointTreeWithSize()
{
    _pointTree = new Octree();
    _pointTree->findSpaceBorders(_points);
    _pointTree->constructWithMinSize(_size, _points);
}

/**
 * @brief MeshReconstructor::computePlanes
 */
void MeshReconstructor::computeCentroidsAndTangentPlanes()
{
    std::cout << "computing centroids" << std::endl;

    // simply calling Octree::findKNearestNeighbours for each point in _mesh
    std::vector<Vertex *>::iterator pointIterator;
    std::vector<std::pair<Vertex *, float>> neighbours;
    int id = 0;
    glm::vec3 centroid(0.f);
    glm::mat3x3 covarianceMatrix(0.f);
    Plane *p;
    Matrix3f eigen;
    Matrix3f ev;

    for (pointIterator = _points.begin(); pointIterator < _points.end(); pointIterator++)
    {
        // clearing previously retrieved neighbours
        neighbours.clear();

        // retrieving the k nearest neighbours of ref
        neighbours = _pointTree->findKNeartestNeighbours((*pointIterator), _k);

        // computing the centroid
        for (unsigned int i = 0; i < neighbours.size(); i++)
        {
            centroid = centroid + neighbours[i].first->getPosition();
        }

        // storing the centroid in _centroids
        centroid = centroid / (float)_k;
        _centroids.push_back(new Vertex(id++, centroid.x, centroid.y, centroid.z));
        //std::cout << "point : " << (*pointIterator)->getPosition().x << ", " << (*pointIterator)->getPosition().y << ", " << (*pointIterator)->getPosition().z << std::endl;
        //std::cout << "centroid : " << centroid.x << ", " << centroid.y << ", " << centroid.z << std::endl;

        // computing the covariance matrix
        for (unsigned int i = 0; i < neighbours.size(); i++)
        {
            // outer product to form the covariance matrix !
            glm::mat3x3 cv;
            glm::vec3 v = neighbours[i].first->getPosition() - centroid;
            cv[0].x = v.x * v.x;
            cv[0].y = v.x * v.y;
            cv[0].z = v.x * v.z;

            cv[1].x = v.y * v.x;
            cv[1].y = v.y * v.y;
            cv[1].z = v.y * v.z;

            cv[2].x = v.z * v.x;
            cv[2].y = v.z * v.y;
            cv[2].z = v.z * v.z;

            covarianceMatrix = covarianceMatrix + cv;
        }

        // using Eigen lib and Principal Component Analysis to find the normal and local plane in Oi (centroid)
        eigen << covarianceMatrix[0].x, covarianceMatrix[0].y, covarianceMatrix[0].z,
                covarianceMatrix[1].x, covarianceMatrix[1].y, covarianceMatrix[1].z,
                covarianceMatrix[2].x, covarianceMatrix[2].y, covarianceMatrix[2].z;

        SelfAdjointEigenSolver<Matrix3f> eigensolver(eigen);
        //std::cout << "eigenvectors : " << eigensolver.eigenvectors() << std::endl;
        ev = eigensolver.eigenvectors();

        // storing the eigenvectors in a Plane
        glm::vec3 temp;
        p = new Plane();
        temp = glm::vec3(ev(0,0), ev(0,1), ev(0,2));
        p->setEigenvector1(temp);
        temp = glm::vec3(ev(1,0), ev(1,1), ev(1,2));
        p->setEigenvector2(temp);
        temp = glm::vec3(ev(2,0), ev(2,1), ev(2,2));
        p->setEigenvector3(temp);
        _planes.push_back(p);
    }
}

/**
 * @brief MeshReconstructor::buildCentroidTreeWithIterations
 */
void MeshReconstructor::buildCentroidTreeWithIterations()
{
    _centroidTree = new Octree();
    _centroidTree->findSpaceBorders(_centroids);
    _centroidTree->constructWithIterations(_iterations, _centroids);
}

/**
 * @brief MeshReconstructor::buildCentroidTreeWithSize
 */
void MeshReconstructor::buildCentroidTreeWithSize()
{
    _centroidTree = new Octree();
    _centroidTree->findSpaceBorders(_centroids);
    _centroidTree->constructWithMinSize(_size, _centroids);
}

/**
 * @brief MeshReconstructor::reorientateTangentPlanes
 */
void MeshReconstructor::reorientateTangentPlanes()
{
    // 5 steps are necessary
    // 1) compute a full connected graph (distance between centroids)
    // 2) compute the Euclidian Minimum Spanning Tree of this graph
    // 3) enrich the EMST with the "degree of closeness of tangent planes"
    //    by remplacing the cost of present edges
    //    and adding those of centroid neighbourhood
    //    This is the Riemannian Graph
    // 4) compute the Minimum Spanning Tree of this Riemannian Graph
    // 5) traverse this MST and reorient the planes based on the cost of edges

    // Step 1 : Full connected graph (Euclidian distance on centroids)
    std::cout << "building a full Euclidian connected graph..." << std::endl;
    _graph->buildEuclidianGraph(_centroids, _planes);
    std::cout << "done" << std::endl;

    // Step 2 : Euclidian Minimum Spanning Tree
    std::cout << "building Euclidian Minimum Spanning Tree..." << std::endl;
    Graph *emst = _graph->buildMinimumSpanningTree();
    _graph->clearNodes();
    _graph->clearEdges();
    _graph = emst;
    std::cout << "done" << std::endl;

    // Step 3 : Riemannian Graph
    std::cout << "building Riemannian graph..." << std::endl;
    _graph->enhanceToRiemannianGraph(_k, _centroidTree);
    std::cout << "done" << std::endl;

    // Step 4 : Minimum Spanning Tree
    std::cout << "building a Minimum Spanning Tree..." << std::endl;
    Graph *mst = _graph->buildMinimumSpanningTree();
    _graph->clearNodes();
    _graph->clearEdges();
    _graph = mst;
    std::cout << "done" << std::endl;

    // Step 5 : reorienting planes
    std::cout << "traversing the MST and reorientating planes..." << std::endl;
    _graph->traverseDepthFirstAndReorientate();
    std::cout << "done" << std::endl;
}

/**
 * @brief MeshReconstructor::__findNearestTangentPlaneAsCentroid
 * @param p point of reference
 * @return a pointer to the closest tangent plane from p
 */
Vertex *MeshReconstructor::__findNearestTangentPlaneAsCentroid(Vertex *p)
{
    // we consider that the closest tangent plane from p is the tangent plane
    // whom associated centroid is the closest to p
    // thus using the centroid octree to quicken the search

    // first searching in which leaf of the centroid octree is p
    Octree *t = _centroidTree;
    while (!t->leaf())
    {
        if (Octree::isInside(p, t->getLowerNE()))
        {
            t = t->getLowerNE();
        }
        else if (Octree::isInside(p, t->getLowerNW()))
        {
            t = t->getLowerNW();
        }
        else if (Octree::isInside(p, t->getLowerSE()))
        {
            t = t->getLowerSE();
        }
        else if (Octree::isInside(p, t->getLowerSW()))
        {
            t = t->getLowerSW();
        }
        else if (Octree::isInside(p, t->getUpperNE()))
        {
            t = t->getUpperNE();
        }
        else if (Octree::isInside(p, t->getUpperNW()))
        {
            t = t->getUpperNW();
        }
        else if (Octree::isInside(p, t->getUpperSE()))
        {
            t = t->getUpperSE();
        }
        else if (Octree::isInside(p, t->getUpperSW()))
        {
            t = t->getUpperSW();
        }
        else
        {
            std::cerr << "MeshReconstructor::__findNearestPlane error point is not inside one of the 8 children !" << std::endl;
            break;
        }
    }

    // finding the closest centroid in this leaf
    float dist;
    float distMin = 10000.f;
    Vertex *closestCentroid;
    std::vector<Vertex *>::iterator it;

    for (it = t->getPoints().begin(); it != t->getPoints().end(); it++)
    {
        dist = Vertex::distance3(p->getPosition(), (*it)->getPosition());
        if (dist < distMin)
        {
            distMin = dist;
            closestCentroid = (*it);
        }
    }

    return closestCentroid;
}

/**
 * @brief MeshReconstructor::__signedDistanceToClosestTangentPlane
 * @param p point of reference
 * @return the signed distance to the closest tangent plane from p
 */
float MeshReconstructor::__signedDistanceToClosestTangentPlane(Vertex *p)
{
    // retrieving the closest tangent plane with __findNearestPlane
    Vertex *centroid = __findNearestTangentPlaneAsCentroid(p);
    // and its associated tangent plane
    Plane *tp = _planes[centroid->getId()];

    // computing z the projection of p into it's closest tangent plane
    // z <- oi - (((p - oi).ni) * ni)
    glm::vec3 z = centroid->getPosition() - ((glm::dot((p->getPosition() - centroid->getPosition()), tp->getEigenvector3())) * tp->getEigenvector3());

    // if distance between z and oi is less than p + d we return the distance
    // else we set the distance as undefined using infinity to represent "undefined" state
    float dist = Vertex::distance3(z, centroid->getPosition());
    if (dist < _p + _d)
        return glm::dot((p->getPosition() - centroid->getPosition()), tp->getEigenvector3());
    else
        return std::numeric_limits<float>::infinity();
}


/**
 * @brief MeshReconstructor::polygonise
 * @param grid
 * @param triangles
 * @return the number of triangles created
 * table and function adapted from Paul Bourke
 * http://paulbourke.net/geometry/polygonise/
 */
int MeshReconstructor::polygonise(GRIDCELL grid, TRIANGLE *triangles)
{
    int i;
    int ntriang;
    int cubeindex;
    Vertex vertlist[12];

    int edgeTable[256]={
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };


    int triTable[256][16] =
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};


    /*
     * Determine the index into the edge table which
     * tells us which vertices are inside of the surface
    */

    cubeindex = 0;
    if (grid.val[0] < _isolevel) cubeindex |= 1;
    if (grid.val[1] < _isolevel) cubeindex |= 2;
    if (grid.val[2] < _isolevel) cubeindex |= 4;
    if (grid.val[3] < _isolevel) cubeindex |= 8;
    if (grid.val[4] < _isolevel) cubeindex |= 16;
    if (grid.val[5] < _isolevel) cubeindex |= 32;
    if (grid.val[6] < _isolevel) cubeindex |= 64;
    if (grid.val[7] < _isolevel) cubeindex |= 128;

    /* Cube is entirely in/out of the surface */
    if (edgeTable[cubeindex] == 0)
        return(0);

    /* Find the vertices where the surface intersects the cube */
    if (edgeTable[cubeindex] & 1)
        vertlist[0] = vertexInterpolate(grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
    if (edgeTable[cubeindex] & 2)
        vertlist[1] = vertexInterpolate(grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
    if (edgeTable[cubeindex] & 4)
          vertlist[2] = vertexInterpolate(grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
    if (edgeTable[cubeindex] & 8)
          vertlist[3] = vertexInterpolate(grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
    if (edgeTable[cubeindex] & 16)
          vertlist[4] = vertexInterpolate(grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
    if (edgeTable[cubeindex] & 32)
          vertlist[5] = vertexInterpolate(grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
    if (edgeTable[cubeindex] & 64)
          vertlist[6] = vertexInterpolate(grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
    if (edgeTable[cubeindex] & 128)
          vertlist[7] = vertexInterpolate(grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
    if (edgeTable[cubeindex] & 256)
          vertlist[8] = vertexInterpolate(grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
    if (edgeTable[cubeindex] & 512)
          vertlist[9] = vertexInterpolate(grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
    if (edgeTable[cubeindex] & 1024)
          vertlist[10] = vertexInterpolate(grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
    if (edgeTable[cubeindex] & 2048)
          vertlist[11] = vertexInterpolate(grid.p[3],grid.p[7],grid.val[3],grid.val[7]);

    /* Create the triangle */
    ntriang = 0;
    for (i = 0; triTable[cubeindex][i] != -1; i += 3)
    {
         triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i  ]];
         triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i+1]];
         triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i+2]];
         ntriang++;
    }

    return(ntriang);
}

/**
 * @brief MeshReconstructor::vertexInterpolate
 * @param p1 vertex position 1
 * @param p2 vertex position 2
 * @param valp1 value associated with vertex 1
 * @param valp2 value associated with vertex 2
 * @return a new vertex being the linear interpolation of p1 and p2
 * function adapted from Paul Bourke
 * http://paulbourke.net/geometry/polygonise/
 */
Vertex MeshReconstructor::vertexInterpolate(Vertex p1, Vertex p2, double valp1, double valp2)
{
    double mu;
    Vertex p;

    if (fabs(_isolevel - valp1) < 0.00001)
        return(p1);
    if (fabs(_isolevel - valp2) < 0.00001)
        return(p2);
    if (fabs(valp1 - valp2) < 0.00001)
        return(p1);

    mu = (_isolevel - valp1) / (valp2 - valp1);

    float x = p1.getX() + mu * (p2.getX() - p1.getX());
    float y = p1.getY() + mu * (p2.getY() - p1.getY());
    float z = p1.getZ() + mu * (p2.getZ() - p1.getZ());
    p.setX(x);
    p.setY(y);
    p.setZ(z);

    return(p);
}

/**
 * @brief MeshReconstructor::createIsosurface
 */
void MeshReconstructor::createIsosurface()
{
    // convention of vertex index used is the following  :
    // 0 : lower south east vertex
    // 1 : lower north east vertex
    // 2 : lower north west vertex
    // 3 : lower south west vertex
    // 4 : upper south east vertex
    // 5 : upper north east vertex
    // 6 : upper north west vertex
    // 7 : upper south west vertex

    glm::vec3 lowerEnglobingVertex = _centroidTree->getBorderLowerSW();
    glm::vec3 upperEnglobingVertex = _centroidTree->getBorderUpperNE();

    std::cout << "Marching cubes algorithm..." << std::endl;

    GRIDCELL cell;
    float minX = lowerEnglobingVertex.x;
    float minY = lowerEnglobingVertex.y;
    float minZ = lowerEnglobingVertex.z;
    float maxX = lowerEnglobingVertex.x + _subdivisionFactor;
    float maxY = lowerEnglobingVertex.y + _subdivisionFactor;
    float maxZ = lowerEnglobingVertex.z + _subdivisionFactor;

    // for each cube in the surface englobing cube
    for (float x = lowerEnglobingVertex.x + _subdivisionFactor; x < upperEnglobingVertex.x; x += _subdivisionFactor)
    {
        for (float y = lowerEnglobingVertex.y + _subdivisionFactor; y < upperEnglobingVertex.y; y += _subdivisionFactor)
        {
            for (float z = lowerEnglobingVertex.z + _subdivisionFactor; z < upperEnglobingVertex.z; z += _subdivisionFactor)
            {
                // vertices
                Vertex *v0 = new Vertex(0, minX, maxY, minZ);
                Vertex *v1 = new Vertex(0, maxX, maxY, minZ);
                Vertex *v2 = new Vertex(0, maxX, minY, minZ);
                Vertex *v3 = new Vertex(0, minX, minY, minZ);
                Vertex *v4 = new Vertex(0, minX, maxY, maxZ);
                Vertex *v5 = new Vertex(0, maxX, maxY, maxZ);
                Vertex *v6 = new Vertex(0, maxX, minY, maxZ);
                Vertex *v7 = new Vertex(0, minX, maxY, maxZ);

                // filling the position of each 8 vertices of the cube
                cell.p[0] = Vertex(0, minX, maxY, minZ);
                cell.p[1] = Vertex(0, maxX, maxY, minZ);
                cell.p[2] = Vertex(0, maxX, minY, minZ);
                cell.p[3] = Vertex(0, minX, minY, minZ);
                cell.p[4] = Vertex(0, minX, maxY, maxZ);
                cell.p[5] = Vertex(0, maxX, maxY, maxZ);
                cell.p[6] = Vertex(0, maxX, minY, maxZ);
                cell.p[7] = Vertex(0, minX, maxY, maxZ);

                // filling the value of each 8 vertices of the cube
                cell.val[0] = MeshReconstructor::__signedDistanceToClosestTangentPlane(v0);
                cell.val[1] = MeshReconstructor::__signedDistanceToClosestTangentPlane(v1);
                cell.val[2] = MeshReconstructor::__signedDistanceToClosestTangentPlane(v2);
                cell.val[3] = MeshReconstructor::__signedDistanceToClosestTangentPlane(v3);
                cell.val[4] = MeshReconstructor::__signedDistanceToClosestTangentPlane(v4);
                cell.val[5] = MeshReconstructor::__signedDistanceToClosestTangentPlane(v5);
                cell.val[6] = MeshReconstructor::__signedDistanceToClosestTangentPlane(v6);
                cell.val[7] = MeshReconstructor::__signedDistanceToClosestTangentPlane(v7);

                TRIANGLE res[5];
                int n = polygonise(cell, res);
                std::cout << n << " triangles created" << std::endl;

                // next cube position
                minX += _subdivisionFactor;
                minY += _subdivisionFactor;
                minZ += _subdivisionFactor;
                maxX += _subdivisionFactor;
                maxY += _subdivisionFactor;
                maxZ += _subdivisionFactor;
            }
        }
    }
}
