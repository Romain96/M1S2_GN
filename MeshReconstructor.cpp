/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <set>
#include <iostream>

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
    _points(vertices)
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
 * @brief MeshReconstructor::buildGraph
 */
void MeshReconstructor::buildGraph()
{
    _graph = Graph();
    //_graph.buildGraph(_k, _centroidTree, _centroids, _planes);
    _graph.buildGraphFull(_k, _centroidTree, _centroids, _planes);
}

/**
 * @brief MeshReconstructor::reorientateTangentPlanes
 */
void MeshReconstructor::reorientateTangentPlanes()
{
    // first compute the MST from the graph
    std::cout << "building MST" << std::endl;
    _mst = _graph.buildMinimumSpanningTree();
    std::cout << "MST built" << std::endl;

    // then reorientate the planes
}

/**
 * @brief MeshReconstructor::__findNearestPlane
 * @param p point of reference
 * @return a pointer to the closest tangent plane from p
 */
Plane *MeshReconstructor::__findNearestPlane(Vertex *p)
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

    // we return the Plane at index centroid->id since
    // Planes and Centroids are indexed in the same order
    return _planes[closestCentroid->getId()];
}
