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
 * @param m mesh (point cloud)
 */
MeshReconstructor::MeshReconstructor(Mesh m) :
    _pointTree(nullptr),
    _centroidTree(nullptr),
    _centroids(),
    _planes()
{
    // calling a method to build the pointTree
    _pointTree = new Octree();
    _pointTree->findSpaceBorders(m.getVertices());
    _pointTree->constructWithIterations(2, m.getVertices());

    // calling a method to build the centroids and tangent planes
    computeCentroidsAndTangentPlanes(m);

    // building the centroidTree
    _centroidTree = new Octree();
    _centroidTree->findSpaceBorders(_centroids);
    _centroidTree->constructWithIterations(2, _centroids);
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
std::vector<Plane>& MeshReconstructor::getPlanes()
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
void MeshReconstructor::setPlanes(std::vector<Plane> &planes)
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
 * @brief MeshReconstructor::computePlanes
 */
void MeshReconstructor::computeCentroidsAndTangentPlanes(Mesh m)
{
    std::cout << "computing centroids" << std::endl;

    // simply calling Octree::findKNearestNeighbours for each point in _mesh
    std::vector<Vertex *>::iterator pointIterator;
    std::vector<std::pair<Vertex *, float>> neighbours;
    int id = 0;

    for (pointIterator = m.getVertices().begin(); pointIterator < m.getVertices().end(); pointIterator++)
    {
        neighbours.clear();
        glm::vec3 centroid(0.f);
        glm::mat3x3 covarianceMatrix(0.f);
        Plane p;
        Matrix3f eigen;
        Matrix3f ev;

        // retrieving the k nearest neighbours of ref
        neighbours = _pointTree->findKNeartestNeighbours((*pointIterator), _k);

        // computing the centroid
        for (unsigned int i = 0; i < neighbours.size(); i++)
        {
            centroid = centroid + neighbours[i].first->getPosition();
        }

        // storing the centroid in _centroids
        centroid = (1.f / neighbours.size()) * centroid;
        _centroids.push_back(new Vertex(id++, centroid.x, centroid.y, centroid.z));
        //std::cout << "centroid : " << centroid.x << ", " << centroid.y << ", " << centroid.z << std::endl;

        // computing the covariance matrix
        for (unsigned int i = 0; i < neighbours.size(); i++)
        {
            // outer product to form the covariance matrix !
            glm::mat3x3 cv;
            glm::vec3 v(neighbours[i].first->getPosition() - centroid);
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
        temp = glm::vec3(ev(0,0), ev(0,1), ev(0,2));
        p.setEigenvector1(temp);
        temp = glm::vec3(ev(1,0), ev(1,1), ev(1,2));
        p.setEigenvector2(temp);
        temp = glm::vec3(ev(2,0), ev(2,1), ev(2,2));
        p.setEigenvector3(temp);
        _planes.push_back(p);
    }
}
