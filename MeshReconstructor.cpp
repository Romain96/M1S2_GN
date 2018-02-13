/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <set>

#include "MeshReconstructor.h"

#include "Octree.h"

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
MeshReconstructor::MeshReconstructor(Mesh m)
{
    _mesh = m;
    // creating the tree
    _tree = new Octree();
    _tree->findSpaceBorders(_mesh.getVertices());
    _tree->constructWithIterations(2, _mesh.getVertices());
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

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------

/**
 * @brief MeshReconstructor::computePlanes
 */
void MeshReconstructor::computePlanes()
{
    std::cout << "computing centroids" << std::endl;

    // simply calling Octree::findKNearestNeighbours for each point in _mesh
    std::vector<Vertex *>::iterator pointIterator;
    std::vector<std::pair<Vertex *, float>> neighbours;

    for (pointIterator = _mesh.getVertices().begin(); pointIterator < _mesh.getVertices().end(); pointIterator++)
    {
        neighbours.clear();
        glm::vec3 centroid(0.f);
        glm::mat3x3 covarianceMatrix(0.f);

        // retrieving the k nearest neighbours of ref
        neighbours = _tree->findKNeartestNeighbours((*pointIterator), _k);

        // computing the centroid
        for (unsigned int i = 0; i < neighbours.size(); i++)
        {
            centroid = centroid + neighbours[i].first->getPosition();
        }

        // storing the centroid in _centroids
        centroid = (1.f / neighbours.size()) * centroid;
        _centroids.push_back(centroid);
        std::cout << "centroid : " << centroid.x << ", " << centroid.y << ", " << centroid.z << std::endl;

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

    }
}
