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
        neighbours = _tree->findKNeartestNeighbours((*pointIterator), _k);

        glm::vec3 centroid;
        for (int i = 0; i < neighbours.size(); i++)
            centroid = centroid + neighbours[i].first->getPosition();
        centroid = (1.f / neighbours.size()) * centroid;

        _centroids.push_back(centroid);

        std::cout << "centroid : " << centroid.x << ", " << centroid.y << ", " << centroid.z << std::endl;
    }



}
