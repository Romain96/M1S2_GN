/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <vector>

#include "Octree.h"
#include "Plane.h"
#include "Graph.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

/**
 * @brief Node::Node Node parametrized constructor
 * @param plane
 */
Node::Node(Plane *plane) :
    _plane(plane)
{
    // nothing
}

/**
 * @brief Edge::Edge Edge parametrized constructor
 * @param left "left" Node
 * @param right "right" Node
 * @param weight weight assiciated with the Edge
 */
Edge::Edge(Node *left, Node *right, float weight) :
    _left(left),
    _right(right),
    _weight(weight)
{
    // nothing
}

/**
 * @brief Graph::Graph Graph parametrized constructor
 * @param planes list of all planes
 * @param centroids list of all associated centroids
 */
Graph::Graph(std::vector<Plane>& planes, std::vector<glm::vec3>& centroids)
{
    // construction of an Octree on the centroids
    _tree = Octree();
    _tree.findSpaceBorders(centroids);
    _tree.constructWithIterations(2, centroids);

    // building edges
    // TODO
}

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------
