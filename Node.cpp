/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include "Vertex.h"
#include "Plane.h"
#include "Node.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief Node::Node default Node constructor
 */
Node::Node() :
    _centroid(nullptr),
    _plane(nullptr)
{
    // nothing
}

/**
 * @brief Node::Node parametrized Node constructor
 * @param centroid centroid of the Node
 * @param plane tangent plane (+ normal) associated with the centroid
 */
Node::Node(Vertex *centroid, Plane *plane) :
    _centroid(centroid),
    _plane(plane)
{
    // nothing
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Node::getCentroid
 * @return a pointer to the centroid of the Node
 */
Vertex *Node::getCentroid()
{
    return _centroid;
}

/**
 * @brief Node::getPlane
 * @return a pointer to the tangent plane of the Node
 */
Plane *Node::getPlane()
{
    return _plane;
}

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Node::setCentroid
 * @param v new centroid
 */
void Node::setCentroid(Vertex *v)
{
    _centroid = v;
}

/**
 * @brief Node::setPlane
 * @param p new tangent plane
 */
void Node::setPlane(Plane *p)
{
    _plane = p;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------
