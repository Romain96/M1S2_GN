/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include "Edge.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief Edge::Edge default Edge constructor
 */
Edge::Edge() :
    _left(nullptr),
    _right(nullptr),
    _weight(0)
{
    // nothing
}

/**
 * @brief Edge::Edge parametrized Edge constructor
 * @param left
 * @param right
 * @param weight
 */
Edge::Edge(Node *left, Node *right, float weight) :
    _left(left),
    _right(right),
    _weight(weight)
{
    // nothing
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Edge::getLeftNode
 * @return a pointer to the left Node
 */
Node *Edge::getLeftNode()
{
    return _left;
}

/**
 * @brief Edge::getRightNode
 * @return a pointer to the right Node
 */
Node *Edge::getRightNode()
{
    return _right;
}

/**
 * @brief Edge::getWeight
 * @return the weight of the Edge
 */
float Edge::getWeight()
{
    return _weight;
}

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Edge::setLeftNode
 * @param l new left Node
 */
void Edge::setLeftNode(Node *l)
{
    _left = l;
}

/**
 * @brief Edge::setRightNode
 * @param r new right Node
 */
void Edge::setRightNode(Node *r)
{
    _right = r;
}

/**
 * @brief Edge::setWeight
 * @param w new Edge weight
 */
void Edge::setWeight(float w)
{
    _weight = w;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------
