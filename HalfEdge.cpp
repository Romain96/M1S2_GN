/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include "Vertex.h"
#include "Face.h"
#include "HalfEdge.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief HalfEdge::HalfEdge constructs a new Half Edge from an origin Vertex and a destination one
 * @param origin pointer to the Vertex which is the origin of the Half Edge
 * @param destination pointer to the Vertex which is the destination of the Half Edge
 */
HalfEdge::HalfEdge(Vertex *origin, Vertex *destination) :
    _origin(origin),
    _destination(destination),
    _face(nullptr),
    _previousHE(nullptr),
    _nextHE(nullptr),
    _oppositeHE(nullptr)
{
    // nothing
}

/**
 * @brief HalfEdge::HalfEdge constructs a new Half Edge with every informations
 * @param origin pointer to the origin Vertex
 * @param destination pointer to the Vertex which is the destination of the Half Edge
 * @param face pointer to the Face in which the Half Edge belongs
 * @param previous pointer to the previous Half Edge
 * @param next pointer to the next Half Edge
 * @param opposite pointer the the opposite Half Edge
 */
HalfEdge::HalfEdge(Vertex *origin, Vertex *destination, Face *face, HalfEdge *previous, HalfEdge *next, HalfEdge *opposite) :
    _origin(origin),
    _destination(destination),
    _face(face),
    _previousHE(previous),
    _nextHE(next),
    _oppositeHE(opposite)
{
    // nothing
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief HalfEdge::getOrigin
 * @return a pointer to the origin Vertex
 */
Vertex *HalfEdge::getOrigin()
{
    return _origin;
}

/**
 * @brief HalfEdge::getDestination
 * @return a pointer to the destination Vertex
 */
Vertex *HalfEdge::getDestination()
{
    return _destination;
}

/**
 * @brief HalfEdge::getFace
 * @return a pointer to the Face in which the Half Edge belongs
 */
Face *HalfEdge::getFace()
{
    return _face;
}

/**
 * @brief HalfEdge::getPrevious
 * @returna pointer to the previous Half Edge
 */
HalfEdge *HalfEdge::getPrevious()
{
    return _previousHE;
}

/**
 * @brief HalfEdge::getNext
 * @return a pointer to the next Half Edge
 */
HalfEdge *HalfEdge::getNext()
{
    return _nextHE;
}

/**
 * @brief HalfEdge::getOpposite
 * @returna pointer to the opposite Half Edge
 */
HalfEdge *HalfEdge::getOpposite()
{
    return _oppositeHE;
}

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief HalfEdge::setOrigin
 * @param origin new pointer to the origin Vertex
 */
void HalfEdge::setOrigin(Vertex *origin)
{
    _origin = origin;
}

/**
 * @brief HalfEdge::setDestination
 * @param destination new pointer to the destination vertex
 */
void HalfEdge::setDestination(Vertex *destination)
{
    _destination = destination;
}

/**
 * @brief HalfEdge::setFace
 * @param face new pointer to the Face in which the Half Edge belongs
 */
void HalfEdge::setFace(Face *face)
{
    _face = face;
}

/**
 * @brief HalfEdge::setPrevious
 * @param previous new pointer to the previous Half Edge
 */
void HalfEdge::setPrevious(HalfEdge *previous)
{
    _previousHE = previous;
}

/**
 * @brief HalfEdge::setNext
 * @param next new pointer to the newt Half Edge
 */
void HalfEdge::setNext(HalfEdge *next)
{
    _nextHE = next;
}

/**
 * @brief HalfEdge::setOpposite
 * @param opposite new pointer to the opposite Half Edge
 */
void HalfEdge::setOpposite(HalfEdge *opposite)
{
    _oppositeHE = opposite;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------
