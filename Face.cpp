/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include "HalfEdge.h"
#include "Face.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief Face::Face constructs a new Face without an incident Face
 * @param id unique id of the Face
 */
Face::Face(int id) :
    _id(id), _incidentHE(nullptr)
{
    // nothing
}

/**
 * @brief Face::Face constructs a new Face with an incident Half Edge
 * @param id unique id of the face
 * @param incidentHE pointer to one of the incident Half Edges belonging to the Face
 */
Face::Face(int id, HalfEdge *incidentHE) :
    _id(id), _incidentHE(incidentHE)
{
    // nothing
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Face::getId
 * @return the unique id of the Face
 */
int& Face::getId()
{
    return _id;
}

/**
 * @brief Face::getIncidentHE
 * @return a pointer to one of the incident Half Edges belonging to the Face
 */
HalfEdge *Face::getIncidentHE()
{
    return _incidentHE;
}

/**
 * @brief Face::getVerticesNumber
 * @return the number of vertices in the Face
 */
int Face::getVerticesNumber()
{
    return _vertices.size();
}

/**
 * @brief Face::getVertex
 * @param i index of the vertex to retrieve
 * @return a pointer to the Vertex at position i
 */
Vertex *Face::getVertex(int i)
{
    return _vertices[i];
}

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Face::setId
 * @param id new unique id of the Face
 */
void Face::setId(int &id)
{
    _id = id;
}

/**
 * @brief Face::setIncidentHE
 * @param incidentHE new pointer to one of the incident Half Edges to the Face
 */
void Face::setIncidentHE(HalfEdge *incidentHE)
{
    _incidentHE = incidentHE;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------

/**
 * @brief Face::addVertex
 * @param vertex vertex to add to the Face
 */
void Face::addVertex(Vertex *vertex)
{
    _vertices.push_back(vertex);
}
