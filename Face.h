#ifndef FACE_H
#define FACE_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <vector>

#include "Vertex.h"
#include "HalfEdge.h"

// forward declaration
class Vertex;
class HalfEdge;

/**
 * @brief The Face class
 */
class Face
{
protected:
    int _id;
    HalfEdge *_incidentHE;
    std::vector<Vertex *> _vertices;

public:
    // constructor
    Face(int id);
    Face(int id, HalfEdge *incidentHE);

    // getters
    int& getId();
    HalfEdge *getIncidentHE();
    int getVerticesNumber();
    Vertex *getVertex(int i);

    // setters
    void setId(int& id);
    void setIncidentHE(HalfEdge *incidentHE);

    // methods
    void addVertex(Vertex *vertex);
};

#endif // FACE_H
