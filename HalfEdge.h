#ifndef HALFEDGE_H
#define HALFEDGE_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include "Vertex.h"
#include "Face.h"

// forward declaration
class Vertex;
class Face;

/**
 * @brief The HalfEdge class
 */
class HalfEdge
{
protected:
    // origin of half edge : vertex
    Vertex *_origin;
    // destination of Half Edge : vertex
    Vertex *_destination;
    // face in which the half edge belongs
    Face *_face;
    // previous half edge (from current HE)
    HalfEdge *_previousHE;
    // next half edge (from current HE)
    HalfEdge *_nextHE;
    // opposite half edge (to the current HE)
    HalfEdge *_oppositeHE;

public:
    // constructor(s)
    HalfEdge(Vertex *origin, Vertex *destination);
    HalfEdge(Vertex *origin, Vertex *destination, Face *face, HalfEdge *previous, HalfEdge *next, HalfEdge *opposite);

    // getters
    Vertex *getOrigin();
    Vertex *getDestination();
    Face *getFace();
    HalfEdge *getPrevious();
    HalfEdge *getNext();
    HalfEdge *getOpposite();

    // setters
    void setOrigin(Vertex *origin);
    void setDestination(Vertex *destination);
    void setFace(Face *face);
    void setPrevious(HalfEdge *previous);
    void setNext(HalfEdge *next);
    void setOpposite(HalfEdge *opposite);

    // methods
};

#endif // HALFEDGE_H
