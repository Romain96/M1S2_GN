#ifndef NODE_H
#define NODE_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include "Plane.h"
#include "Vertex.h"

/**
 * @brief The Node class
 */
class Node
{
protected:
    // used to store a centroid
    Vertex *_centroid;
    // an it's associated plane (+ normal)
    Plane *_plane;

public:
    // constructor(s)
    Node();
    Node(Vertex *centroid, Plane *plane);

    // getter(s)
    Vertex *getCentroid();
    Plane *getPlane();

    // setter(s)
    void setCentroid(Vertex *v);
    void setPlane(Plane *p);

    // method(s)
};

#endif // NODE_H
