#ifndef EDGE_H
#define EDGE_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include "Node.h"

/**
 * @brief The Edge class
 */
class Edge
{
protected:
    // an Edge is between two nodes and has a weight
    Node *_left;
    Node *_right;
    float _weight;

public:
    // constructor(s)
    Edge();
    Edge(Node *left, Node *right, float weight);

    // getter(s)
    Node *getLeftNode();
    Node *getRightNode();
    float getWeight();

    // setter(s)
    void setLeftNode(Node *l);
    void setRightNode(Node *r);
    void setWeight(float w);

    // method(s)
};

#endif // EDGE_H
