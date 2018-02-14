#ifndef GRAPH_H
#define GRAPH_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// STL
#include <vector>

// GLM
#include "glm/glm.hpp"
#include "glm/vec3.hpp"

#include "Plane.h"
#include "Octree.h"

/**
 * @brief The Node class
 */
class Node
{
protected:
    Plane *_plane;

public:
    // constructor(s)
    Node(Plane *plane);

    // getter(s)
    Plane *getPlane();

    // setter(s)
    void setPlane(Plane *plane);

    // method(s)
};

/**
 * @brief The Edge class
 */
class Edge
{
protected:
    Node *_left;
    Node *_right;
    float _weight;

public:
    // constructor(s)
    Edge(Node *left, Node *right, float weight);

    // getter(s)
    Node *getLeftNode();
    Node *getRightNode();
    float getWeight();

    // setter(s)
    void setLeftNode(Node *left);
    void setRightNode(Node *right);
    void setWeight(float weight);

    // method(s)

};

/**
 * @brief The Graph class
 */
class Graph
{
protected:
    // just storing a vector of Edges
    std::vector<Edge *> _edges;

    // and keeping an Octree on the centroids
    Octree _tree;

public:
    // constructor(s)
    Graph(std::vector<Plane>& planes, std::vector<glm::vec3>& centroids);

    // getter(s)
    Edge *getEdge(int i);

    // setter(s)
    void setEdge(Edge *edge);

    // method(s)
    // TODO Kruskal to compute Minimal Spanning Tree

private:
    void __buildEdges();
};

#endif // GRAPH_H
