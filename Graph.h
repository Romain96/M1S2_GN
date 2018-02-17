#ifndef GRAPH_H
#define GRAPH_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// STL
#include <map>
#include <vector>
#include <iostream>

#include "Node.h"
#include "Edge.h"
#include "Octree.h"
#include "Vertex.h"

/**
 * @brief The Graph class
 */
class Graph
{
protected:
    // storing all nodes (same number as centroids or tangent planes)
    std::vector<Node *> _nodes;
    // storing all edges (number unknown)
    std::vector<Edge *> _edges;

public:
    // constructor(s)
    Graph();

    // getter(s)
    Node *getNodeAtIndex(unsigned int i);
    Edge *getEdgeAtIndex(unsigned int i);
    std::vector<Node *>& getNodes();
    std::vector<Edge *>& getEdges();

    // setter(s)
    void setNodeAtIndex(Node *n, unsigned int i);
    void setEdgeAtIndex(Edge *e, unsigned int i);
    void setNodes(std::vector<Node *>& nodes);
    void setEdges(std::vector<Edge *>& edges);

    // method(s)
    void buildGraph(Octree *t, std::vector<Vertex *>& centroids, std::vector<Plane *>& planes);
    void buildMinimumSpanningTree();

private:
    // internal methods
    void __buildNodes(std::vector<Vertex *>& centroids, std::vector<Plane *>& planes);
};

#endif // GRAPH_H
