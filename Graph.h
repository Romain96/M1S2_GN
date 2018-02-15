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

#include "Node.h"
#include "Edge.h"
#include "Octree.h"
#include "Vertex.h"

/**
 * @brief The EdgeComparator class comparator class for Edges (map)
 */
class EdgeComparator
{
    bool operator() (std::pair<int, Edge *>e1, std::pair<int, Edge *> e2)
    {
        return e1.second->getWeight() < e2.second->getWeight();
    }
};

/**
 * @brief The Graph class
 */
class Graph
{
protected:
    // storing all nodes (same number as centroids or tangent planes)
    std::vector<Node *> _nodes;
    // storing all edges (number unknown)
    std::map<int, Edge *, EdgeComparator> _edges;

public:
    // constructor(s)
    Graph();

    // getter(s)
    Node *getNodeAtIndex(unsigned int i);
    Edge *getEdgeAtIndex(unsigned int i);
    std::vector<Node *>& getNodes();
    std::map<int, Edge *, EdgeComparator>& getEdges();

    // setter(s)
    void setNodeAtIndex(Node *n, unsigned int i);
    void setEdgeAtIndex(Edge *e, unsigned int i);
    void setNodes(std::vector<Node *>& nodes);
    void setEdges(std::map<int, Edge *, EdgeComparator>& edges);

    // method(s)
    void buildGraph(Octree *t, std::vector<Vertex *>& centroids, std::vector<Plane *>& planes);
    void buildMinimumSpanningTree();
};

#endif // GRAPH_H
