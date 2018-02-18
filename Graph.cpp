/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// STL
#include <vector>
#include <map>
#include <algorithm>
#include <bits/stdc++.h>

#include "Graph.h"
#include "DisjointSets.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief Graph::Graph default Graph constructor
 */
Graph::Graph() :
    _nodes(),
    _edges()
{
    // nothing
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Graph::getNodeAtIndex
 * @param i index
 * @return a pointer to the node at index i or nullptr if this node doesn't exist
 */
Node *Graph::getNodeAtIndex(unsigned int i)
{
    if (i >= _nodes.size())
        return nullptr;
    else
        return _nodes[i];
}

/**
 * @brief Graph::getEdgeAtIndex
 * @param i index
 * @return a pointer to the Edge at index i or nullptr if this Edge doesn't exist
 */
Edge *Graph::getEdgeAtIndex(unsigned int i)
{
    if (i >= _edges.size())
            return nullptr;
    else
            return _edges[i];
}

/**
 * @brief Graph::getNodes
 * @return a vector containing all Nodes
 */
std::vector<Node *>& Graph::getNodes()
{
    return _nodes;
}

/**
 * @brief Graph::getEdges
 * @return a map containing all Edges
 */
std::vector<Edge *>& Graph::getEdges()
{
    return _edges;
}


//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Graph::setNodeAtIndex
 * @param n new Node
 * @param i index
 */
void Graph::setNodeAtIndex(Node *n, unsigned int i)
{
    // only if index is not out of range
    if (i < _nodes.size())
        _nodes[i] = n;
}

/**
 * @brief Graph::setEdgeAtIndex
 * @param e new Edge
 * @param i index
 */
void Graph::setEdgeAtIndex(Edge *e, unsigned int i)
{
    // only if index is not out of range
    if (i < _edges.size())
        _edges[i] = e;
}

/**
 * @brief Graph::setNodes
 * @param nodes a vector containing all Nodes
 */
void Graph::setNodes(std::vector<Node *> &nodes)
{
    _nodes = nodes;
}

/**
 * @brief Graph::setEdges
 * @param edges a map containing all Edges
 */
void Graph::setEdges(std::vector<Edge *> &edges)
{
    _edges = edges;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------

/**
 * @brief Graph::addNode
 * @param n new Node
 */
void Graph::addNode(Node *n)
{
    _nodes.push_back(n);
}

/**
 * @brief Graph::addEdge
 * @param e new Edge
 */
void Graph::addEdge(Edge *e)
{
    _edges.push_back(e);
}

/**
 * @brief Graph::buildGraph
 * @param k number of neighbours (kNearestNeighbours search)
 * @param t Octree built upon centroids
 * @param centroids list of all centroids
 * @param planes list of all associated tangent planes
 */
void Graph::buildGraph(int k, Octree *t, std::vector<Vertex *> &centroids, std::vector<Plane *> &planes)
{
    // building all nodes before the edges
    __buildNodes(centroids, planes);

    std::vector<Node *>::iterator nodeIterator = _nodes.begin();
    std::vector<std::pair<Vertex *, float>> neighbours;
    Edge *e;
    float weight;

    // For each node we find the k nearest neighbours
    while (nodeIterator != _nodes.end())
    {
        // retrieving the k nearest neighbours of the centroid
        neighbours.clear();
        neighbours = t->findKNeartestNeighbours((*nodeIterator)->getCentroid(), k);

        // for each neighbours which doesn't have been treated
        for (unsigned int i = 0; i < neighbours.size(); i++)
        {
            // if centroid ID is inferior to the current treated centroid then
            // the edge between it and the current one has already been created
            if (neighbours[i].first->getId() > (*nodeIterator)->getCentroid()->getId())
            {
                // the weight is 1 - |ni . nj| where
                // ni is the normal of the tangent plane of the currently treated node
                // nj is the normal of the tangent plane of the currently treated neighbour node
                weight = 1.f - fabs(glm::dot((*nodeIterator)->getPlane()->getEigenvector3(),
                                             _nodes[neighbours[i].first->getId()]->getPlane()->getEigenvector3()));
                e = new Edge((*nodeIterator), _nodes[neighbours[i].first->getId()], weight);

                _edges.push_back(e);
            }
        }

        nodeIterator++;
    }

    std::cout << _edges.size() << " edges created" << std::endl;

    /*
    for (unsigned int i = 0; i < _edges.size(); i++)
    {
        std::cout << "edge " << _edges[i]->getLeftNode()->getCentroid()->getId() << " to "
                  << _edges[i]->getRightNode()->getCentroid()->getId() << " with weight"
                  << _edges[i]->getWeight() << std::endl;
    }
    */
}

/**
 * @brief Graph::buildMinimumSpanningTree builds a MST using Kruskal's algorithm
 */
void Graph::buildMinimumSpanningTree()
{
    // objective is to minimize the weight of the graph
    float mstWeight = 0.f;

    // sorting edges in increasing weight order
    std::sort(_edges.begin(), _edges.end());

    // creating disjoint sets
    DisjointSets ds(_nodes.size(), _nodes);

    std::vector<Edge *>::iterator edgeIterator;
    for (edgeIterator = _edges.begin(); edgeIterator != _edges.end(); edgeIterator++)
    {
        Node *u = (*edgeIterator)->getLeftNode();
        Node *v = (*edgeIterator)->getRightNode();

        Node *set_u = ds.DS_find(u);
        Node *set_v = ds.DS_find(v);

        if (set_u != set_v)
        {
            // updating MST weight
            mstWeight += (*edgeIterator)->getWeight();

            // merge two sets
            ds.DS_merge(set_u, set_v);
        }
    }
}

//-----------------------------------------------------------------------------
// Internal Method(s)
//-----------------------------------------------------------------------------

/**
 * @brief Graph::__buildNodes
 * @param centroids list of all centroids
 * @param planes list of all tangent planes
 */
void Graph::__buildNodes(std::vector<Vertex *> &centroids, std::vector<Plane *> &planes)
{
    // just creating one node per centroid/tangent plane
    // centroids and tangent planes are stored in the same order points were
    // so the resulting node list will be ordered as so

    Node *n;
    std::vector<Vertex *>::iterator centroidIterator = centroids.begin();
    std::vector<Plane *>::iterator planeIterator = planes.begin();
    std::cout << "building nodes" << std::endl;

    while (centroidIterator != centroids.end() && planeIterator != planes.end())
    {
        n = new Node((*centroidIterator),(*planeIterator));
        _nodes.push_back(n);

        centroidIterator++;
        planeIterator++;
    }

    std::cout << _nodes.size() << " nodes created" << std::endl;
}
