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
    _edges(),
    _root(nullptr)
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

/**
 * @brief Graph::getRoot
 * @return a pointer to the root Node
 */
Node *Graph::getRoot()
{
    return _root;
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

/**
 * @brief Graph::setRoot
 * @param root pointer to new root Node
 */
void Graph::setRoot(Node *root)
{
    _root = root;
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
 * @brief Graph::addEdgeWithRedundancyTest
 * @param e new Edge to add
 */
void Graph::addEdgeWithRedundancyTest(Edge *e)
{
    // iterating on each edges and addint the edge if
    // and only if it isn't already in the list of edges
    std::vector<Edge *>::iterator it;

    for (it = _edges.begin(); it != _edges.end(); it++)
    {
        // we have to test both sides since the graph is undirected
        if ((e->getLeftNode()->getCentroid()->getId() == (*it)->getLeftNode()->getCentroid()->getId() &&
                e->getRightNode()->getCentroid()->getId() == (*it)->getRightNode()->getCentroid()->getId())
                ||
                (e->getLeftNode()->getCentroid()->getId() == (*it)->getRightNode()->getCentroid()->getId() &&
                 e->getRightNode()->getCentroid()->getId() == (*it)->getLeftNode()->getCentroid()->getId()))
        {
                return;
        }
    }

    // if we are here, the edge is not already in the list so we add it
    _edges.push_back(e);
}

/**
 * @brief Graph::clearNodes
 */
void Graph::clearNodes()
{
    _nodes.clear();
}

/**
 * @brief Graph::clearEdges
 */
void Graph::clearEdges()
{
    _edges.clear();
}

/**
 * @brief Graph::buildEuclidianGraph builds an Euclidian graph (based on Euclidian distance)
 * @param centroids list of all computed centroids
 * @param planes list of all computed tangent planes
 */
void Graph::buildEuclidianGraph(std::vector<Vertex *>& centroids, std::vector<Plane *>& planes)
{
    // building all nodes before the edges
    __buildNodes(centroids, planes);

    std::vector<Node *>::iterator nodeIterator1;
    std::vector<Node *>::iterator nodeIterator2;
    Edge *e;
    float weight;

    // foe each centroid
    for (nodeIterator1 = _nodes.begin(); nodeIterator1 != _nodes.end(); nodeIterator1++)
    {
        nodeIterator2 = nodeIterator1;
        nodeIterator2++;
        // for each centroid > referencial centroid (otherwise already treated and added an edge)
        while (nodeIterator2 != _nodes.end())
        {  
            // compute weight between two nodes and create an edge (Euclidian distance between centroids)
            weight = Vertex::distance3((*nodeIterator1)->getCentroid()->getPosition(),
                                       (*nodeIterator2)->getCentroid()->getPosition());

            e = new Edge((*nodeIterator1), (*nodeIterator2), weight);
            _edges.push_back(e);

            nodeIterator2++;
        }
    }

    std::cout << _edges.size() << " edges created with full graph" << std::endl;
}

/**
 * @brief The sortCompare struct for std:sort used in buildMinimumSpanningTree
 */
struct sortCompare
{
    bool operator() (Edge *e1, Edge *e2)
    {
        return e1->getWeight() < e2->getWeight();
    }
};

/**
 * @brief Graph::buildMinimumSpanningTree builds a MST using Kruskal's algorithm
 * @return a new Graph which is a MST
 */
Graph *Graph::buildMinimumSpanningTree()
{
    // building a new Graph and adding each nodes of the current one
    Graph *test = new Graph();
    for (unsigned int i = 0; i < _nodes.size(); i++)
        test->addNode(_nodes[i]);

    // sorting edges in increasing weight order
    std::sort(_edges.begin(), _edges.end(), sortCompare());

    // creating disjoint sets
    DisjointSets ds(_nodes.size(), _nodes);
    int rootId = 0;

    std::vector<Edge *>::iterator edgeIterator;
    for (edgeIterator = _edges.begin(); edgeIterator != _edges.end(); edgeIterator++)
    {
        Node *u = (*edgeIterator)->getLeftNode();
        Node *v = (*edgeIterator)->getRightNode();

        Node *uRep = ds.find_set(u);
        Node *vRep = ds.find_set(v);

        if (uRep != vRep)
        {
            // adding this edge to the MST
            test->addEdge((*edgeIterator));

            // merge two sets
            rootId = ds.union_set(uRep, vRep);
        }
    }
    test->setRoot(_nodes[rootId]);

    std::cout << "MST root is " << rootId << std::endl;
    std::cout << "MST has " << test->getNodes().size() << " nodes" << std::endl;
    std::cout << "MST has " << test->getEdges().size() << " edges" << std::endl;
    return test;
}

/**
 * @brief Graph::enhanceToRiemannianGraph
 * @param k number of neighbours used in k neighbours search
 * @param t Octree on centroids used in k neighbours search
 */
void Graph::enhanceToRiemannianGraph(int k, Octree *t)
{
    // the objective is to add to the MST
    // edges between centroids/tangent planes who are
    // k neighbours

    std::vector<std::pair<Vertex *, float>> neighbours;
    std::vector<Node *>::iterator nodeIt;
    std::vector<Edge *>::iterator edgeIt;

    // first and foremost we need to change the weight in the MST
    // since the information is not the same
    // right now the Euclidian distance is represented
    // but we want now a degree of closeness in term of tangent planes
    for (edgeIt = _edges.begin(); edgeIt != _edges.end(); edgeIt++)
    {
        float weight = 1.f - fabs(glm::dot( (*edgeIt)->getLeftNode()->getPlane()->getEigenvector3(),
                                      (*edgeIt)->getRightNode()->getPlane()->getEigenvector3()));
        (*edgeIt)->setWeight(weight);
    }

    // for each node we find the k nearest neighbours
    for (nodeIt = _nodes.begin(); nodeIt != _nodes.end(); nodeIt++)
    {
        // retrieving the k nearest neighbours
        neighbours.clear();
        neighbours = t->findKNeartestNeighbours((*nodeIt)->getCentroid(), k);

        // for each neighbours we create an edge
        // the weight is already computed
        for (unsigned int i = 0; i < neighbours.size(); i++)
        {
            float weight = 1.f - fabs(glm::dot((*nodeIt)->getPlane()->getEigenvector3(),
                                               _nodes[neighbours[i].first->getId()]->getPlane()->getEigenvector3()));
            Edge *e = new Edge((*nodeIt), _nodes[neighbours[i].first->getId()], weight);

            // adding the edge with a redundancy test
            addEdgeWithRedundancyTest(e);
        }
    }
    std::cout << "Riemannian graph has " << _edges.size() << " edges" << std::endl;
}

/**
 * @brief Graph::traverseDepthFirstAndReorientate
 */
void Graph::traverseDepthFirstAndReorientate()
{
    // the objective is to traverse the MST in depth first order
    // and propagate the "right" orientation

    // initializing orientation for largest Z component plane
    // and rooting the tree at this node
    __findLargestZNormalAndRootTreeAtNode();

    // initial recursive call
    __depthFirstTraversingAndReorientation(_root, _root);
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
    _nodes.clear();

    while (centroidIterator != centroids.end() && planeIterator != planes.end())
    {
        n = new Node((*centroidIterator),(*planeIterator));
        _nodes.push_back(n);

        centroidIterator++;
        planeIterator++;
    }

    std::cout << _nodes.size() << " nodes created" << std::endl;
}

/**
 * @brief Graph::__findChildrenOfNode
 * @param n Node
 * @return a list of all Nodes that are children of n
 */
std::vector<Node *>& Graph::__findChildrenOfNode(Node *n)
{
    static std::vector<Node *> children;
    children.clear();

    std::vector<Edge *>::iterator edgeIt;

    // finding all edges from in which n is the origin
    for (edgeIt = _edges.begin(); edgeIt != _edges.end(); edgeIt++)
    {
        if ((*edgeIt)->getLeftNode() == n)
        {
            children.push_back((*edgeIt)->getRightNode());
        }
    }

    return children;
}

/**
 * @brief Graph::__depthFirstTraversingAndReorientation
 * @param parent pointer to parent Node (caller of recursive function)
 * @param current pointer to currently treated Node
 */
void Graph::__depthFirstTraversingAndReorientation(Node *parent, Node *current)
{
    std::cerr << "__dftar " << parent->getCentroid()->getId() << ", " << current->getCentroid()->getId() << std::endl;
    // if Ni . Nj < 0 Nj is replaced by -Nj where
    // * Ni is the Normal of the tangent plane associated with the Node i
    // * Nj is the Normal of the tangent plane associated with the Node j
    if (glm::dot(parent->getPlane()->getEigenvector3(), current->getPlane()->getEigenvector3()) < 0.f)
    {
        glm::vec3 normal = -current->getPlane()->getEigenvector3();
        current->getPlane()->setEigenvector3(normal);
    }

    // retrieving each children of the current Node
    std::vector<Node *> children = __findChildrenOfNode(current);
    std::vector<Node *>::iterator childrenIt;

    // recursively calling the function in depth first order for each children
    for (childrenIt = children.begin(); childrenIt != children.end(); childrenIt++)
    {
        __depthFirstTraversingAndReorientation(current, (*childrenIt));
    }
}

/**
 * @brief Graph::__findLargestZNormalAndRootTreeAtNode
 */
void Graph::__findLargestZNormalAndRootTreeAtNode()
{
    // first we find the plane who has the largest z component normal (eigenvector 3)
    // let's call it maxZ
    std::vector<Node *>::iterator nodeIt;
    Node *maxZ = nullptr;
    float maxZValue = 0.f;  // as absolute value !

    for (nodeIt = _nodes.begin(); nodeIt != _nodes.end(); nodeIt++)
    {
        if (fabs((*nodeIt)->getPlane()->getEigenvector3().z) > maxZValue)
        {
            maxZ = (*nodeIt);
            maxZValue = fabs((*nodeIt)->getPlane()->getEigenvector3().z);
        }
    }

    std::cout << "largest Z is node " << maxZ->getCentroid()->getId() << std::endl;

    // if z is pointing towards -Z then we make it pointing towards +Z
    if (maxZ->getPlane()->getEigenvector3().z < 0.f)
    {
        glm::vec3 normal = maxZ->getPlane()->getEigenvector3();
        normal.z = -normal.z;
        maxZ->getPlane()->setEigenvector3(normal);
    }

    // if maxZ plane has a parent in the MST
    // we delete the edge parent-maxZ
    // we root the tree at maxZ by creating the edge maxZ-root
    // where root is the current root of the MST
    std::vector<Edge *>::iterator edgeIt;

    for (edgeIt = _edges.begin(); edgeIt != _edges.end(); edgeIt++)
    {
        if( (*edgeIt)->getRightNode() == maxZ)
        {
            std::cerr << "swaping " << (*edgeIt)->getLeftNode()->getCentroid()->getId()
                      << " and " << _root->getCentroid()->getId() << std::endl;
            // rooting at maxZ
            float weight = 1.f - fabs(glm::dot(maxZ->getPlane()->getEigenvector3(),
                                      (*edgeIt)->getLeftNode()->getPlane()->getEigenvector3()));
            Edge *e = new Edge(maxZ, (*edgeIt)->getLeftNode(), weight);
            addEdge(e);

            // removing parent-maxZ edge
            _edges.erase(edgeIt);

            // updating the root
            setRoot(maxZ);
            break;
        }
    }

    std::cout << "root is now " << _root->getCentroid()->getId() << std::endl;
}
