/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// STL
#include <iostream>

#include "Node.h"
#include "DisjointSets.h"

/**
 * @brief DisjointSets::DisjointSets
 * @param n number of nodes in graph
 * @param nodes list of all nodes
 */
DisjointSets::DisjointSets(int n, std::vector<Node *> &nodes)
{
    this->n = n;
    // initially all vertices are in
    // different sets and have a rank of 0
    for (int i = 0; i < n; i++)
    {
        std::pair<Node *, int> p(nodes[i], 0);
        _element.push_back(p);
    }
}

/**
 * @brief DisjointSets::find_set
 * @param u node
 * @return the parent node of u
 */
Node *DisjointSets::find_set(Node *u)
{
    if (u == _element[u->getCentroid()->getId()].first)
        return u;
    else
        return find_set(_element[u->getCentroid()->getId()].first);

    /*
    if (u != parent[u->getCentroid()->getId()])
        parent[u->getCentroid()->getId()] = find_set(parent[u->getCentroid()->getId()]);
    return parent[u->getCentroid()->getId()];*/
}

/**
 * @brief DisjointSets::union_set
 * @param x node x
 * @param y node y
 * @return id of root node
 */
int DisjointSets::union_set(Node *x, Node *y)
{
    Node *xRep = find_set(x);
    Node *yRep = find_set(y);

    if (xRep != yRep)
    {
        // rank of x > rank of y
        if (_element[xRep->getCentroid()->getId()].second > _element[yRep->getCentroid()->getId()].second)
        {
            _element[yRep->getCentroid()->getId()].first = _element[xRep->getCentroid()->getId()].first;
            return xRep->getCentroid()->getId();
        }
        // rank of x < rank of y
        else if (_element[xRep->getCentroid()->getId()].second <_element[yRep->getCentroid()->getId()].second)
        {
            _element[xRep->getCentroid()->getId()].first = _element[yRep->getCentroid()->getId()].first;
            return yRep->getCentroid()->getId();
        }
        // rank x == rank y
        else
        {
            // rooting in x, increaing rank of x
            _element[yRep->getCentroid()->getId()].first = _element[xRep->getCentroid()->getId()].first;
            _element[xRep->getCentroid()->getId()].second += 1;
            return xRep->getCentroid()->getId();
        }
    }
}

