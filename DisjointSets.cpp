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
    parent = new Node *[n - 1];

    // initially all vertices are in
    // different sets and have a rank of 0
    for (int i = 0; i < n; i++)
    {
        parent[i] = nodes[i];
    }
}

/**
 * @brief DisjointSets::find_set
 * @param u node
 * @return the parent node of u
 */
Node *DisjointSets::find_set(Node *u)
{
    if (u == parent[u->getCentroid()->getId()])
        return parent[u->getCentroid()->getId()];
    else
        return find_set(parent[u->getCentroid()->getId()]);

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
    /*
    Node *xroot = find_set(x);
    Node *yroot = find_set(y);

    if (rnk[xroot->getCentroid()->getId()] < rnk[yroot->getCentroid()->getId()])
        parent[xroot->getCentroid()->getId()] = yroot;
    else if (rnk[xroot->getCentroid()->getId()] > rnk[yroot->getCentroid()->getId()])
        parent[yroot->getCentroid()->getId()] = xroot;
    else    // rank of x == rank of y
    {
        parent[yroot->getCentroid()->getId()] = xroot;
        rnk[xroot->getCentroid()->getId()]++;
    }
    */
    parent[x->getCentroid()->getId()] = parent[y->getCentroid()->getId()];
    return parent[y->getCentroid()->getId()]->getCentroid()->getId();
}

