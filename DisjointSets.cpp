/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

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
    parent = new Node *[n + 1];
    rnk = new int[n + 1];

    // initially all vertices are in
    // different sets and have a rank of 0
    for (int i = 0; i <= n; i++)
    {
        rnk[i] = 0;
        parent[i] = nodes[i];
    }
}

/**
 * @brief DisjointSets::DS_find
 * @param u node
 * @return the parent node of u
 */
Node *DisjointSets::DS_find(Node *u)
{
    if (u != parent[u->getCentroid()->getId()])
        parent[u->getCentroid()->getId()] = DS_find(parent[u->getCentroid()->getId()]);
    return parent[u->getCentroid()->getId()];
}

/**
 * @brief DisjointSets::DS_merge
 * @param x node x
 * @param y node y
 */
void DisjointSets::DS_merge(Node *x, Node *y)
{
    Node *px = DS_find(x);
    Node *py = DS_find(y);

    if (rnk[px->getCentroid()->getId()] > rnk[py->getCentroid()->getId()])
        parent[py->getCentroid()->getId()] = px;
    else
        parent[px->getCentroid()->getId()] = py;

    if (rnk[px->getCentroid()->getId()] == rnk[py->getCentroid()->getId()])
        rnk[py->getCentroid()->getId()]++;
}

