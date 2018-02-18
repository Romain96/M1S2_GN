#ifndef DISJOINTSETS_H
#define DISJOINTSETS_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include "Node.h"

struct DisjointSets
{
    // list of all parents
    Node **parent;
    // list of all ranks
    int *rnk;
    // number of nodes in graph
    int n;

    // constructor
    DisjointSets(int n, std::vector<Node *>& nodes);

    // find parent of Node u
    Node *DS_find(Node *u);

    // merge
    void DS_merge(Node *x, Node *y);
};

#endif // DISJOINTSETS_H
