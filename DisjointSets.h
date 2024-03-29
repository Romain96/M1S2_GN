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
    // list of all parents and their respective ranks
    std::vector<std::pair<Node *, int>> _element;
    // number of nodes in graph
    int n;    

    // constructor
    DisjointSets(int n, std::vector<Node *>& nodes);

    // find parent of Node u
    Node *find_set(Node *u);

    // merge
    int union_set(Node *x, Node *y);
};

#endif // DISJOINTSETS_H
