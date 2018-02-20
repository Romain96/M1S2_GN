#ifndef TREE_H
#define TREE_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// STL
#include "vector"

/**
 * @brief The TreeNode class
 */
class Tree
{
protected:
    // a node has a pointer to his father
    Tree *_father;
    // and has a certain number of children (edges are implict)
    std::vector<Tree *>_children;

public:
    // constructor
    Tree();

    // getter(s)
    Tree *getFather();
    Tree *getChildAtIndex(unsigned int i);
    std::vector<Tree *>& getChildren();

    // setter(s)
    void setFather(Tree *t);
    void setChildAtIndex(Tree *t, unsigned int i);
    void setChildren(std::vector<Tree *>& children);

    // method(s)
    void addChild(Tree *t);
};

#endif // TREE_H
