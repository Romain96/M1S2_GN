/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// glm
#include "glm/glm.hpp"
#include "glm/vec3.hpp"

// stl
#include <iostream>
#include <vector>
#include <utility>
#include <queue>

#include "Vertex.h"
#include "Octree.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

float Octree::_minSize = 0.f;
int Octree::_depth = 10;
bool Octree::_depthTest = true;
bool Octree::_sizeTest = false;

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief Octree::Octree default constructor (useful for root)
 */
Octree::Octree() :
    _borderLowerNW(glm::vec3(0.f, 0.f, 0.f)),
    _borderLowerNE(glm::vec3(0.f, 0.f, 0.f)),
    _borderLowerSW(glm::vec3(0.f, 0.f, 0.f)),
    _borderLowerSE(glm::vec3(0.f, 0.f, 0.f)),
    _borderUpperNW(glm::vec3(0.f, 0.f, 0.f)),
    _borderUpperNE(glm::vec3(0.f, 0.f, 0.f)),
    _borderUpperSW(glm::vec3(0.f, 0.f, 0.f)),
    _borderUpperSE(glm::vec3(0.f, 0.f, 0.f)),
    _isLeaf(true),
    _father(nullptr),
    _lowerNW(nullptr),
    _lowerNE(nullptr),
    _lowerSW(nullptr),
    _lowerSE(nullptr),
    _upperNW(nullptr),
    _upperNE(nullptr),
    _upperSW(nullptr),
    _upperSE(nullptr),
    _points()
{
    // nothing
}

/**
 * @brief Octree::Octree boundaries parametrized constructor (useful for creating children from parent)
 * @param lowerNW lower north west boundary point
 * @param lowerNE lower north east boundary point
 * @param lowerSW lower north west boundary point
 * @param lowerSE lower north east boundary point
 * @param upperNW upper south west boundary point
 * @param upperNE upper south east boundary point
 * @param upperSW upper south west boundary point
 * @param upperSE upper south east boundary point
 */
Octree::Octree(glm::vec3 &lowerNW, glm::vec3 &lowerNE, glm::vec3 &lowerSW, glm::vec3 &lowerSE,
               glm::vec3 &upperNW, glm::vec3 &upperNE, glm::vec3 &upperSW, glm::vec3 &upperSE) :
    _borderLowerNW(lowerNW),
    _borderLowerNE(lowerNE),
    _borderLowerSW(lowerSW),
    _borderLowerSE(lowerSE),
    _borderUpperNW(upperNW),
    _borderUpperNE(upperNE),
    _borderUpperSW(upperSW),
    _borderUpperSE(upperSE),
    _isLeaf(true),
    _father(nullptr),
    _lowerNW(nullptr),
    _lowerNE(nullptr),
    _lowerSW(nullptr),
    _lowerSE(nullptr),
    _upperNW(nullptr),
    _upperNE(nullptr),
    _upperSW(nullptr),
    _upperSE(nullptr),
    _points()
{
    /*std::cout << "Octree dim (diagonal) is "
              << _borderLowerSW.x << "," << _borderLowerSW.y << "," << _borderLowerSW.z
              << " to "
              << _borderUpperNE.x << "," << _borderUpperNE.y << "," << _borderUpperNE.z << std::endl;*/
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Octree::getBorderLowerNW
 * @return the point which lies at the extreme lower north west of the node boundary
 */
glm::vec3& Octree::getBorderLowerNW()
{
    return _borderLowerNW;
}

/**
 * @brief Octree::getBorderLowerNE
 * @return the point which lies at the extreme lower north east of the node boundary
 */
glm::vec3& Octree::getBorderLowerNE()
{
    return _borderLowerNE;
}

/**
 * @brief Octree::getBorderLowerSW
 * @return the point which lies at the extreme lower south west of the node boundary
 */
glm::vec3& Octree::getBorderLowerSW()
{
    return _borderLowerSW;
}

/**
 * @brief Octree::getBorderLowerSE
 * @return the point which lies at the extreme lower south east of the node boundary
 */
glm::vec3& Octree::getBorderLowerSE()
{
    return _borderLowerSE;
}

/**
 * @brief Octree::getBorderUpperNW
 * @return the point which lies at the extreme upper north west of the node boundary
 */
glm::vec3& Octree::getBorderUpperNW()
{
    return _borderUpperNW;
}

/**
 * @brief Octree::getBorderUpperNE
 * @return the point which lies at the extreme upper north east of the node boundary
 */
glm::vec3& Octree::getBorderUpperNE()
{
    return _borderUpperNE;
}

/**
 * @brief Octree::getBorderUpperSW
 * @return the point which lies at the extreme upper south west of the node boundary
 */
glm::vec3& Octree::getBorderUpperSW()
{
    return _borderUpperSW;
}

/**
 * @brief Octree::getBorderUpperSE
 * @return the point which lies at the extreme upper south east of the node boundary
 */
glm::vec3& Octree::getBorderUpperSE()
{
    return _borderUpperSE;
}

/**
 * @brief Octree::getFather
 * @return a pointer to the father node
 */
Octree *Octree::getFather()
{
    return _father;
}

/**
 * @brief Octree::getLowerNW
 * @return the partition which lies at the lower north west position in the current cube
 */
Octree *Octree::getLowerNW()
{
    return _lowerNW;
}

/**
 * @brief Octree::getLowerNE
 * @return the partition which lies at the lower north east position in the current cube
 */
Octree *Octree::getLowerNE()
{
    return _lowerNE;
}

/**
 * @brief Octree::getLowerSW
 * @return the partition which lies at the lower south west position in the current cube
 */
Octree *Octree::getLowerSW()
{
    return _lowerSW;
}

/**
 * @brief Octree::getLowerSE
 * @return the partition which lies at the lower south east position in the current cube
 */
Octree *Octree::getLowerSE()
{
    return _lowerSE;
}

/**
 * @brief Octree::getUpperNW
 * @return the partition which lies at the upper north west position in the current cube
 */
Octree *Octree::getUpperNW()
{
    return _upperNW;
}

/**
 * @brief Octree::getUpperNE
 * @return the partition which lies at the upper north east position in the current cube
 */
Octree *Octree::getUpperNE()
{
    return _upperNE;
}

/**
 * @brief Octree::getUpperSW
 * @return the partition which lies at the upper south west position in the current cube
 */
Octree *Octree::getUpperSW()
{
    return _upperSW;
}

/**
 * @brief Octree::getUpperSE
 * @return the partition which lies at the upper south east position in the current cube
 */
Octree *Octree::getUpperSE()
{
    return _upperSE;
}

/**
 * @brief Octree::getPoints
 * @return the list of points belonging to a node of the tree (only leraves have points)
 */
std::vector<Vertex *>& Octree::getPoints()
{
    return _points;
}

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Octree::setBorderLowerNW
 * @param v new lower north west boundary of the node
 */
void Octree::setBorderLowerNW(glm::vec3& v)
{
    _borderLowerNW = v;
}

/**
 * @brief Octree::setBorderLowerNE
 * @param v new lower north east boundary of the node
 */
void Octree::setBorderLowerNE(glm::vec3& v)
{
    _borderLowerNE = v;
}

/**
 * @brief Octree::setBorderLowerSW
 * @param v new lower south west boundary of the node
 */
void Octree::setBorderLowerSW(glm::vec3& v)
{
    _borderLowerSW = v;
}

/**
 * @brief Octree::setBorderLowerSE
 * @param v new lower south east boundary of the node
 */
void Octree::setBorderLowerSE(glm::vec3& v)
{
    _borderLowerSE = v;
}

/**
 * @brief Octree::setBorderUpperNW
 * @param v new upper north west boundary of the node
 */
void Octree::setBorderUpperNW(glm::vec3& v)
{
    _borderUpperNW = v;
}

/**
 * @brief Octree::setBorderUpperNE
 * @param v new upper north east boundary of the node
 */
void Octree::setBorderUpperNE(glm::vec3& v)
{
    _borderUpperNE = v;
}

/**
 * @brief Octree::setBorderUpperSW
 * @param v new upper south west boundary of the node
 */
void Octree::setBorderUpperSW(glm::vec3& v)
{
    _borderUpperSW = v;
}

/**
 * @brief Octree::setBorderUpperSE
 * @param v new upper south east boundary of the node
 */
void Octree::setBorderUpperSE(glm::vec3& v)
{
    _borderUpperSE = v;
}

/**
 * @brief Octree::setFather
 * @param t new pointer to father node
 */
void Octree::setFather(Octree *t)
{
    _father = t;
}

/**
 * @brief Octree::setLowerNW
 * @param t new lower north west partition
 */
void Octree::setLowerNW(Octree *t)
{
    if (_lowerNW != nullptr)
        delete _lowerNW;
    _lowerNW = t;
}

/**
 * @brief Octree::setLowerNE
 * @param t new lower north east partition
 */
void Octree::setLowerNE(Octree *t)
{
    if (_lowerNE != nullptr)
        delete _lowerNE;
    _lowerNE = t;
}

/**
 * @brief Octree::setLowerSW
 * @param t new lower south west partition
 */
void Octree::setLowerSW(Octree *t)
{
    if (_lowerSW != nullptr)
        delete _lowerSW;
    _lowerSW = t;
}

/**
 * @brief Octree::setLowerSE
 * @param t new lower south east partition
 */
void Octree::setLowerSE(Octree *t)
{
    if (_lowerSE != nullptr)
        delete _lowerSE;
    _lowerSE = t;
}

/**
 * @brief Octree::setUpperNW
 * @param t new upper north west partition
 */
void Octree::setUpperNW(Octree *t)
{
    if (_upperNW != nullptr)
        delete _upperNW;
    _upperNW = t;
}

/**
 * @brief Octree::setUpperNE
 * @param t new upper north east partition
 */
void Octree::setUpperNE(Octree *t)
{
    if (_upperNE != nullptr)
        delete _upperNE;
    _upperNE = t;
}

/**
 * @brief Octree::setUpperSW
 * @param t new upper south west partition
 */
void Octree::setUpperSW(Octree *t)
{
    if (_upperSW != nullptr)
        delete _upperSW;
    _upperSW = t;
}

/**
 * @brief Octree::setUpperSE
 * @param t new upper south east partition
 */
void Octree::setUpperSE(Octree *t)
{
    if (_upperSE != nullptr)
        delete _upperSE;
    _upperSE = t;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------

/**
 * @brief Octree::addPoint adds a point to a node
 * @param v the point to add
 */
void Octree::addPoint(Vertex *v)
{
    _points.push_back(v);
}

/**
 * @brief Octree::isInside tests is a point is inside a node
 * @param v point to test
 * @param t node of the octree
 * @return true if v is inside the subregion defined by t
 */
bool Octree::isInside(Vertex *v, Octree *t)
{
    if (v->getPosition().x >= t->getBorderLowerSW().x && v->getPosition().x <= t->getBorderLowerSE().x &&
                v->getPosition().y >= t->getBorderLowerSE().y && v->getPosition().y <= t->getBorderLowerNE().y &&
                v->getPosition().z >= t->getBorderLowerSE().z && v->getPosition().z <= t->getBorderUpperSE().z)
        return true;
    else
        return false;
}

/**
 * @brief Octree::minDistance
 * @param v point to test
 * @param t node of the octree
 * @return the minimal distance between v and a node of the octree
 */
float Octree::minDistance(Vertex *v, Octree *t)
{
    // computing the distance between v and each of the 8 vertices of the octree cube
    float min = 10000.f;    // ~ infinity
    float dist;

    dist = Vertex::distance3(v->getPosition(), t->getBorderLowerNE());
    if (dist < min)
        min = dist;
    dist = Vertex::distance3(v->getPosition(), t->getBorderLowerNW());
    if (dist < min)
        min = dist;
    dist = Vertex::distance3(v->getPosition(), t->getBorderLowerSE());
    if (dist < min)
        min = dist;
    dist = Vertex::distance3(v->getPosition(), t->getBorderLowerSW());
    if (dist < min)
        min = dist;

    dist = Vertex::distance3(v->getPosition(), t->getBorderUpperNE());
    if (dist < min)
        min = dist;
    dist = Vertex::distance3(v->getPosition(), t->getBorderUpperNW());
    if (dist < min)
        min = dist;
    dist = Vertex::distance3(v->getPosition(), t->getBorderUpperSE());
    if (dist < min)
        min = dist;
    dist = Vertex::distance3(v->getPosition(), t->getBorderUpperSW());
    if (dist < min)
        min = dist;

    return dist;
}

/**
 * @brief Octree::setLeaf warning deletes each children
 */
void Octree::setLeaf()
{
    _isLeaf = true;

    // deleting chidren recursively by calling Octree::deleteFromNode
    Octree::deleteFromNode(_lowerNW);
    Octree::deleteFromNode(_lowerNE);
    Octree::deleteFromNode(_lowerSW);
    Octree::deleteFromNode(_lowerSE);
    Octree::deleteFromNode(_upperNW);
    Octree::deleteFromNode(_upperNE);
    Octree::deleteFromNode(_upperSW);
    Octree::deleteFromNode(_upperSE);

    _lowerNW = nullptr;
    _lowerNE = nullptr;
    _lowerSW = nullptr;
    _lowerSE = nullptr;
    _upperNW = nullptr;
    _upperNE = nullptr;
    _upperSW = nullptr;
    _upperSE = nullptr;
}

/**
 * @brief Octree::setNotLeaf warning the 8 pointers to children will still be null !
 */
void Octree::setNotLeaf()
{
    _isLeaf = false;
}

/**
 * @brief Octree::leaf
 * @return true if the current node is a leaf, false otherwise
 */
bool Octree::leaf()
{
    return _isLeaf;
}

/**
 * @brief Octree::deleteFromNode
 * @param t node to delete (all children of that node will be deleted as well)
 */
void Octree::deleteFromNode(Octree *t)
{
    if (t == nullptr)
    {
        return;
    }
    else if (t->leaf())
    {
        delete t;
    }
    else
    {
        // deleting recursively each 8 children then the current node
        Octree::deleteFromNode(t->getLowerNW());
        Octree::deleteFromNode(t->getLowerNE());
        Octree::deleteFromNode(t->getLowerSW());
        Octree::deleteFromNode(t->getLowerSE());
        Octree::deleteFromNode(t->getUpperNW());
        Octree::deleteFromNode(t->getUpperNE());
        Octree::deleteFromNode(t->getUpperSW());
        Octree::deleteFromNode(t->getUpperSE());

        delete t;
    }
}

/**
 * @brief Octree::findSpaceBorders find the space borders of the starting cube (root)
 * @param vertices list of all vertices (point cloud in that case)
 */
void Octree::findSpaceBorders(std::vector<Vertex *> &vertices)
{
    float minX = 0.f;
    float maxX = 0.f;
    float minY = 0.f;
    float maxY = 0.f;
    float minZ = 0.f;
    float maxZ = 0.f;

    std::vector<Vertex *>::iterator vertexIterator;
    glm::vec3 pos;

    // for each vertex of vertices
    for (vertexIterator = vertices.begin(); vertexIterator != vertices.end(); vertexIterator++)
    {
        pos = (*vertexIterator)->getPosition();

        // finding min & max for all XYZ coordinates
        if(pos.x < minX)
            minX = pos.x;
        if(pos.x > maxX)
            maxX = pos.x;
        if(pos.y < minY)
            minY = pos.y;
        if(pos.y > maxY)
            maxY = pos.y;
        if(pos.z < minZ)
            minZ = pos.z;
        if(pos.z > maxZ)
            maxZ = pos.z;
    }

    // debug
    std::cout << "Space borders (outer root cube) :" << std::endl;
    std::cout << "X = [" << minX << "," << maxX << "]" << std::endl;
    std::cout << "Y = [" << minY << "," << maxY << "]" << std::endl;
    std::cout << "Z = [" << minZ << "," << maxZ << "]" << std::endl;

    // setting the 8 points constraining the root cube
    _borderLowerNW = glm::vec3(minX, maxY, minZ);
    _borderLowerNE = glm::vec3(maxX, maxY, minZ);
    _borderLowerSW = glm::vec3(minX, minY, minZ);
    _borderLowerSE = glm::vec3(maxX, minY, minZ);
    _borderUpperNW = glm::vec3(minX, maxY, maxZ);
    _borderUpperNE = glm::vec3(maxX, maxY, maxZ);
    _borderUpperSW = glm::vec3(minX, minY, maxZ);
    _borderUpperSE = glm::vec3(maxX, minY, maxZ);
}

/**
 * @brief Octree::constructWithMinSize recursively construct the octree until the leaves have the given size
 * @param size minimum size the cube should reach
 * @param vertices list of all vertices (point cloud)
 */
void Octree::constructWithMinSize(float size, std::vector<Vertex *> &vertices)
{
    _minSize = size;
    _depthTest = false;
    _sizeTest = true;

    Octree::__buildOctreeNode(this, 0, vertices);
}

/**
 * @brief Octree::constructWithIterations recursively construct the octree for a depth of k
 * @param k depth of the tree (number of subdivisions of space)
 * @param vertices list of all vertices (point cloud)
 */
void Octree::constructWithIterations(int k, std::vector<Vertex *> &vertices)
{
    _depth = k;
    _depthTest = true;
    _sizeTest = false;

    Octree::__buildOctreeNode(this, 0, vertices);
}

// custom compare function for the priority queue in findKNearestNeighbours method
class Compare
{
public:
    bool operator() (std::pair<Vertex *, float> v1, std::pair<Vertex *, float> v2)
    {
        return v1.second > v2.second;
    }
};

/**
 * @brief Octree::findKNeartestNeighbours returns the k nearest neighbours of a point
 * @param ref the point of reference
 * @param k the number of neighbours to find
 * @return a vector containing the k nearest neighbours and their distances to v
 */
std::vector<std::pair<Vertex *, float>>& Octree::findKNeartestNeighbours(Vertex *ref, int k)
{
    // here we will store the neighbours
    static std::vector<std::pair<Vertex *, float>> neighbours;
    neighbours.clear();
    // here we will class the neighbours in increasing relative distance to ref
    std::priority_queue<std::pair<Vertex *, float>, std::vector<std::pair<Vertex *, float>>, Compare> pq;

    // first we need to find the subregion in which v belongs to
    Octree *t = this;

    /*std::cout << "finding " << k << " nearest neighbours of " << ref->getPosition().x
              << ", " << ref->getPosition().y << ", " << ref->getPosition().z << std::endl;*/

    while (!t->leaf())
    {
        if (Octree::isInside(ref, t->getLowerNE()))
        {
            t = t->getLowerNE();
        }
        else if (Octree::isInside(ref, t->getLowerNW()))
        {
            t = t->getLowerNW();
        }
        else if (Octree::isInside(ref, t->getLowerSE()))
        {
            t = t->getLowerSE();
        }
        else if (Octree::isInside(ref, t->getLowerSW()))
        {
            t = t->getLowerSW();
        }
        else if (Octree::isInside(ref, t->getUpperNE()))
        {
            t = t->getUpperNE();
        }
        else if (Octree::isInside(ref, t->getUpperNW()))
        {
            t = t->getUpperNW();
        }
        else if (Octree::isInside(ref, t->getUpperSE()))
        {
            t = t->getUpperSE();
        }
        else if (Octree::isInside(ref, t->getUpperSW()))
        {
            t = t->getUpperSW();
        }
        else
        {
            std::cerr << "Octree::findKNearestNeighbour error point is not iside one of the 8 children !" << std::endl;
            break;
        }
    }

    // then we place each points other than v in a priority queue according to their distances to v
    std::vector<Vertex *>::iterator it;
    std::vector<Vertex *> points = t->getPoints();
    std::pair<Vertex *, float> elem;
    float distance = 0.f;

    for (it = points.begin(); it < points.end(); it++)
    {
        distance = Vertex::distance3(ref->getPosition(), (*it)->getPosition());
        elem.first = (*it);
        elem.second = distance;

        // we shouldn't included ref himself !
        if ((*it)->getId() != ref->getId())
            pq.push(elem);
    }

    // placing the k smallest values in neighbours
    int nb = 0;
    while (nb < k && !pq.empty())
    {
        elem = pq.top();
        pq.pop();
        neighbours.push_back(elem);
        nb++;
    }

    return neighbours;
}

//-----------------------------------------------------------------------------
// Iternal method(s)
//-----------------------------------------------------------------------------

/**
 * @brief Octree::__buildOctreeNode builds the children of the current node
 * @param t node to build
 * @param depth
 * @param vertices
 */
void Octree::__buildOctreeNode(Octree *t, int depth, std::vector<Vertex *>& vertices)
{
    // halting condition is either a min size or a max depth to reach
    if (_depthTest && depth >= _depth)
    {
        t->setLeaf();
        // finding all points in the current subregion
        Octree::__findPointsInRegion(t ,vertices);
        return;
    }

    if (_sizeTest && Vertex::distance3(t->getBorderUpperNE(), t->getBorderUpperNW()) <= _minSize)
        return;

    // just dividing the region in 8 subregions (4 lower, 4 upper)
    t->setNotLeaf();

    // center of cube
    glm::vec3 center = (t->getBorderLowerSW() + t->getBorderUpperNE()) /2.f; // = middle of diagonal

    // center of each face (determined by center of diagonals)
    glm::vec3 centerUpperFace = (t->getBorderUpperNE() + t->getBorderUpperSW()) / 2.f;
    glm::vec3 centerLowerFace = (t->getBorderLowerNE() + t->getBorderLowerSW()) / 2.f;
    glm::vec3 centerNorthFace = (t->getBorderUpperNE() + t->getBorderLowerNW()) / 2.f;
    glm::vec3 centerSouthFace = (t->getBorderUpperSE() + t->getBorderLowerSW()) / 2.f;
    glm::vec3 centerWestFace = (t->getBorderUpperNW() + t->getBorderLowerSW()) / 2.f;
    glm::vec3 centerEastFace = (t->getBorderUpperNE() + t->getBorderLowerSE()) / 2.f;

    // center of each edge
    glm::vec3 centerNorthWestEdge = (t->getBorderUpperNW() + t->getBorderLowerNW()) / 2.f;
    glm::vec3 centerNorthEastEdge = (t->getBorderUpperNE() + t->getBorderLowerNE()) / 2.f;
    glm::vec3 centerSouthWestEdge = (t->getBorderUpperSW() + t->getBorderLowerSW()) / 2.f;
    glm::vec3 centerSouthEastEdge = (t->getBorderUpperSE() + t->getBorderLowerSE()) / 2.f;

    glm::vec3 centerLowerNorthEdge = (t->getBorderLowerNW() + t->getBorderLowerNE()) / 2.f;
    glm::vec3 centerLowerSouthEdge = (t->getBorderLowerSW() + t->getBorderLowerSE()) / 2.f;
    glm::vec3 centerLowerEastEdge = (t->getBorderLowerNE() + t->getBorderLowerSE()) / 2.f;
    glm::vec3 centerLowerWestEdge = (t->getBorderLowerNW() + t->getBorderLowerSW()) / 2.f;

    glm::vec3 centerUpperNorthEdge = (t->getBorderUpperNW() + t->getBorderUpperNE()) / 2.f;
    glm::vec3 centerUpperSouthEdge = (t->getBorderUpperSW() + t->getBorderUpperSE()) / 2.f;
    glm::vec3 centerUpperEastEdge = (t->getBorderUpperNE() + t->getBorderUpperSE()) / 2.f;
    glm::vec3 centerUpperWestEdge = (t->getBorderUpperNW() + t->getBorderUpperSW()) / 2.f;

    // creating the 8 subregions
    t->setUpperNW(new Octree(centerNorthWestEdge, centerNorthFace, centerWestFace, center,
                             t->getBorderUpperNW(), centerUpperNorthEdge, centerUpperWestEdge, centerUpperFace));

    t->setUpperNE(new Octree(centerNorthFace, centerNorthEastEdge, center, centerEastFace,
                             centerUpperNorthEdge, t->getBorderUpperNE(), centerUpperFace, centerUpperEastEdge));

    t->setUpperSW(new Octree(centerWestFace, center, centerSouthWestEdge, centerSouthFace,
                             centerUpperWestEdge, centerUpperFace, t->getBorderUpperSW(), centerUpperSouthEdge));

    t->setUpperSE(new Octree(center, centerEastFace, centerSouthFace, centerSouthEastEdge,
                             centerUpperFace, centerUpperEastEdge, centerUpperSouthEdge, t->getBorderUpperSE()));

    t->setLowerNW(new Octree(t->getBorderLowerNW(), centerLowerNorthEdge, centerLowerWestEdge, centerLowerFace,
                             centerNorthWestEdge, centerNorthFace, centerWestFace, center));

    t->setLowerNE(new Octree(centerLowerNorthEdge, t->getBorderLowerNE(), centerLowerFace, centerLowerEastEdge,
                             centerNorthFace, centerNorthEastEdge, center, centerEastFace));

    t->setLowerSW(new Octree(centerLowerWestEdge, centerLowerFace, t->getBorderLowerSW(), centerLowerSouthEdge,
                             centerWestFace, center, centerSouthWestEdge, centerSouthFace));

    t->setLowerSE(new Octree(centerLowerFace, centerLowerEastEdge, centerLowerSouthEdge, t->getBorderLowerSE(),
                             center, centerEastFace, centerSouthFace, centerSouthEastEdge));

    // setting the current node as the father of each 8 children
    t->getUpperNE()->setFather(t);
    t->getUpperNW()->setFather(t);
    t->getUpperSE()->setFather(t);
    t->getUpperSW()->setFather(t);
    t->getLowerNE()->setFather(t);
    t->getLowerNW()->setFather(t);
    t->getLowerSE()->setFather(t);
    t->getLowerSW()->setFather(t);

    // recurvely calling the Octree creation on its 8 children
    Octree::__buildOctreeNode(t->getLowerNW(), depth + 1, vertices);
    Octree::__buildOctreeNode(t->getLowerNE(), depth + 1, vertices);
    Octree::__buildOctreeNode(t->getLowerSW(), depth + 1, vertices);
    Octree::__buildOctreeNode(t->getLowerSE(), depth + 1, vertices);
    Octree::__buildOctreeNode(t->getUpperNW(), depth + 1, vertices);
    Octree::__buildOctreeNode(t->getUpperNE(), depth + 1, vertices);
    Octree::__buildOctreeNode(t->getUpperSW(), depth + 1, vertices);
    Octree::__buildOctreeNode(t->getUpperSE(), depth + 1, vertices);
}

/**
 * @brief Octree::__findPointsInRegion find all points belonging to a node (leaf)
 * @param t node of the Octree in which to store all points found
 * @param vertives list of all vertices (point cloud)
 */
void Octree::__findPointsInRegion(Octree *t, std::vector<Vertex *>& vertives)
{
    // should only be done for leaves
    if (!t->leaf())
        return;

    std::vector<Vertex *>::iterator vertexIterator;
    glm::vec3 pos;
    int nb = 0;

    for (vertexIterator = vertives.begin(); vertexIterator != vertives.end(); vertexIterator++)
    {
        pos = (*vertexIterator)->getPosition();

        if (pos.x >= t->getBorderLowerSW().x && pos.x <= t->getBorderLowerSE().x &&
                pos.y >= t->getBorderLowerSE().y && pos.y <= t->getBorderLowerNE().y &&
                pos.z >= t->getBorderLowerSE().z && pos.z <= t->getBorderUpperSE().z)
        {
            nb++;
            t->addPoint((*vertexIterator));
        }
    }

   //std::cout << "found " << nb << " points" << std::endl;
}
