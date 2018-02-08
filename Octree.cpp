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

#include "Vertex.h"
#include "Octree.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief Octree::Octree default constructor (useful for root)
 */
Octree::Octree() :
    _borderLowerNW(),
    _borderLowerNE(glm::vec3(0.f, 0.f, 0.f)),
    _borderLowerSW(glm::vec3(0.f, 0.f, 0.f)),
    _borderLowerSE(glm::vec3(0.f, 0.f, 0.f)),
    _borderUpperNW(glm::vec3(0.f, 0.f, 0.f)),
    _borderUpperNE(glm::vec3(0.f, 0.f, 0.f)),
    _borderUpperSW(glm::vec3(0.f, 0.f, 0.f)),
    _borderUpperSE(glm::vec3(0.f, 0.f, 0.f)),
    _isLeaf(true),
    _lowerNW(nullptr),
    _lowerNE(nullptr),
    _lowerSW(nullptr),
    _lowerSE(nullptr),
    _upperNW(nullptr),
    _upperNE(nullptr),
    _upperSW(nullptr),
    _upperSE(nullptr)
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
    _lowerNW(nullptr),
    _lowerNE(nullptr),
    _lowerSW(nullptr),
    _lowerSE(nullptr),
    _upperNW(nullptr),
    _upperNE(nullptr),
    _upperSW(nullptr),
    _upperSE(nullptr)
{

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
    _borderUpperSW = v;
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
