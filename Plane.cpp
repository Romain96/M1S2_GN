/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// glm
#include "glm/glm.hpp"
#include "glm/vec3.hpp"

#include "Plane.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

Plane::Plane() :
    _eigenvector1(glm::vec3(0.f)),
    _eigenvector2(glm::vec3(0.f)),
    _eigenvector3(glm::vec3(0.f))
{
    // nothing
}

/**
 * @brief Plane::Plane constructs a new plane with ACP computed eigenvectors
 * @param ev1 eigenvector 1
 * @param ev2 eigenvector 2
 * @param ev3 eigenvector 3
 */
Plane::Plane(glm::vec3 &ev1, glm::vec3 &ev2, glm::vec3 &ev3) :
    _eigenvector1(ev1),
    _eigenvector2(ev2),
    _eigenvector3(ev3)
{
    // nothing
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Plane::getEigenvector1
 * @return the first eigenvector
 */
glm::vec3& Plane::getEigenvector1()
{
    return _eigenvector1;
}

/**
 * @brief Plane::getEigenvector2
 * @return the second eigenvector
 */
glm::vec3& Plane::getEigenvector2()
{
    return _eigenvector2;
}

/**
 * @brief Plane::getEigenvector3
 * @return the third eigenvector
 */
glm::vec3& Plane::getEigenvector3()
{
    return _eigenvector3;
}

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Plane::setEigenvector1
 * @param ev1 new first eigenvector
 */
void Plane::setEigenvector1(glm::vec3 &ev1)
{
    _eigenvector1 = ev1;
}

/**
 * @brief Plane::setEigenvector2
 * @param ev2 new second eigenvector
 */
void Plane::setEigenvector2(glm::vec3 &ev2)
{
    _eigenvector2 = ev2;
}

/**
 * @brief Plane::setEigenvector3
 * @param ev3 new third eigenvector
 */
void Plane::setEigenvector3(glm::vec3 &ev3)
{
    _eigenvector3 = ev3;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------
