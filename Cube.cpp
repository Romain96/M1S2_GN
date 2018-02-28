/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// STL
#include <limits>

// GLM
#include "glm/glm.hpp"

#include "Cube.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief Cube::Cube unparametrized Cube constructor
 */
Cube::Cube()
{
    // creating the array of positions
    _positions = new glm::vec3[8];

    // filling with (0,0,0)
    glm::vec3 pos(0.f, 0.f, 0.f);
    for (unsigned int i = 0; i < 8; i++)
        _positions[i] = pos;

    // creating the array of values
    _values = new float[8];

    // filling with "undefined" (infinity)
    for (unsigned int i = 0; i < 8; i++)
        _positions[i] = std::numeric_limits<float>::infinity();
}

/**
 * @brief Cube::Cube parametrized Cube constructor (with positions)
 * @param lnw lower north west position
 * @param lne lower north east position
 * @param lse lower south east position
 * @param lsw lower south west position
 * @param unw upper north west position
 * @param une upper north east position
 * @param use upper south east position
 * @param usw upper south west position
 */
Cube::Cube(glm::vec3 &lnw, glm::vec3 &lne, glm::vec3 &lse, glm::vec3 &lsw,
           glm::vec3 &unw, glm::vec3 &une, glm::vec3 &use, glm::vec3 &usw) :
{
    // creating the array of positions
    _positions = new glm::vec3[8];

    // filling with values
    _positions[0] = lnw;
    _positions[1] = lne;
    _positions[2] = lse;
    _positions[3] = lsw;

    _positions[4] = unw;
    _positions[5] = une;
    _positions[6] = use;
    _positions[7] = usw;

    // creating the array of values
    _values = new float[8];

    // filling with "undefined" (infinity)
    for (unsigned int i = 0; i < 8; i++)
        _positions[i] = std::numeric_limits<float>::infinity();
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Cube::getLNWPos
 * @return the lower north west vertex position
 */
glm::vec3& Cube::getLNWPos()
{
    return _positions[0];
}

/**
 * @brief Cube::getLNEPos
 * @return the lower north east vertex position
 */
glm::vec3& Cube::getLNEPos()
{
    return _positions[1];
}

/**
 * @brief Cube::getLSEPos
 * @return the lower south east vertex position
 */
glm::vec3& Cube::getLSEPos()
{
    return _positions[2];
}

/**
 * @brief Cube::getLSWPos
 * @return the lower south west vertex position
 */
glm::vec3& Cube::getLSWPos()
{
    return _positions[3];
}

/**
 * @brief Cube::getUNWPos
 * @return the upper north west vertex position
 */
glm::vec3& Cube::getUNWPos()
{
    return _positions[4];
}

/**
 * @brief Cube::getUNEPos
 * @return the upper north east vertex position
 */
glm::vec3& Cube::getUNEPos()
{
    return _positions[5];
}

/**
 * @brief Cube::getUSEPos
 * @return the upper south east vertex position
 */
glm::vec3& Cube::getUSEPos()
{
    return _positions[6];
}

/**
 * @brief Cube::getUSWPos
 * @return the upper south west vertex position
 */
glm::vec3& Cube::getUSWPos()
{
    return _positions[7];
}

//-------------------------------------

/**
 * @brief Cube::getLNWVal
 * @return the lower north west vertex value/sign
 */
float Cube::getLNWVal()
{
    return _values[0];
}

/**
 * @brief Cube::getLNEVal
 * @return the lower north east vertex value/sign
 */
float Cube::getLNEVal()
{
    return _values[1];
}

/**
 * @brief Cube::getLSEVal
 * @return the lower south east vertex value/sign
 */
float Cube::getLSEVal()
{
    return _values[2];
}

/**
 * @brief Cube::getLSWVal
 * @return the lower south west vertex value/sign
 */
float Cube::getLSWVal()
{
    return _values[3];
}

/**
 * @brief Cube::getUNWVal
 * @return the upper north west vertex value/sign
 */
float Cube::getUNWVal()
{
    return _values[4];
}

/**
 * @brief Cube::getUNEVal
 * @return the upper north east vertex value/sign
 */
float Cube::getUNEVal()
{
    return _values[5];
}

/**
 * @brief Cube::getUSEVal
 * @return the upper south east vertex value/sign
 */
float Cube::getUSEVal()
{
    return _values[6];
}

/**
 * @brief Cube::getUSWVal
 * @return the upper south west vertex value/sign
 */
float Cube::getUSWVal()
{
    return _values[7];
}

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Cube::setLNWPos
 * @param lnw new lower north west vertex position
 */
void Cube::setLNWPos(glm::vec3 &lnw)
{
    _positions[0] = lnw;
}

/**
 * @brief Cube::setLNEPos
 * @param lne new lower north east vertex position
 */
void Cube::setLNEPos(glm::vec3 &lne)
{
    _positions[1] = lne;
}

/**
 * @brief Cube::setLSEPos
 * @param lse new lower south east vertex position
 */
void Cube::setLSEPos(glm::vec3 &lse)
{
    _positions[2] = lse;
}

/**
 * @brief Cube::setLSWPos
 * @param lsw new lower south west vertex position
 */
void Cube::setLSWPos(glm::vec3 &lsw)
{
    _positions[3] = lsw;
}

/**
 * @brief Cube::setUNWPos
 * @param unw new upper north west vertex position
 */
void Cube::setUNWPos(glm::vec3 &unw)
{
    _positions[4] = unw;
}

/**
 * @brief Cube::setUNEPos
 * @param une new upper north east vertex position
 */
void Cube::setUNEPos(glm::vec3 &une)
{
    _positions[5] = une;
}

/**
 * @brief Cube::setUSEPos
 * @param use new upper south east vertex position
 */
void Cube::setUSEPos(glm::vec3 &use)
{
    _positions[6] = use;
}

/**
 * @brief Cube::setUSWPos
 * @param usw new upper south west vertex position
 */
void Cube::setUSWPos(glm::vec3 &usw)
{
    _positions[7] = usw;
}

//-------------------------------------


/**
 * @brief Cube::setLNWVal
 * @param lnw new lower north west vertex value/sign
 */
void Cube::setLNWVal(float lnw)
{
    _values[0] = lnw;
}

/**
 * @brief Cube::setLNEVal
 * @param lne new lower north east vertex value/sign
 */
void Cube::setLNEVal(float lne)
{
    _values[1] = lne;
}

/**
 * @brief Cube::setLSEVal
 * @param lse new lower south east vertex value/sign
 */
void Cube::setLSEVal(float lse)
{
    _values[2] = lse;
}

/**
 * @brief Cube::setLSWVal
 * @param lsw new lower south west vertex value/sign
 */
void Cube::setLSWVal(float lsw)
{
    _values[3] = lsw;
}

/**
 * @brief Cube::setUNWVal
 * @param unw new upper north west vertex value/sign
 */
void Cube::setUNWVal(float unw)
{
    _values[4] = unw;
}

/**
 * @brief Cube::setUNEVal
 * @param une new upper north east vertex value/sign
 */
void Cube::setUNEVal(float une)
{
    _values[5] = une;
}

/**
 * @brief Cube::setUSEVal
 * @param use new upper south east vertex value/sign
 */
void Cube::setUSEVal(float use)
{
    _values[6] = use;
}

/**
 * @brief Cube::setUSWVal
 * @param usw new upper south west vertex value/sign
 */
void Cube::setUSWVal(float usw)
{
    _values[7] = usw;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------
