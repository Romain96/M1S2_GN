#ifndef CUBE_H
#define CUBE_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// GLM
#include "glm/glm.hpp"

/**
 * @brief The Cube class used in Marching cubes
 */
class Cube
{
protected:
    // a cube has 8 vertices to represent its boundaries
    // standard used :
    // U = upper | L = lower
    // N = north | S = south | E = east | W = west
    glm::vec3 _positions[8];

    // and has 8 associated values (or at least signs) for each vertex
    float _values[8];

public:
    // constructor(s)
    Cube();
    Cube(glm::vec3& lnw, glm::vec3& lne, glm::vec3& lse, glm::vec3& lsw,
         glm::vec3& unw, glm::vec3& une, glm::vec3& use, glm::vec3& usw);

    // getter(s)

    // for positions
    glm::vec3& getLNWPos();
    glm::vec3& getLNEPos();
    glm::vec3& getLSEPos();
    glm::vec3& getLSWPos();

    glm::vec3& getUNWPos();
    glm::vec3& getUNEPos();
    glm::vec3& getUSEPos();
    glm::vec3& getUSWPos();

    // for values/signs
    float getLNWVal();
    float getLNEVal();
    float getLSEVal();
    float getLSWVal();

    float getUNWVal();
    float getUNEVal();
    float getUSEVal();
    float getUSWVal();

    // setter(s)

    // for positions
    void setLNWPos(glm::vec3& lnw);
    void setLNEPos(glm::vec3& lne);
    void setLSEPos(glm::vec3& lse);
    void setLSWPos(glm::vec3& lsw);

    void setUNWPos(glm::vec3& unw);
    void setUNEPos(glm::vec3& une);
    void setUSEPos(glm::vec3& use);
    void setUSWPos(glm::vec3& usw);

    // for values/signs
    void setLNWPos(float lnw);
    void setLNEPos(float lne);
    void setLSEPos(float lse);
    void setLSWPos(float lsw);

    void setUNWPos(float unw);
    void setUNEPos(float une);
    void setUSEPos(float use);
    void setUSWPos(float usw);

    // method(s)

};

#endif // CUBE_H
