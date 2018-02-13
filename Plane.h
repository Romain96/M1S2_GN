#ifndef PLANE_H
#define PLANE_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// glm
#include "glm/glm.hpp"
#include "glm/vec3.hpp"

/**
 * @brief The Plane class plane represented by three vectors (eigenvectors)
 */
class Plane
{
protected:
    glm::vec3 _eigenvector1;
    glm::vec3 _eigenvector2;
    glm::vec3 _eigenvector3;

public:
    //constructor(s)
    Plane(glm::vec3& ev1, glm::vec3& ev2, glm::vec3& ev3);

    // getter(s)
    glm::vec3& getEigenvector1();
    glm::vec3& getEigenvector2();
    glm::vec3& getEigenvector3();

    // setter(s)
    void setEigenvector1(glm::vec3& ev1);
    void setEigenvector2(glm::vec3& ev2);
    void setEigenvector3(glm::vec3& ev3);

    // methods
};

#endif // PLANE_H
