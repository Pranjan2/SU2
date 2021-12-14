/*!
* \file precice.cpp
* \brief Adapter class for coupling SU2 with preCICE for FSI.
* \author Alexander Rusch
*/
#include <iostream>
#include <string.h>
#include <stdio.h>
#include "precice.hpp"

//cice::Precice( const string& preciceConfigurationFileName, int solverProcessIndex, int solverProcessSize, CGeometry**** geometry_container, CSolver**** solver_container, CConfig** config_container, CVolumetricMovement** grid_movement )
//Precice::Precice( int solverProcessIndex, int solverProcessSize):
//(
//procid = solverProcessIndex;
//numprocs = solverProcessSize;
//)

void Precice::check()
{
    std::cout << " Precice being called " << std::endl;
}