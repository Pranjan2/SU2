/*!
* \file precice.hpp
* \brief Header of adapter class for coupling SU2 with preCICE for FSI.
* \author Alexander Rusch
*/


#include <string.h>
#include <stdlib.h>

#include "../../include/SU2_CFD.hpp"

class Precice {
private:
 
  int procid, numprocs;

  public:

  //Precice( const string& preciceConfigurationFileName, int solverProcessIndex, int solverProcessSize, CGeometry*** geometry_container, CSolver**** solver_container, CConfig** config_container, CVolumetricMovement** grid_movement );
  // Override default constructor for class Precice
  Precice(const std::string& preciceConfigurationFileName, int solverProcessIndex, int solverProcessSize, CGeometry**** geometry_container)
  {
    procid = solverProcessIndex;
    numprocs = solverProcessSize;
  }

  /* Method to see if this adapter file is being called from SU2 */

  void check();


};

