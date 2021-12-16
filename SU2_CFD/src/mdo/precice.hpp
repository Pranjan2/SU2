/*!
* \file precice.hpp
* \brief Header of adapter class for coupling SU2 with preCICE for FSI.
* \author Alexander Rusch
*/


#include <string.h>
#include <stdlib.h>

#include "../../include/SU2_CFD.hpp"

#include "precice/SolverInterface.hpp"

using namespace precice;
using namespace std;

class Precice {
private:
  /* Local process and total process count */
  int solverProcessIndex, solverProcessSize;

  /* Coupling interface object */
  SolverInterface solverInterface;

  /* Fluid mesh and boundary information */
  CGeometry**** geometry_container;
  
  /* Fluid solution information */
  CSolver***** solver_container;

  /* Configuration information */
  CConfig** config_container;

  /* Fluid mesh elasticity information */
  CVolumetricMovement*** grid_movement;

  int nDim;    // Dimension of the problem
  unsigned long *vertexSize; //Number of nodes at the wet surface (for each wet surface)
  short *valueMarkerWet; //List of wet surface marker values
  int **vertexIDs;
  int *forceID;
  int *displDeltaID;
  double *forces;
  double *displacements;
  double *displacementDeltas;
  //const std::string& coric;
  //const std::string& cowic;
  bool processWorkingOnWetSurface;
  bool verbosityLevel_high;
  unsigned long globalNumberWetSurfaces; //Number of wet surfaces of the whole problem
  unsigned long localNumberWetSurfaces; //Number of wet surfaces, which this process is actually working on
  short *indexMarkerWetMappingLocalToGlobal; //Mapping relations for wet surfaces
  
  //Variables for storing the old state to reset to in case of an implicit simulation
  int nPoint; //Overall number of nodes of the problem
  int nVar; //Number of variables of the problem
  double **Coord_Saved, **Coord_n_Saved, **Coord_n1_Saved, **Coord_p1_Saved, **GridVel_Saved, ***GridVel_Grad_Saved;
  double dt_savedState;
  bool StopCalc_savedState;
  double **solution_Saved, **solution_time_n_Saved, **solution_time_n1_Saved;

  public:

  /*--------------- Override default constructor for class Precice ------*/
  Precice(const std::string& preciceConfigurationFileName, int solverProcessIndex, int solverProcessSize, CConfig** config_container, CGeometry**** geometry_container, CSolver***** solver_container, CVolumetricMovement*** grid_movement)
  {
    solverProcessIndex = solverProcessIndex;
    solverProcessSize = solverProcessSize;
    //solverInterface= precice::SolverInterface("SU2_CFD", preciceConfigurationFileName, solverProcessIndex, solverProcessSize);
    //solverInterface = precice::SolverInterface( "SU2_CFD", preciceConfigurationFileName, solverProcessIndex, solverProcessSize );

    //coric = precice::constants::actionReadIterationCheckpoint();
    //cowic = precice::constants::actionWriteIterationCheckpoint();
  //solverInterface( "SU2_CFD", preciceConfigurationFileName, solverProcessIndex, solverProcessSize ),
    /* Get dimension of the problem */
    nDim = geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
  
    /* Set the relevant containers */
    solver_container = solver_container;
    geometry_container = geometry_container;
    config_container = config_container;
    grid_movement = grid_movement;

    /* Initialize the coupling datasets to NULL */
    vertexIDs = NULL;
    forceID = NULL;
    displDeltaID = NULL;
    forces = NULL;
    displacementDeltas = NULL;

    /* Initialize the variables for implicit coupling */

    processWorkingOnWetSurface = true;
    verbosityLevel_high = config_container[ZONE_0]->GetpreCICE_VerbosityLevel_High();
    globalNumberWetSurfaces = config_container[ZONE_0]->GetpreCICE_NumberWetSurfaces();
    localNumberWetSurfaces = 0;
    valueMarkerWet = NULL;
    vertexSize = NULL;
    indexMarkerWetMappingLocalToGlobal = NULL;



    /* Get the total number of nodes in the domain */
    nPoint = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();

    /* Get the total number of variables in the domain */
    nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Initialize and allocate memory for arrays required for implicit coupling */
    Coord_Saved            = NULL;
    Coord_n_Saved          = NULL;
    Coord_n1_Saved         = NULL;
    Coord_p1_Saved         = NULL;
    GridVel_Saved          = NULL;
    GridVel_Grad_Saved     = NULL;
    dt_savedState          = 0;
    StopCalc_savedState    = false;
    solution_Saved         = NULL;
    solution_time_n_Saved  = NULL;
    solution_time_n1_Saved = NULL;

    Coord_Saved = new double*[nPoint];
    Coord_n_Saved = new double*[nPoint];
    Coord_n1_Saved = new double*[nPoint];
    Coord_p1_Saved = new double*[nPoint];
    GridVel_Saved = new double*[nPoint];
    GridVel_Grad_Saved = new double**[nPoint];
    solution_Saved = new double*[nPoint];
    solution_time_n_Saved = new double*[nPoint];
    solution_time_n1_Saved = new double*[nPoint];

    for (int iPoint = 0; iPoint < nPoint; iPoint++) 
    {
      Coord_Saved[iPoint] = new double[nDim];
      Coord_n_Saved[iPoint] = new double[nDim];
      Coord_n1_Saved[iPoint] = new double[nDim];
      Coord_p1_Saved[iPoint] = new double[nDim];
      GridVel_Saved[iPoint] = new double[nDim];
      GridVel_Grad_Saved[iPoint] = new double*[nDim];
      for (int iDim = 0; iDim < nDim; iDim++) 
      {
        GridVel_Grad_Saved[iPoint][iDim] = new double[nDim];
      }
      solution_Saved[iPoint] = new double[nVar];
      solution_time_n_Saved[iPoint] = new double[nVar];
      solution_time_n1_Saved[iPoint] = new double[nVar];
    }
  }


  
  /*--------- Override default destructor for class Precice -------*/
  ~Precice()
  {
    for (int i = 0; i < localNumberWetSurfaces; i++) 
    {
      if (vertexIDs[i] != NULL)
        {
          delete [] vertexIDs[i];
        }
    }

    if (vertexIDs != NULL) 
    {
      delete [] vertexIDs;
    }
    if (forceID != NULL) 
    {
      delete [] forceID;
    }
    if (displDeltaID != NULL) 
    {
      delete [] displDeltaID;
    }
    if (valueMarkerWet != NULL) 
    {
      delete [] valueMarkerWet;
    }
    if (vertexSize != NULL) 
    {
      delete [] vertexSize;
    }
    if (indexMarkerWetMappingLocalToGlobal != NULL) 
    {
      delete [] indexMarkerWetMappingLocalToGlobal;
    }

    for (int iPoint = 0; iPoint < nPoint; iPoint++) 
    {
      if (Coord_Saved[iPoint] != NULL) 
      {
        delete [] Coord_Saved[iPoint];
      }
      if (Coord_n_Saved[iPoint] != NULL) 
      {
        delete [] Coord_n_Saved[iPoint];
      }
      if (Coord_n1_Saved[iPoint] != NULL) 
      {
        delete [] Coord_n1_Saved[iPoint];
      }
      if (Coord_p1_Saved[iPoint] != NULL) 
      {
        delete [] Coord_p1_Saved[iPoint];
      }
      if (GridVel_Saved[iPoint] != NULL) 
      {
        delete [] GridVel_Saved[iPoint];
      }
      for (int iDim = 0; iDim < nDim; iDim++) 
      {
        if (GridVel_Grad_Saved[iPoint][iDim] != NULL) 
        {
          delete [] GridVel_Grad_Saved[iPoint][iDim];
        }
      }
      if (GridVel_Grad_Saved[iPoint] != NULL) 
      {
        delete [] GridVel_Grad_Saved[iPoint];
      }
      if (solution_Saved[iPoint] != NULL) 
      {
        delete [] solution_Saved[iPoint];
      }
      if (solution_time_n_Saved[iPoint] != NULL) 
      {
        delete [] solution_time_n_Saved[iPoint];
      }
      if (solution_time_n1_Saved[iPoint] != NULL) 
      {
        delete [] solution_time_n1_Saved[iPoint];
      }
    }
    if (Coord_Saved != NULL) 
    {
      delete [] Coord_Saved;
    }
    if (Coord_n_Saved != NULL) 
    {
      delete [] Coord_n_Saved;
    }
    if (Coord_n1_Saved != NULL)
    {
      delete [] Coord_n1_Saved;
    }
    if (Coord_p1_Saved != NULL) 
    {
      delete [] Coord_p1_Saved;
    }
    if (GridVel_Saved != NULL) 
    {
      delete [] GridVel_Saved;
    }
    if (GridVel_Grad_Saved != NULL) 
    {
      delete [] GridVel_Grad_Saved;
    }
    if (solution_Saved != NULL) 
    {
      delete [] solution_Saved;
    }
    if (solution_time_n_Saved != NULL) 
    {
      delete [] solution_time_n_Saved;
    }
    if (solution_time_n1_Saved != NULL) 
    {
      delete [] solution_time_n1_Saved;
    }
  }



  /* Method to see if this adapter file is being called from SU2 */

  void check();

  /*!
  * \brief Fully initializes preCICE to be used.
  *
  * Preconditions:
  * - configure() has been called successfully.
  *
  * Postconditions:
  * - Communication to the coupling partner/s is setup.
  * - Meshes are created from geometries or sent/received to/from coupling partners.
  * - If the solver is not starting the simulation, coupling data is received
  *   from the coupling partner's first computation.
  * - The length limitation of the first solver timestep is computed and returned.
  *
  * \return Maximum length of first timestep to be computed by the solver.
  */
  double initialize();


/*!
  * \brief Advances preCICE after the solver has computed one timestep.
  *
  * Preconditions:
  * - initialize() has been called successfully.
  * - The solver calling advance() has computed one timestep.
  * - The solver has read and written all coupling data.
  * - The solver has the length of the timestep used available to be passed to
  *   preCICE.
  * - finalize() has not yet been called.
  *
  * Postconditions:
  * - Coupling data values specified to be exchanged in the configuration are
  *   exchanged with coupling partner/s.
  * - Coupling scheme state (computed time, computed timesteps, ...) is updated.
  * - For staggered coupling schemes, the coupling partner has computed one
  *   timestep/iteration with the coupling data written by this participant
  *   available.
  * - The coupling state is printed.
  * - Meshes with data are exported to files if specified in the configuration.
  * - The length limitation of the next solver timestep is computed and returned.
  *
  * \param[in] computedTimestepLength - Length of timestep computed by solver.
  * \return Maximum length of next timestep to be computed by solver.
  */
  double advance( double computedTimestepLength );

  /*!
  * \brief Checks wether preCICE coupling is still ongoing or not.
  *
  * \return True if coupling is still ongoing or false if coupling has ended.
  */
  bool isCouplingOngoing();

  /*!
  * \brief Checks wether an action (read or write iteration checkpoint) is required or not.
  *
  * \param[in] action - Reference to the strings coric or cowic, which determine wether read (coric) or write (cowic) action is checked.
  * \return True if an action is required or false if no action is required.
  */
  bool isActionRequired( const string& action );

  /*!
  * \brief Get reference to the string cowic which determines wether an iteration checkpoint needs to be written or not.
  *
  * \return The aforementioned reference to the cowic string.
  */
  const string& getCowic();

  /*!
  * \brief Get reference to the string coric which determines wether an iteration checkpoint needs to be read or not.
  *
  * \return The aforementioned reference to the coric string.
  */
  const string& getCoric();

  /*!
  * \brief Save the current state before the next SU2 run
  *
  * \param[in] StopCalc - Pointer to the bool StopCalc, which determines wether SU2 finishes computation after the current run or not. This variable must be passed so that its value can be stored and reloaded later.
  * \param[in] double - Pointer to the double dt, which is the current SU2 time step size. This variable must be passed so that its value can be stored and reloaded later.
  */
  void saveOldState( bool *StopCalc, double *dt );

  /*!
  * \brief Reload the state before the current SU2 run
  *
  * \param[in] StopCalc - Pointer to the bool StopCalc, which determines wether SU2 finishes computation after the current run or not. This variable must be passed so that its value can be reloaded from the saved state.
  * \param[in] double - Pointer to the double dt, which is the current time step size. This variable must be passed so that its value can be reloaded from the saved state.
  */
  void reloadOldState( bool *StopCalc, double *dt );

  /*!
  * \brief Finalizes preCICE.
  *
  * Preconditions:
  * - initialize() has been called successfully.
  * - isCouplingOngoing() has returned false.
  *
  * Postconditions:
  * - Communication channels are closed.
  * - Meshes and data are deallocated.
  */
  void finalize();


};

