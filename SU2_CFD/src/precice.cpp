/*!
* \file precice.cpp
* \brief Adapter class for coupling SU2 with preCICE for FSI.
* \author Prateek Ranjan ( adapted from Rush )
*/
#include <iostream>
#include <string.h>
#include <stdio.h>

#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../include/variables/CEulerVariable.hpp"
#include "../include/solvers/CEulerSolver.hpp"
#include "../include/solvers/CFVMFlowSolverBase.inl"


#include "../include/precice.hpp"


/*---Main class description -----*/

  Precice::Precice(const std::string& preciceConfigurationFileName, int solverProcessIndex, int solverProcessSize, CConfig** config_container, CGeometry**** geometry_container, CSolver***** solver_container, CVolumetricMovement*** grid_movement)
  :coric(precice::constants::actionReadIterationCheckpoint()), cowic(precice::constants::actionWriteIterationCheckpoint()), solverProcessIndex(solverProcessIndex), solverProcessSize(solverProcessSize),solverInterface("SU2_CFD", preciceConfigurationFileName, solverProcessIndex, solverProcessSize),config_container(config_container),geometry_container(geometry_container)

{
    /* Get dimension of the problem */
    nDim = geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
  

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


/*---Class destructor -----*/  
Precice::~Precice(void)
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


void Precice::check()
{
    std::cout << " Precice being called " << std::endl;
}


double Precice::initialize()
{
  /* Check for dimensional consistency between SU2 and .xml file */

  if (solverInterface.getDimensions() != nDim)
  {
    std::cout << "Fluid Domain dimension mismatch " <<std::endl;
    std::cout << "Dimension in SU2: " << nDim << std::endl;
    std::cout << "Dimension in MDA config file: " << solverInterface.getDimensions() << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << " Dimension check complete " << std::endl;
  
  /* Check for total number of aero-elastic interfaces */
  
  if ( globalNumberWetSurfaces < 1)
  {
    std::cout << "Invalid number of aero-elastic interfaces." << std::endl;
    std::cout << "Number of aero-elastic surfaces: " << globalNumberWetSurfaces << std::endl;
    exit(EXIT_FAILURE);
  }

  else
  {

    /* The next three are scalars ( single interface ) or vector pointers (multiple interfaces) */
    meshID = new int[globalNumberWetSurfaces];
    forceID = new int[globalNumberWetSurfaces];
    displDeltaID = new int[globalNumberWetSurfaces];

    /* Get meshIDs from preCICE */
    for ( int i = 0; i < globalNumberWetSurfaces; i++)
    {
      meshID[i] = solverInterface.getMeshID("SU2_Mesh" + to_string(i));
    }
  }

  std::cout << " Mesh iD allocation complete " << std::endl;
   
  /*Determine the number of wet surfaces, that this process is working on, then loop over this number for all respective preCICE-related tasks */
  for (int i = 0; i < globalNumberWetSurfaces; i++) 
  {       
    if (config_container[ZONE_0]->GetMarker_All_TagBound(config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() + to_string(i)) == -1) 
    {
      std::cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Does not work on " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << i << endl;
    } 
    else 
    {
      localNumberWetSurfaces++;
    }
  }
  std::cout << " Finished determining number of wet surfaces " << std::endl;


  if (localNumberWetSurfaces < 1) 
  {
    std::cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Does not work on the wet surface at all." << endl;
    processWorkingOnWetSurface = false;
  }

  if (processWorkingOnWetSurface) 
  {
    std::cout << " Discretizing Aero-elastic Interface! \n";

    /*--Store the wet surface marker values in an array, which has the size equal to the number of wet surfaces actually being worked on by this process*/
    valueMarkerWet = new short[localNumberWetSurfaces];
    indexMarkerWetMappingLocalToGlobal = new short[localNumberWetSurfaces];
    int j = 0;
    for (int i = 0; i < globalNumberWetSurfaces; i++) 
    {
      if (config_container[ZONE_0]->GetMarker_All_TagBound(config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() + to_string(i)) != -1) 
      {
        valueMarkerWet[j] = config_container[ZONE_0]->GetMarker_All_TagBound(config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() + to_string(i));
        indexMarkerWetMappingLocalToGlobal[j] = i;
        j++;
      }
    }
    vertexIDs = new int*[localNumberWetSurfaces];
  }

  std::cout << " Finished setting indices " << std::endl;
   


  if (processWorkingOnWetSurface) 
  {
    vertexSize = new unsigned long[localNumberWetSurfaces];

    for (int i = 0; i < localNumberWetSurfaces; i++) 
    {
      vertexSize[i] = geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[valueMarkerWet[i]];
      std::cout << " Vertex Size: " << vertexSize[i] << std::endl;
      /*--- coordinates of all nodes at the wet surface ---*/
      double coupleNodeCoord[vertexSize[i]][nDim]; 

      /*--- variable for storing the node indices - one at the time ---*/
      unsigned long iNode;  


      //Loop over the vertices of the (each) boundary
      for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) 
      {
        //Get node number (= index) to vertex (= node) of the aeroelastic interface
        iNode = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNode();
   
       // vector = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetCoord(iNode);

        /*---Get coordinates for nodes at aero-elastic interface --*/
        for (int iDim = 0; iDim < nDim; iDim++) 
        {  
         // coupleNodeCoord[iVertex][iDim] = vector[iDim];
         coupleNodeCoord[iVertex][iDim] = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetCoord(iNode, iDim);
        }

        /* -- Print detailed surface information --*/
        //std::cout << " Value Markerwet: " << valueMarkerWet[i] << std::endl;
        std::cout << "Vertex: " << iVertex << " Node index: " << iNode << " x: " << coupleNodeCoord[iVertex][0] << " y: " << coupleNodeCoord[iVertex][1] <<" z: " << coupleNodeCoord[iVertex][2] << std::endl; 

        /* -- Extract the exture node indices vector --*/
       /// unsigned long nodeVertex[vertexSize[i]];
      //  nodeVertex[iVertex] = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNode(); /*--- Store all nodes (indices) in a vector ---*/
      //  std::cout << nodeVertex[iVertex] << std::endl;
      }


      std::cout << " Finished acquiring node IDs and coordinates" << std::endl;

      /*---------preCICE Internal Calculations -----*/
      
      //preCICE conform the coordinates of vertices (= points = nodes) at wet surface
      double coords[vertexSize[i]*nDim];
      for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) 
      {
        for (int iDim = 0; iDim < nDim; iDim++) 
        {
          coords[iVertex*nDim + iDim] = coupleNodeCoord[iVertex][iDim];
        }
      }

      //preCICE internal
      vertexIDs[i] = new int[vertexSize[i]];

      solverInterface.setMeshVertices(meshID[indexMarkerWetMappingLocalToGlobal[i]], vertexSize[i], coords, vertexIDs[i]);


      forceID[indexMarkerWetMappingLocalToGlobal[i]] = solverInterface.getDataID("Forces" + to_string(indexMarkerWetMappingLocalToGlobal[i]), meshID[indexMarkerWetMappingLocalToGlobal[i]]);


      displDeltaID[indexMarkerWetMappingLocalToGlobal[i]] = solverInterface.getDataID("DisplacementDeltas" + to_string(indexMarkerWetMappingLocalToGlobal[i]), meshID[indexMarkerWetMappingLocalToGlobal[i]]);
    }

    for (int i = 0; i < globalNumberWetSurfaces; i++) 
    {
      bool flag = false;
      for (int j = 0; j < localNumberWetSurfaces; j++) 
      {
        if (indexMarkerWetMappingLocalToGlobal[j] == i) 
        {
          flag = true;
        }
      }
      if (!flag) 
      {
        solverInterface.setMeshVertices(meshID[i], 0, NULL, NULL);
        forceID[i] = solverInterface.getDataID("Forces" + to_string(i), meshID[i]);
        displDeltaID[i] = solverInterface.getDataID("DisplacementDeltas" + to_string(i), meshID[i]);
      }
    }
  }

  /* Obtain Timestep across the interface */

  double precice_dt;

  std::cout << " Initialize the interface " << std::endl;

  precice_dt = solverInterface.initialize();

  std::cout << " Interface initialization complete " << std::endl;

  return precice_dt;

}

double Precice::advance( double computedTimestepLength )
{

   int procid = solverProcessIndex;

   if ( processWorkingOnWetSurface)
   {

     /*---Get total number of markers---*/
      unsigned short Markers = config_container[ZONE_0]->GetnMarker_All();

     /*---Identify the marker for the FSI Interface ---*/
      string FSI_NAME = config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName();

     /*--- Get the marker ID for FSI surface---*/
      short int FSI_ID = config_container[ZONE_0]->GetMarker_All_TagBound(FSI_NAME+to_string(0));

     /*---Number of vertices on FSI surface---*/

      unsigned long FSI_nVert = geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[FSI_ID];

     /* Tw-dimensional array consisting of all tractions ---*/

      std::cout << " Registering forces ..." << std::endl;

      double FSI_Trac[Markers][FSI_nVert][nDim];
     
      std::cout << " # of vertices on FSI Surface: " << FSI_nVert << std::endl;

     /*--- Begin Local Calculations here---*/

      const double *Velocity_ND, *Velocity_Real;
      double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;
      double factor;

    // Check whether the problem is viscous
      bool viscous_flow = config_container[ZONE_0]->GetViscous();

    // Parameters for the calculations
      double Pn = 0.0;
      double auxForce[3] = {1.0, 0.0, 0.0};

      unsigned short iMarker;
      unsigned long iVertex, iPoint;
      const double* iNormal;

      double Pressure_Inf = config_container[ZONE_0]->GetPressure_FreeStreamND();

      Velocity_Real = config_container[ZONE_0]->GetVelocity_FreeStream();
      Density_Real  = config_container[ZONE_0]->GetDensity_FreeStream();

      Velocity_ND = config_container[ZONE_0]->GetVelocity_FreeStreamND();
      Density_ND  = config_container[ZONE_0]->GetDensity_FreeStreamND();
    
    
      if (procid == 0)
      {
        std::cout << " P_inf " << Pressure_Inf << std::endl;
      }
    
      Velocity2_Real = GeometryToolbox::SquaredNorm(nDim, Velocity_Real);
      Velocity2_ND   = GeometryToolbox::SquaredNorm(nDim, Velocity_ND);

      factor = Density_Real * Velocity2_Real / ( Density_ND * Velocity2_ND );

      if (procid == 0)
      {
        std::cout << " factor " << factor << std::endl;
      }

   /*---Begin loop through all markers but compute only if the marker index
   matches the index of the wetsurface ---*/

      for ( int iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++)
      {
         /*--- If this is defined as a wall ---*/
        if (!config_container[ZONE_0]->GetSolid_Wall(iMarker)) continue;

        std::cout << " Pass now  at iMarker: " << iMarker << std::endl;

        // Loop over the vertices
        for (iVertex = 0; iVertex < geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker]; iVertex++) 
        {

        // Recover the point index
        iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
      // Get the normal at the vertex: this normal goes inside the fluid domain.
        iNormal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();

      /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
        if (geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint)) 
        {

        // Retrieve the values of pressure
       // Pn = base_nodes->GetPressure(iPoint);
          Pn = 0;
        //std::cout << "Pressure at "<< iPoint << " is " << Pn << std::endl;

        // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
          for (int iDim = 0; iDim < nDim; iDim++)
          {
            auxForce[iDim] = -(Pn-Pressure_Inf)*iNormal[iDim];
          }
        // Calculate tn in the fluid nodes for the viscous term
       /* if (viscous_flow) {
          su2double Viscosity = base_nodes->GetLaminarViscosity(iPoint);
          su2double Tau[3][3];
          CNumerics::ComputeStressTensor(nDim, Tau, base_nodes->GetGradient_Primitive(iPoint)+1, Viscosity);
          for (iDim = 0; iDim < nDim; iDim++) {
            auxForce[iDim] += GeometryToolbox::DotProduct(nDim, Tau[iDim], iNormal);
          }
        }

        */

          double FSI_Trac[FSI_nVert][nDim];

          if ( iMarker == FSI_ID)
          {  
          // Redimensionalize the forces
            for (int iDim = 0; iDim < nDim; iDim++) 
            {
              FSI_Trac[iVertex][iDim] = factor * auxForce[iDim];
            }
          }

          if (procid == 0)
          {
            std::cout << " Vertex Index " << iVertex << "/"<< FSI_nVert << " Traction_x: " << FSI_Trac[iVertex][0] << " Traction_y: " << FSI_Trac[iVertex][1] << " Traction_z: " << FSI_Trac[iVertex][2] << std::endl;
          }
        } 
      }
    }
  
  
  /*---Define else condition here */
  
   /*---Advance ends here ---*/
   return 0;
}

bool Precice::isCouplingOngoing()
{
  return solverInterface.isCouplingOngoing();
}

bool Precice::isActionRequired( const string& action )
{
  return solverInterface.isActionRequired(action);
}

const string& Precice::getCowic()
{
  return cowic;
}

const string& Precice::getCoric()
{
  return coric;
}

void Precice::finalize()
{
  solverInterface.finalize();
}



