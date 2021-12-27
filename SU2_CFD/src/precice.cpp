/*!
* \file precice.cpp
* \brief Adapter class for coupling SU2 with preCICE for FSI.
* \author Prateek Ranjan ( adapted from Rush )
*/
#include <iostream>
#include <string.h>
#include <stdio.h>

#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"


#include "../include/precice.hpp"


/*---Main class description -----*/

  Precice::Precice(const std::string& preciceConfigurationFileName, int solverProcessIndex, int solverProcessSize, CConfig** config_container, CGeometry**** geometry_container, CSolver***** solver_container, CVolumetricMovement*** grid_movement)
  :coric(precice::constants::actionReadIterationCheckpoint()), cowic(precice::constants::actionWriteIterationCheckpoint()), solverProcessIndex(solverProcessIndex), solverProcessSize(solverProcessSize),solverInterface("SU2_CFD", preciceConfigurationFileName, solverProcessIndex, solverProcessSize),config_container(config_container),geometry_container(geometry_container)

{
    /* Get dimension of the problem */
    nDim = geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
  
    /* Set the relevant containers */
    solver_container = solver_container;
  //  geometry_container = geometry_container;
   // config_container = config_container;
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

  //double Pn = 0.0;

  //Pn = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->nodes->GetPressure();
  //std::cout << " Pressure at node: " << Pn << std::endl;
  return precice_dt;

}

 double Precice::advance( double computedTimestepLength )
 {
   if (processWorkingOnWetSurface) 
   {
     /*-- If the current process if working on the Aero-elastic surface, flow variables ---*/

     /*-- Get simulation information --*/
     //bool incompressible = (config_container[ZONE_0]->ENUM_REGIME::GetKind_Regime() == INCOMPRESSIBLE);
     bool viscous_flow = config_container[ZONE_0]->GetViscous();

     /*-- Compute variables for re-dimensionalizing forces (ND := Non-Dimensional) ---*/
     double* Velocity_Real = config_container[ZONE_0]->GetVelocity_FreeStream();
     double Density_Real = config_container[ZONE_0]->GetDensity_FreeStream();
     double* Velocity_ND = config_container[ZONE_0]->GetVelocity_FreeStreamND();
     double Density_ND = config_container[ZONE_0]->GetDensity_FreeStreamND();
     double Velocity2_Real = 0.0;  /*--- denotes squared real velocity ---*/
     double Velocity2_ND = 0.0;  /*--- denotes squared non-dimensional velocity ---*/

     /* -- Compute squared values --*/
     for (int iDim = 0; iDim < nDim; iDim++)
     {
       Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
       Velocity2_ND += Velocity_ND[iDim]*Velocity_ND[iDim];
     }

     /* -- Compute factors for re-dimensionalizing forces --*/
     double factorForces = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);

     /* -- Loop over each local aero-elastic surface to perform operations --*/
     for (int i = 0; i < localNumberWetSurfaces; i++)
     {
       /* -- Define temporary variables for calculations --*/
       unsigned long nodeVertex[vertexSize[i]];
       double normalsVertex[vertexSize[i]][nDim];
       double normalsVertex_Unit[vertexSize[i]][nDim];
       double Area;
       double Pn = 0.0;  /*--- denotes pressure at a node ---*/
       double Pinf = 0.0;  /*--- denotes environmental (farfield) pressure ---*/
       //double** Grad_PrimVar = NULL; /*--- denotes (u.A. velocity) gradients needed for computation of viscous forces ---*/
       //double Viscosity = 0.0;
       //double Tau[3][3];
       //double TauElem[3];
       double forces_su2[vertexSize[i]][nDim];  /*--- forces will be stored such, before converting to simple array ---*/

       /* -- Loop over vertices of coupled boundary -- */
       for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) 
       {
         //Get node number (= index) to vertex (= node)
         nodeVertex[iVertex] = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNode(); /*--- Store all nodes (indices) in a vector ---*/
         // Get normal vector
         for (int iDim = 0; iDim < nDim; iDim++)
         {
           normalsVertex[iVertex][iDim] = (geometry_container[ZONE_0][INST_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNormal())[iDim];
         }
         // Unit normals
         Area = 0.0;
         for (int iDim = 0; iDim < nDim; iDim++) 
         {
           Area += normalsVertex[iVertex][iDim]*normalsVertex[iVertex][iDim];
         }
         Area = sqrt(Area);

         for (int iDim = 0; iDim < nDim; iDim++) 
         {
           normalsVertex_Unit[iVertex][iDim] = normalsVertex[iVertex][iDim]/Area;
         }

         /* -- Get the values of pressure and viscosity --*/

         // Nodes pointer in SU2 corresponds to all nodes in the domain. Extract a set of nodes ( interface nodes) from nodes
//         double *IFNodes;
//         IFNodes = nodes[nodeVertex[iVertex]];

         /* -- Compute the pressure at each node at aero-elastic interface -- */
           Pn = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetPressure(nodeVertex[iVertex]);
         //Pn = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->nodes[nodeVertex[iVertex]]->GetPressure(nodeVertex[iVertex]);
//         Pn = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->IFNodes->GetPressure();
           Pinf = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetPressure_Inf();

       //  if (viscous_flow)
      //   {
      //     Viscosity = solver_container[ZONE_0][MESH_0][INST_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(nodeVertex[iVertex]);
      //     Grad_PrimVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetGradient_Primitive(nodeVertex[iVertex]);
            // Viscosity = solver_container[ZONE_0][MESH_0][INST_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(nodeVertex[iVertex]);
//           Grad_PrimVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->IFNodes->GetGradient_Primitive();
//           Viscosity = solver_container[ZONE_0][MESH_0][INST_0][FLOW_SOL]->IFNodes->GetLaminarViscosity();
    //     }

          /* -- Compute the forces_su2 in the nodes for the inviscid term -- */

         for (int iDim = 0; iDim < nDim; iDim++) 
         {
           forces_su2[iVertex][iDim] = -(Pn-Pinf)*normalsVertex[iVertex][iDim];
         }

         /* -- Compute the viscous contributions --*/

         if (viscous_flow)
         {
           double Viscosity = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(nodeVertex[iVertex]);
           double Tau[3][3];
           CNumerics::ComputeStressTensor(nDim, Tau, solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetGradient_Primitive(nodeVertex[iVertex])+1, Viscosity);
           //for (int iDim = 0; iDim < nDim; iDim++) 
          // {
          //  forces_su2[iVertex][iDim] += GeometryToolbox::DotProduct(nDim, Tau[iDim],normalsVertex[iVertex][iDim]);
          // }
         }

         /* -- Compute the forces_su2 in the nodes for the viscous term -- */

        // if (viscous_flow)
        // {
           // Divergence of the velocity
        //   double div_vel = 0.0;
        //   for (int iDim = 0; iDim < nDim; iDim++)
        //   {
        //     div_vel += Grad_PrimVar[iDim+1][iDim];
        //   }

           //if (incompressible)
          // {
          //   div_vel = 0.0;  /*--- incompressible flow is divergence-free ---*/
          // }
        //   for (int iDim = 0; iDim < nDim; iDim++)
        //   {
        //     for (int jDim = 0 ; jDim < nDim; jDim++) 
        //     {
               // Dirac delta
        //       double Delta = 0.0;
        //       if (iDim == jDim)
        //       {
        //         Delta = 1.0;
        //       }

        //       // Viscous stress
        //       Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -2/3*Viscosity*div_vel*Delta;
               // Add Viscous component in the forces_su2 vector --> Units of force (non-dimensional).
        //       forces_su2[iVertex][iDim] += Tau[iDim][jDim]*normalsVertex[iVertex][jDim];
        //     }
        //   }
       //  }

         // Rescale forces_su2 to SI units
         for (int iDim = 0; iDim < nDim; iDim++)
         {
           forces_su2[iVertex][iDim] = forces_su2[iVertex][iDim]*factorForces;
         }
       }

       /* -- Convert the force vector in row-major format --*/
       forces = new double[vertexSize[i]*nDim];

       for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++)
       {
         for (int iDim = 0; iDim < nDim; iDim++) 
         {
           //Do not write forces for duplicate nodes! -> Check wether the color of the node matches the MPI-rank of this process. Only write forces, if node originally belongs to this process.
           if (geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetColor(nodeVertex[iVertex]) == solverProcessIndex) 
           {
             forces[iVertex*nDim + iDim] = forces_su2[iVertex][iDim];
           }
           else
           {
             forces[iVertex*nDim + iDim] = 0;
           }
         }
       }

       solverInterface.writeBlockVectorData(forceID[indexMarkerWetMappingLocalToGlobal[i]], vertexSize[i], vertexIDs[i], forces);

       if (forces != NULL)
       {
         delete [] forces;
       }
     }

     /* -- Advance solver interface --*/

     double max_precice_dt;
     max_precice_dt = solverInterface.advance( computedTimestepLength );

     for (int i = 0; i < localNumberWetSurfaces; i++) 
     {
       //4. Read displacements/displacementDeltas
       double displacementDeltas_su2[vertexSize[i]][nDim]; /*--- displacementDeltas will be stored such, before converting to simple array ---*/
       displacementDeltas = new double[vertexSize[i]*nDim];
       solverInterface.readBlockVectorData(displDeltaID[indexMarkerWetMappingLocalToGlobal[i]], vertexSize[i], vertexIDs[i], displacementDeltas);

         //5. Set displacements/displacementDeltas
    
       //convert displacementDeltas into displacementDeltas_su2
       for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) 
       {
         for (int iDim = 0; iDim < nDim; iDim++)
         {
           displacementDeltas_su2[iVertex][iDim] = displacementDeltas[iVertex*nDim + iDim];
         }
       }
       if (displacementDeltas != NULL) 
       {
         delete [] displacementDeltas;
       }

       //Set change of coordinates (i.e. displacementDeltas)
       for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) 
       {
                  geometry_container[ZONE_0][INST_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->SetVarCoord(displacementDeltas_su2[iVertex]);
       }
     }
    
     return max_precice_dt;
   }


   else
   {
     double max_precice_dt;
     max_precice_dt = solverInterface.advance( computedTimestepLength );
     return max_precice_dt;
   }


}



void Precice::finalize()
{
  solverInterface.finalize();
}


