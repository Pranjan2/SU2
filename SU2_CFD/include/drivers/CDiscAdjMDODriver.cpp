/*!
 * \file CDiscAdjMDODriver.cpp
 * \brief The main subroutines for driving coupled MDO adjoints
 * \author P. Ranjan
 * \version 7.2.0 "Columbia"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/drivers/CDiscAdjMDODriver.hpp"
#include "../../include/output/tools/CWindowingTools.hpp"
#include "../../include/output/COutputFactory.hpp"
#include "../../include/output/COutputLegacy.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/iteration/CIterationFactory.hpp"
#include "../../include/iteration/CTurboIteration.hpp"
#include "../../../Common/include/toolboxes/CQuasiNewtonInvLeastSquares.hpp"

CDiscAdjMDODriver::CDiscAdjMDODriver(char* confFile,
                                                   unsigned short val_nZone,
                                                   SU2_Comm MPICommunicator) : CSinglezoneDriver(confFile,
                                                                                                 val_nZone,
                                                                                                 MPICommunicator) 
{

    /*--- Store the number of internal iterations that will be run by the adjoint solver ---*/
    nAdjoint_Iter = config_container[ZONE_0]->GetnInner_Iter();

    /*--- Store the pointers ---*/
    config      = config_container[ZONE_0];
    iteration   = iteration_container[ZONE_0][INST_0];
    solver      = solver_container[ZONE_0][INST_0][MESH_0];
    numerics    = numerics_container[ZONE_0][INST_0][MESH_0];
    geometry    = geometry_container[ZONE_0][INST_0][MESH_0];
    integration = integration_container[ZONE_0][INST_0];

    /*--- Store the recording state ---*/
    RecordingState = RECORDING::CLEAR_INDICES;

    switch (config->GetKind_Solver()) 
    {

        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
        case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:

        if (rank == MASTER_NODE)
        {
            cout << "Direct MDA iteration for Euler/RANS equations." << endl;
        }

        /*--- Create an iteration instance for the current main solver ---*/
        direct_iteration = CIterationFactory::CreateIteration(EULER, config);

        if (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) 
        {
            /*--- Create an output class for the current COMPRESSIBLE main solver ---*/
            direct_output = COutputFactory::CreateOutput(EULER, config, nDim);
        }
        else 
        { 
            /*--- Create an output class for the current INCOMPRESSIBLE main solver ---*/
            direct_output =  COutputFactory::CreateOutput(INC_EULER, config, nDim); 
        }

        /*---Define fluid state as MAIN variables---*/
        MainVariables = RECORDING::SOLUTION_VARIABLES;
    
        if (config->GetDeform_Mesh()) 
        {
            SecondaryVariables = RECORDING::MESH_DEFORM;
        }

        else 
        { 
        /*---Declare mesh coordinates as secondary variables---*/
        SecondaryVariables = RECORDING::MESH_COORDS;
        }

        /*---Define the main solver---*/
        /*---Position of the continuous adjoint flow solution in the solver container array---*/
        MainSolver = ADJFLOW_SOL;
        break;
    }

    if ( rank == MASTER_NODE)
    {
        std::cout << "Solver-specific initialization complete!" <<std::endl;
    }

    direct_output->PreprocessHistoryOutput(config, false);
}

/*---Class destructor---*/
CDiscAdjMDODriver::~CDiscAdjMDODriver(void) 
{

  delete direct_iteration;
  delete direct_output;

}

void CDiscAdjMDODriver::Preprocess(unsigned long TimeIter) 
{

    if ( rank == MASTER_NODE)
    {
      std::cout << "Calling Preprocess() from MDO driver()"<<std::endl;
    }

    /*---Get the time counter---*/
    config_container[ZONE_0]->SetTimeIter(TimeIter);

    /*--- Preprocess the adjoint iteration ---*/
    iteration->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
                        solver_container, numerics_container, config_container,
                        surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- For the adjoint iteration we need the derivatives of the iteration function with
   *--- respect to the conservative variables. Since these derivatives do not change in the steady state case
   *--- we only have to record if the current recording is different from the main variables. ---*/

  if (RecordingState != MainVariables)
  {
    /*---Record the main conservative variables for the steady-state equation---*/
    /*---Initialize the tape, identify input, one-step direct analysis, identify output, close tape---*/
    MainRecording();
  }

}

void CDiscAdjMDODriver::Run() 
{
  
  if ( rank == MASTER_NODE)
  {
    std::cout << "Calling Run() from MDO driver"<<std::endl;
  }
  CQuasiNewtonInvLeastSquares<passivedouble> fixPtCorrector;

  if (config->GetnQuasiNewtonSamples() > 1) 
  {
    fixPtCorrector.resize(config->GetnQuasiNewtonSamples(),
                          geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(),
                          GetTotalNumberOfVariables(ZONE_0,true),
                          geometry_container[ZONE_0][INST_0][MESH_0]->GetnPointDomain());

    if (TimeIter != 0)
    { 
      GetAllSolutions(ZONE_0, true, fixPtCorrector);
    }
  }

  /*---Begin the adjoint-iteration loop---*/
    for (auto Adjoint_Iter = 0ul; Adjoint_Iter < nAdjoint_Iter; Adjoint_Iter++) 
    {

        /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
        *--- of the previous iteration. The values are passed to the AD tool.
        *--- Issues with iteration number should be dealt with once the output structure is in place. ---*/

        config->SetInnerIter(Adjoint_Iter);

        iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

        /*--- Initialize the adjoint of the objective function with 1.0. ---*/
        SetAdj_ObjFunction();

        /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/
        AD::ComputeAdjoint();


        /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
        iteration->IterateDiscAdj(geometry_container, solver_container,
                              config_container, ZONE_0, INST_0, false);

        /*--- Monitor the pseudo-time ---*/
        StopCalc = iteration->Monitor(output_container[ZONE_0], integration_container, geometry_container,
                                  solver_container, numerics_container, config_container,
                                  surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

        /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/
        AD::ClearAdjoints();

        /*--- Output files for steady state simulations. ---*/
        if (!config->GetTime_Domain()) 
        {
        iteration->Output(output_container[ZONE_0], geometry_container, solver_container,
                        config_container, Adjoint_Iter, false, ZONE_0, INST_0);
        }

        if (StopCalc) break;

        /*--- Correct the solution with the quasi-Newton approach, for current iteration ---*/
        if (fixPtCorrector.size()) 
        {
        /*---Stores the adjoint variabels at each node of the mesh ---*/  
        GetAllSolutions(ZONE_0, true, fixPtCorrector.FPresult());
        SetAllSolutions(ZONE_0, true, fixPtCorrector.compute());
        }
    }
}


void CDiscAdjMDODriver::Postprocess() 
{
  switch(config->GetKind_Solver())
  {
    case DISC_ADJ_EULER :     case DISC_ADJ_NAVIER_STOKES :     case DISC_ADJ_RANS :
    case DISC_ADJ_INC_EULER : case DISC_ADJ_INC_NAVIER_STOKES : case DISC_ADJ_INC_RANS :
    case DISC_ADJ_HEAT :

      /*--- Compute the geometrical sensitivities ---*/
      SecondaryRecording();
      break;
  }

}

/*---SetRecording() is called by MainRecording()---*/
void CDiscAdjMDODriver::SetRecording(RECORDING kind_recording)
{
  /*---Clean tape and adjoints---*/
  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution---*/
  /*---This calls multiple cascading functions in CDiscIteration and CDiscFluidSolver*cpps---*/
  /*---This function should know what is the kind_recording---*/
  iteration->SetRecording(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  /*---Enable recording and register input of the iteration --- */
  if (kind_recording != RECORDING::CLEAR_INDICES)
  {
    /*---Begin tape recording on tape---*/
    AD::StartRecording();

    if (rank == MASTER_NODE && kind_recording == MainVariables) 
    {
      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << "Direct iteration to store the primal computational graph on tape." << endl;
      cout << "Compute residuals to check the convergence of the direct problem." << endl;
    }

    /*---Identify the input variables---*/
    /*---Register solution + Register variables---*/
    iteration->RegisterInput(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);
    /*---Next step from here is forward analysis---*/
  }

  /*--- Set the dependencies of the iteration ---*/
  iteration->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_0,
                             INST_0, kind_recording);

  /*--- Do one iteration of the direct solver ---*/
  if ( rank == MASTER_NODE)
  {
    std::cout << "Running one iteration of the solver in AD driver " << std::endl;
  }
  /*---kind_recording : clear_indices, solution_variables, mesh_coords, mesh_Deform, solution+mesh*/
    DirectRun(kind_recording);
  
  /*--- Identify the kind of recording the tape currently holds ---*/
  RecordingState = kind_recording;

  /*--- Register Output of the iteration ---*/
  iteration->RegisterOutput(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Extract the objective function and register it as an output --- */
  SetObjFunction();

    /*---Print AD performance metrics---*/
  if (kind_recording != RECORDING::CLEAR_INDICES && config_container[ZONE_0]->GetWrt_AD_Statistics()) {
    if (rank == MASTER_NODE) AD::PrintStatistics();
#ifdef CODI_REVERSE_TYPE
    if (size > SINGLE_NODE) {
      su2double myMem = AD::getGlobalTape().getTapeValues().getUsedMemorySize(), totMem = 0.0;
      SU2_MPI::Allreduce(&myMem, &totMem, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      if (rank == MASTER_NODE) {
        cout << "MPI\n";
        cout << "-------------------------------------\n";
        cout << "  Total memory used      :  " << totMem << " MB\n";
        cout << "-------------------------------------\n" << endl;
      }
    }
#endif
  }
  /*---Stop tape recording---*/
  AD::StopRecording();

}

/*---Seed the adjoint value of the objective functions---*/
void CDiscAdjMDODriver::SetAdj_ObjFunction()
{
  const auto IterAvg_Obj = config->GetIter_Avg_Objective();
  su2double seeding = 1.0;

  CWindowingTools windowEvaluator = CWindowingTools();
  /*---Set seeding for unsteady adjoint calculation----*/
  if (config->GetTime_Marching() != TIME_MARCHING::STEADY)
  {
    if (TimeIter < IterAvg_Obj){
      /*--- Default behavior (in case no specific window is chosen) is to use Square-Windowing, i.e. the numerator equals 1.0 ---*/
      seeding = windowEvaluator.GetWndWeight(config->GetKindWindow(),TimeIter, IterAvg_Obj-1)/ (static_cast<su2double>(IterAvg_Obj));
    }
    else 
    {
      /*---Set seeding for steady adjoint calculation---*/  
      seeding = 0.0;
    }
  }

   /*---Set the seeding values---*/
  if (rank == MASTER_NODE)
  {
    SU2_TYPE::SetDerivative(ObjFunc, SU2_TYPE::GetValue(seeding));
  } 
  else 
  {
    SU2_TYPE::SetDerivative(ObjFunc, 0.0);
  }
}

/*---Sets the objective function (Called by Recording())---*/
void CDiscAdjMDODriver::SetObjFunction()
{

  bool heat         = (config->GetWeakly_Coupled_Heat());
  bool turbo        = (config->GetBoolTurbomachinery());
  ObjFunc = 0.0;

  direct_output->SetHistory_Output(geometry, solver, config,
                                   config->GetTimeIter(),
                                   config->GetOuterIter(),
                                   config->GetInnerIter());


  /*--- Specific scalar objective functions ---*/
  switch (config->GetKind_Solver()) 
  {
    case DISC_ADJ_INC_EULER:       case DISC_ADJ_INC_NAVIER_STOKES:      case DISC_ADJ_INC_RANS:
    case DISC_ADJ_EULER:           case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
    case DISC_ADJ_FEM_EULER:       case DISC_ADJ_FEM_NS:                 case DISC_ADJ_FEM_RANS:

    /*---Set the total combo objective function
    solver[FLOW_SOL]->SetTotal_ComboObj(0.0);

    /*--- Surface based obj. function ---*/

    solver[FLOW_SOL]->Evaluate_ObjFunc(config);
    ObjFunc += solver[FLOW_SOL]->GetTotal_ComboObj();
    if (heat)
    {
      if (config->GetKind_ObjFunc() == TOTAL_HEATFLUX) 
      {
        ObjFunc += solver[HEAT_SOL]->GetTotal_HeatFlux();
      }
      else if (config->GetKind_ObjFunc() == AVG_TEMPERATURE) 
      {
        ObjFunc += solver[HEAT_SOL]->GetTotal_AvgTemperature();
      }
    }

    /*--- This calls to be moved to a generic framework at a next stage         ---*/
    /*--- Some things that are currently hacked into output must be reorganized ---*/
    if (turbo)
    {

      solver[FLOW_SOL]->SetTotal_ComboObj(0.0);
      output_legacy->ComputeTurboPerformance(solver[FLOW_SOL], geometry, config);

      unsigned short nMarkerTurboPerf = config->GetnMarker_TurboPerformance();
      unsigned short nSpanSections = config->GetnSpanWiseSections();

      switch (config_container[ZONE_0]->GetKind_ObjFunc())
      {
      case ENTROPY_GENERATION:
        solver[FLOW_SOL]->AddTotal_ComboObj(output_legacy->GetEntropyGen(nMarkerTurboPerf-1, nSpanSections));
        break;
      case FLOW_ANGLE_OUT:
        solver[FLOW_SOL]->AddTotal_ComboObj(output_legacy->GetFlowAngleOut(nMarkerTurboPerf-1, nSpanSections));
        break;
      case MASS_FLOW_IN:
        solver[FLOW_SOL]->AddTotal_ComboObj(output_legacy->GetMassFlowIn(nMarkerTurboPerf-1, nSpanSections));
        break;
      default:
        break;
      }

      ObjFunc = solver[FLOW_SOL]->GetTotal_ComboObj();

    }
    break;
  }

  if (rank == MASTER_NODE)
  {
    AD::RegisterOutput(ObjFunc);
  }

}

void CDiscAdjMDODriver::DirectRun(RECORDING kind_recording)
{

  /*--- Mesh movement ---*/
  direct_iteration->SetMesh_Deformation(geometry_container[ZONE_0][INST_0], solver, numerics, config, kind_recording);

  /*--- Zone preprocessing ---*/
  direct_iteration->Preprocess(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Iterate the direct solver ---*/
  direct_iteration->Iterate(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Postprocess the direct solver ---*/
  direct_iteration->Postprocess(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Print the direct residual to screen ---*/
  Print_DirectResidual(kind_recording);

}

void CDiscAdjMDODriver::Print_DirectResidual(RECORDING kind_recording)
{

  /*--- Print the residuals of the direct iteration that we just recorded ---*/
  /*--- This routine should be moved to the output, once the new structure is in place ---*/
  if ((rank == MASTER_NODE) && (kind_recording == MainVariables))
  {
    cout << "Recording computational graph with respect to the fluid state variables" <<endl;
    switch (config->GetKind_Solver()) 
    {

        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
        case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
        case DISC_ADJ_FEM_EULER : case DISC_ADJ_FEM_NS : case DISC_ADJ_FEM_RANS :
        cout << "log10[U(0)]: "   << log10(solver[FLOW_SOL]->GetRes_RMS(0))
           << ", log10[U(1)]: " << log10(solver[FLOW_SOL]->GetRes_RMS(1))
           << ", log10[U(2)]: " << log10(solver[FLOW_SOL]->GetRes_RMS(2)) << "." << endl;
        cout << "log10[U(3)]: " << log10(solver[FLOW_SOL]->GetRes_RMS(3));
        if (geometry->GetnDim() == 3) cout << ", log10[U(4)]: " << log10(solver[FLOW_SOL]->GetRes_RMS(4));
            cout << "." << endl;
        if ( config->GetKind_Turb_Model() != NONE && !config->GetFrozen_Visc_Disc()) 
        {
            cout << "log10[Turb(0)]: "   << log10(solver[TURB_SOL]->GetRes_RMS(0));
            if (solver[TURB_SOL]->GetnVar() > 1) cout << ", log10[Turb(1)]: " << log10(solver[TURB_SOL]->GetRes_RMS(1));
            cout << "." << endl;
        }
      break;

    }

    cout << "-------------------------------------------------------------------------" << endl << endl;
  }

  else if ((rank == MASTER_NODE) && (kind_recording == SecondaryVariables) && (SecondaryVariables != RECORDING::CLEAR_INDICES))
  {
    cout << endl << "Recording the computational graph with respect to the ";
    switch (SecondaryVariables)
    {
      case RECORDING::MESH_COORDS: cout << "mesh coordinates." << endl;    break;
      default:                     cout << "secondary variables." << endl; break;
    }
  }

}

void CDiscAdjMDODriver::MainRecording()
{

  /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with NONE
   *    as argument ensures that all information from a previous recording is removed. ---*/

  SetRecording(RECORDING::CLEAR_INDICES);

  /*--- Store the computational graph of one direct iteration with the solution variables as input. ---*/
  SetRecording(MainVariables);

}

void CDiscAdjMDODriver::SecondaryRecording()
{

  /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with NONE
   *    as argument ensures that all information from a previous recording is removed. ---*/

  SetRecording(RECORDING::CLEAR_INDICES);

  /*--- Store the computational graph of one direct iteration with the secondary variables as input. ---*/
  SetRecording(SecondaryVariables);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the current iteration. The values are passed to the AD tool. ---*/
  /*--- Here, the converged adjoint of the output from mainRecording are set ---*/

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/
  SetAdj_ObjFunction();

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/
  AD::ComputeAdjoint();

  /*--- Extract the computed sensitivity values ---*/
  if (SecondaryVariables == RECORDING::MESH_COORDS) 
  {
    solver[MainSolver]->SetSensitivity(geometry, config);
  }
  else 
  {
    // MESH_DEFORM
    solver[ADJMESH_SOL]->SetSensitivity(geometry, config, solver[MainSolver]);
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

}













