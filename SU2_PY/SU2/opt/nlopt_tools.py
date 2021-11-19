#!/usr/bin/env python

## \file nlopt_tools.py
#  \brief tools for interfacing with nlopt
#  \author Prateek. Ranjan
#  \version 7.2.0 "Blackbird"

import nlopt
from numpy import * 
import sys
from .. import eval as su2eval

def nlopt_optim(project,x0,xb,its,accu):
    """ Runs nlopt implementation with an SU2 project

    Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10

    """

    #handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []

    #Function handle for objective and objective sensitivity
    #func           = obj_f
    #fprime         = obj_df

    # Number of design varibales
    dv_size = project.config['DEFINITION_DV']['SIZE']
    n_dv = sum( dv_size)
    project.n_dv = n_dv

    # Set the initial guess to zero
    if not x0: x0 = [0.0]*n_dv

    # Pre-scale the DV
    dv_scales = project.config['DEFINITION_DV']['SCALE']
    k = 0
    for i, dv_scl in enumerate(dv_scales):
        for j in range(dv_size[i]):
            x0[k] =x0[k]/dv_scl
            k = k + 1

    #Scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = []
    for this_obj in obj.keys():
        obj_scale = obj_scale + [obj[this_obj]['SCALE']]
        
    #Scale accuracy
    eps = 1.0e-20

    lb = [0.0]*n_dv
    ub = [0.0]*n_dv
    for i in range(n_dv):
        lb[i] = xb[i][0]
        ub[i] = xb[i][1]


    print('Check complete \n')    

        












