#!/usr/bin/env python
"""Determine the FGS2 to FGS1 matrix given two sets of SIAF parameters.

Authors
-------
    Tony Sohn
    
Notes
-----
    This is the FGS2 to FGS1 matrix calculation code based on the 
    pysiaf.utils.tools.jwst_fgs_to_fgs_matrix. The pysiaf tool is not 
    able to take inputs, but instead only provides the FGS2_to_FGS1 
    matrix as stored in the LoadsPRD.

"""

import sys
import numpy as np
from pysiaf.utils.rotations import rotate

def fgs2tofgs1_matrix(fgs1_siaf_params, fgs2_siaf_params):
    """

    Parameters
    ----------
    fgs1_siaf_params: array (or tuple) or V2Ref, V3Ref, V3IdlYAngle for FGS1
    fgs2_siaf_params: array (or tuple) or V2Ref, V3Ref, V3IdlYAngle for FGS2

    Returns
    -------
    R_BtoA: ndarray, FGS2 to FGS1 rotation matrix

    """

    FGS1V2 = fgs1_siaf_params[0]
    FGS1V3 = fgs1_siaf_params[1]
    FGS1pa = fgs1_siaf_params[2]

    FGS2V2 = fgs2_siaf_params[0]
    FGS2V3 = fgs2_siaf_params[1]
    FGS2pa = fgs2_siaf_params[2]

    # Form RA = R3.R2.R1
    R1 = rotate(1, -FGS1pa)
    R2 = rotate(2, -FGS1V3 / 3600.0)
    R3 = rotate(3,  FGS1V2 / 3600.0)
    RA = np.dot(R2, R1)
    RA = np.dot(R3, RA)

    # Form RB = R6.R5.R4
    R4 = rotate(1, -FGS2pa)
    R5 = rotate(2, -FGS2V3 / 3600.0)
    R6 = rotate(3,  FGS2V2 / 3600.0)
    RB = np.dot(R5, R4)
    RB = np.dot(R6, RB)

    R_BtoA = np.dot(np.transpose(RA), RB)

    return R_BtoA

