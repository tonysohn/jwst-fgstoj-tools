#!/usr/bin/env python
"""Determine VtoJ SIAF parameters to be used for adding an aperture for FGS1 SIAF.

Authors
-------
    Tony Sohn

Notes
-----
    WARNING: Do NOT use this routine for deriving final values of VtoJ. Use vtoj_siaf.py.
    This is an alternative version of vtoj_siaf.py using scipy.spatial.transform.Rotation.
    The resulting SIAF values and matrix are slightly different from those of vtoj_siaf.py.
"""

import pysiaf
import sys
import numpy as np
from scipy.spatial.transform import Rotation as R

def get_vtoj_siaf_rotation(FGS1_V2Ref, FGS1_V3Ref, FGS1_V3IdlYAngle, FGStoJ_matrix):
    """
     Description
    -----------
    This version is using the scipy.spatial.transform.Rotation
    NOTE: Results are slightly different compared to vtoj_siaf.py

    Input Parameters
    ----------------
    FGS1_V2Ref: V2Ref SIAF parameter for FGS1 in arcsec
    FGS1_V3Ref: V3Ref SIAF parameter for FGS1 in arcsec
    FGS1_V3IdlYAngle: V3IdlYAngle SIAF parameter for FGS1 in degrees
    FGStoJ_matrix: FGStoJ matrix in 3x3 form

    Returns
    -------
    v2ref: V2Ref SIAF parameter corrsponding to pseudo VtoJ aperture in arcsec
    v3ref: V3Ref SIAF parameter corrsponding to pseudo VtoJ aperture in arcsec
    v3angle: V3IdlYAngle SIAF parameter corresponding to pseudo VtoJ aperture in degrees
    VtoJ_matrix: Full 3x3 VtoJ matrix
    """


    R_FGStoJ = FGStoJ_matrix
    R_JtoFGS = R_FGStoJ.T

    aa =  FGS1_V3Ref/3600.
    bb = -FGS1_V2Ref/3600.
    cc =  FGS1_V3IdlYAngle
    r = R.from_euler('xyz', [aa,bb,cc], degrees=True)
    zxy = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    R_FGStoV = np.dot(zxy, r.as_matrix().T)
    R_VtoFGS = R_FGStoV.T
    print('R_FGStoV = \n', R_FGStoV)

    # Check that the input FGStoJ matrix has a form close to:
    # [0 0 1]
    # [1 0 0]
    # [0 1 0]
    if (R_FGStoJ[0,2]<0.95) or (R_FGStoJ[1,0]<0.95) or (R_FGStoJ[2,1]<0.95):
        sys.exit('Input FGS1ics_to_J matrix appears to be incorrect.')

    R_VtoJ = np.dot(R_FGStoJ, R_VtoFGS)
    R_JtoV = R_VtoJ.T

    rr = R.from_matrix(R_VtoJ)
    abc = rr.as_euler('xyz', degrees=True)
    v2ref   = -abc[2]*3600.
    v3ref   =  abc[1]*3600.
    v3angle =  abc[0]

    return v2ref, v3ref, v3angle, R_VtoJ
