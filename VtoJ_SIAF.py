#!/usr/bin/env python
"""Determine VtoJ SIAF parameters to be used for adding an aperture for FGS1 SIAF.

Authors
-------
    Tony Sohn

Notes
-----

"""

import pysiaf
import sys
import numpy as np
from scipy.spatial.transform import Rotation as R

def get_vtoj_siaf(FGS1_V2Ref, FGS1_V3Ref, FGS1_V3IdlYAngle, FGStoJ_matrix):
    """
     Description
    -----------

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
    R_VtoJ: Full 3x3 VtoJ matrix
    """

    R_FGStoJ = FGStoJ_matrix
    R_JtoFGS = R_FGStoJ.T

    # Check that the input FGStoJ matrix has a form close to:
    # [0 0 1]
    # [1 0 0]
    # [0 1 0]
    if (R_FGStoJ[0,2]<0.95) or (R_FGStoJ[1,0]<0.95) or (R_FGStoJ[2,1]<0.95):
        sys.exit('Input FGS1ics_to_J matrix appears to be incorrect.')

    ya = np.radians(FGS1_V3IdlYAngle)
    v2 = np.radians(FGS1_V2Ref/3600.)
    v3 = np.radians(FGS1_V3Ref/3600.)

    a11 = -np.cos(ya)*np.sin(v2) + np.sin(ya)*np.sin(v3)*np.cos(v2)
    a12 =  np.cos(ya)*np.cos(v2) + np.sin(ya)*np.sin(v3)*np.sin(v2)
    a13 = -np.sin(ya)*np.cos(v3)
    a21 = -np.sin(ya)*np.sin(v2) - np.cos(ya)*np.sin(v3)*np.cos(v2)
    a22 =  np.sin(ya)*np.cos(v2) - np.cos(ya)*np.sin(v3)*np.sin(v2)
    a23 =  np.cos(ya)*np.cos(v3)
    a31 =  np.cos(v3)*np.cos(v2)
    a32 =  np.cos(v3)*np.sin(v2)
    a33 =  np.sin(v3)

    R_VtoFGS = np.array([[a11, a12, a13],
                         [a21, a22, a23],
                         [a31, a32, a33]])
    R_FGStoV = R_VtoFGS.T
    print('R_FGStoV = \n', R_FGStoV)

    R_VtoJ = np.dot(R_FGStoJ, R_VtoFGS)
    R_JtoV = R_VtoJ.T

    v2ref   = 3600.*np.degrees(np.arctan2(R_VtoJ[0,1],R_VtoJ[0,0]))
    v3ref   = 3600.*np.degrees(np.arcsin(R_VtoJ[0,2]))
    v3angle = np.degrees(np.arctan2(-R_VtoJ[1,2],R_VtoJ[2,2]))

    return v2ref, v3ref, v3angle, R_VtoJ

def get_vtoj_siaf_rotation(FGS1_V2Ref, FGS1_V3Ref, FGS1_V3IdlYAngle, FGStoJ_matrix):
    """
    Description
    -----------
    This version is using the scipy.spatial.transform.Rotation
    NOTE: Results are slightly different compared to get_vtoj_siaf function above.

    Warning
    -------
    Do NOT use this function for deriving final values of VtoJ.
    Use the get_gtoj_siaf function above instead. This function is provided as
    an alternative version which uses scipy.spatial.transform.Rotation.

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
