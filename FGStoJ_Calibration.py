#!/usr/bin/env python
"""Determine the FGS Ideal-to-J frame alignment matrix

Authors
-------
    Tony Sohn
    Johannes Sahlmann

Notes
-----
    This is a modified version of fgs_to_jframe_calibration.py, which itself
    was directly ported from the matlab script FGStoJ_Calibration_main.m and
    associated functions provided by NGAS (FGS to J-frame Alignment Calibration
    Tool, C. Rui, 09/18/2017)

    The modifications were applied to handle some issues including:
    - Replace q2m function with scipy.spatial.transform.Rotation
    - Allow multiple positions of sources as inputs
    - Only allow a single entry of telemetry data
    - (Later,build in function to read telemetry .csv file and take average values)
"""

import numpy as np
import pysiaf
from astropy import units as u
from astropy import constants as const
import scipy
from scipy.linalg import expm, pinv, pinv2, lstsq
from scipy.spatial.transform import Rotation as R

# eps = the smallest representable float
eps = np.finfo(float).eps

def skew(w):
    """Return skew matrix."""
    return np.array([ [    0, -w[2],  w[1]],
                      [ w[2],     0, -w[0]],
                      [-w[1],  w[0],     0] ])

def fgs_to_jframe(xFGS_meas, yFGS_meas, RA, DEC, PA, qEST, SCVel):

    """Compute the FGS to J frame matrix.

    Input Parameters
    ----------
    xFGS_meas : either float or ndarray with size of mm (mm = # of sources)
        Measured Guide Star centroid x-position in FGS Guider 1 ICS Frame
        ICS = Ideal Coordinate System (undistorted frame)
    yFGS_meas : either float or ndarray with size of mm
        Measured Guide Star centroid y-position in FGS Guider 1 ICS Frame
    RA : either float or ndarray with size of mm
        Guide Star Right Ascension, deg
    DEC : either float or ndarray with size of mm
        Guide Star Declination, deg
    PA : either float or ndarray with size of mm
        Position Angle about Guide Star, deg
    qEST : ndarray with size of 4
        Quaternion that specifies telescope attitude
    SCVel : ndarray with size of 3
        Spacecraft velocity relative to the sun in ECI km/s
    verbose : bool
        verbosity

    Returns
    -------
    FGS1ic_to_J : ndarray
        FGS to J frame alignment matrix, 3x3
    theta : ndarray

    """

    # Hard-coded: Rev. H of ACSK::attDETrHatFGStoJ
    # UPDATE this matrix if necessary. Correct as of May 2021.
    FGS1ics_to_J_ref = np.array([ [-0.000873067342356,  0.003603343153256,  0.999993131751562],
                                  [ 0.999757625459087, -0.021994952998128,  0.000952115833223],
                                  [ 0.021998233057671,  0.999751583464858, -0.003583259602304] ])

    if isinstance(xFGS_meas, list) or isinstance(xFGS_meas, np.ndarray):
        mm = len(xFGS_meas)
        print('Multiple input sources detected.')
        print('Number of the input stars: {}'.format(mm))
    else:
        print('Single input source detected.')
        mm = 1
        xFGS_meas = [xFGS_meas]
        yFGS_meas = [yFGS_meas]
        RA = [RA]
        DEC = [DEC]
        PA = [PA]

    # Loop for each star
    for i in range(mm):

        # ECI -> GS matrix (= attitude matrix of GS)
        R_GStoECI = pysiaf.utils.rotations.attitude(0., 0., RA[i], DEC[i], PA[i])

        xFGS = np.deg2rad(xFGS_meas[i]/3600.)
        yFGS = np.deg2rad(yFGS_meas[i]/3600.)
        uFGS_meas = np.array([xFGS, yFGS, np.sqrt(1. - xFGS ** 2 - yFGS ** 2)]).T

        # ECI -> J matrix
        r = R.from_quat(qEST)
        R_ECItoJ = r.as_matrix().T  #  qEST provides JtoECI, so .T is required here

        # Compute the velocity aberration of the Guide Star in the ECI frame
        uu = R_GStoECI.T[0,:].T
        vel_sc_mag = np.linalg.norm(SCVel)

        if vel_sc_mag > eps: # If the VA correction is non-negligible
            n_vec = SCVel / vel_sc_mag
        else: # Otherwise
            n_vec = np.array([0, 0, 0]).T

        beta = vel_sc_mag / const.c.to('km/s').value
        gamma = 1. / np.sqrt(1. - beta ** 2)

        up_eci = (uu + (gamma - 1) * np.dot(n_vec, uu) * n_vec + beta * gamma * n_vec) \
                 / (gamma * (1 + beta * np.dot(n_vec, uu)))
        up_gs = np.dot(R_GStoECI.T, up_eci)

        # Estimation Algorithm
        uHatFGS = np.dot(FGS1ics_to_J_ref.T, np.dot(R_ECItoJ, np.dot(R_GStoECI, up_gs)))
        du_FGS = uFGS_meas - uHatFGS # Difference b/w measurement and prediction based on old fgs1ics_to_jframe.
        A = skew(uHatFGS)

        print(uHatFGS)
        print()
        print(A)
        print()


        if i == 0:
            AA = A
            UERR = du_FGS
        else:
            AA = np.vstack((AA, A))
            UERR = np.hstack((UERR, du_FGS))

    ###
    ### IMPORTANT NOTE!
    ### While using pinv in Matlab works for solving the overdetermined system,
    ### python numpy/scipy version of pinv or pinv2 does NOT work for some reason.
    ### The solution is significantly different from the Matlab solution.
    ### Based on my quick tests, the solution is diverging.
    ### For this reason, I've looked into using some other methods, and concluded
    ### that the scipy.sparse.linalg.lsqr works closest to Matlab's results.
    ### But then, if I use atol=0 , btol=0 (or smaller than 1e-8), the solution
    ### seems to diverge and I get the same wrong solution as using pinv.
    ###
    if mm == 1:
        theta = np.dot(pinv(AA), UERR)
    else:
        theta, istop = scipy.sparse.linalg.lsqr(AA, UERR)[:2]


    # Finally, apply theta to the old matrix to obtain new matrix
    FGS1ics_to_J = np.dot(expm(-skew(theta)), FGS1ics_to_J_ref.T).T

    return FGS1ics_to_J, theta
