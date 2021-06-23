#!/usr/bin/env python
"""Determine the FGS Ideal-to-J frame alignment matrix using full solution.


Authors
-------
    Tony Sohn

Notes
-----
    This is the FGSics_to_J matrix calculation code using full solution instead of
    calculating and applying theta (dfference b/w predicted vs. measured).
    Because we are deriving a 3x3 matrix that converts a given input 3-element
    vector to an output 3-element vector, this method will NOT work for a single
    source -- there are infinite number of 3x3 matrices that satisfies such condition.
    We need at least 3 sources with (xFGS, yFGS, RA, Dec) to solve the 3x3 matrix.
    Normally, we would have >>3 sources, so the linear equation essentially
    becomes an overestimated system. See comments in the relevant part below.

"""

import pysiaf
import numpy as np
from scipy.linalg import lstsq
from scipy.spatial.transform import Rotation as R
from astropy import units as u
from astropy import constants as const

# eps = the smallest representable float
eps = np.finfo(float).eps

def fgs_to_jframe(xFGS_meas, yFGS_meas, RA, DEC, PA, qEST, SCVel):

    """Compute the FGS to J frame matrix.

    Description
    -----------
    The unit vector of the LOS of the star in the FGS1 ICS frame is given by:
    u_FGS = R_JtoFGS * R_ECItoJ * R_GStoECI * R_GSAPPARENTtoGS * u_GS
    This can be expressed as a matrix equation to derive R_JtoFGS given
    multiple inputs. This python code solves the matrix equation to derive
    R_FGStoJ.


    Input Parameters
    ----------------
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
    FGS1ics_to_J : ndarray
        FGS to J frame alignment matrix, 3x3

    """

    # Hard-coded: Rev. H of ACSK::attDETrHatFGStoJ
    # UPDATE this matrix if necessary. Correct as of May 2021.
    FGS1ics_to_J_ref = np.array([ [-0.000873067342356,  0.003603343153256,  0.999993131751562],
                                  [ 0.999757625459087, -0.021994952998128,  0.000952115833223],
                                  [ 0.021998233057671,  0.999751583464858, -0.003583259602304] ])

    mm = len(xFGS_meas)
    if mm < 3:
        sys.exit('This code requires at least 3 matched sources to derive a 3x3 matrix.')

    # Loop for each star
    for i in range(mm):

        # LOS unit vector of star in FGS1 ICS frame
        xFGS = np.deg2rad(xFGS_meas[i]/3600.)
        yFGS = np.deg2rad(yFGS_meas[i]/3600.)
        uFGS_meas = np.array([xFGS, yFGS, np.sqrt(1. - xFGS ** 2 - yFGS ** 2)])

        # GS -> ECI matrix (= GS attitude matrix defined by RA, DEC, PA)
        R_GStoECI = pysiaf.utils.rotations.attitude(0., 0., RA[i], DEC[i], PA[i])

        # ECI -> J matrix
        r = R.from_quat(qEST)
        R_ECItoJ = r.as_matrix().T # qEST provides JtoECI, so .T is required here

        # Compute the velocity aberration (VA) of the Guide Star in the ECI frame
        uu = R_GStoECI.T[0,:].T
        SCVel_mag = np.linalg.norm(SCVel)

        if SCVel_mag > eps: # If the VA correction is non-negligible
            n_vec = SCVel / SCVel_mag
        else: # Otherwise
            n_vec = np.array([0, 0, 0])

        beta = SCVel_mag / const.c.to('km/s').value
        gamma = 1. / np.sqrt(1. - beta ** 2)

        up_eci = (uu + (gamma - 1) * np.dot(n_vec, uu) * n_vec + beta * gamma * n_vec) \
                 / (gamma * (1 + beta * np.dot(n_vec, uu)))
        up_gs = np.dot(R_GStoECI.T, up_eci)

        A = -np.dot(R_ECItoJ, np.dot(R_GStoECI, up_gs)) # TBD: Figure out why negative is required here
        B =  uFGS_meas

        # Below is a simplified (alternative) way of constructing the transposed matrix AT and BT
        #if i == 0:
        #    AT = A
        #    BT = B
        #else:
        #    AT = np.vstack((AA,A))
        #    BT = np.vstack((BB,B))

        if i == 0:
            AA = np.vstack(A)
            BB = np.vstack(B)
        else:
            AA = np.hstack((AA,np.vstack(A)))
            BB = np.hstack((BB,np.vstack(B)))

    # Now that the matrix equation has been contructed, solve it using least square solver (lstsq).
    # Linear least square is generally used to solve matrix equations A dot x = B where x is unkonwn.
    # However, we have a canse of xA = B here. There is a simple trick to solve this.
    # Take the transpose of the entire system: (x dot A)^T = B^T . Note that the left handside
    # switches position when taking the transpose of matrix product, so the equation now becomes:
    # A.T dot x.T = B.T, and now we can solve this as x.T = lsqtsq(A.T, B.T). Note that we are
    # looking for x.T in the end, we don't need to transpose the outcome back to x.

    AT = AA.T
    BT = BB.T
    R_FGStoJ, _, _, _ = lstsq(AT, BT)

    return R_FGStoJ
