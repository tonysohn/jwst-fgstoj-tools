#!/usr/bin/env python
"""Determine the FGS Ideal-to-J frame alignment matrix.


Authors
-------
    Tony Sohn

Notes
-----
    This is the FGS1ics_to_J matrix calculation code.
"""

import sys
import pysiaf
import numpy as np
from scipy.linalg import expm, lstsq, pinv
from scipy.sparse.linalg import lsqr
from scipy.spatial.transform import Rotation as R
from astropy import units as u
from astropy import constants as const
from astropy.time import Time

# The smallest representable float
eps = np.finfo(float).eps

def skew(w):
    """Return skew matrix."""
    return np.array([ [    0, -w[2],  w[1]],
                      [ w[2],     0, -w[0]],
                      [-w[1],  w[0],     0] ])


def read_telemetry(telemetry_csvfile, average=True, obstime_mjd=None):
    """Read in telemetry file and calculate either average or interpolated
       qEST and SCVel parameters. By default, this routine will return the
       average qEST and SCVel of the entire csv file assuming that the WOW
       Telemetry Report were extracted using correct date/time duration.
       If user wants the interpolated results, the observation time in MJD
       (obstime_mjd) must be provided. This can be easily obtained from the
       FITS file e.g., as follows:

       >>> from jwst import datamodels
       >>> f = datamodels.open('<name of FITS file>')
       >>> t1 = f.meta.exposure.start_time
       >>> t2 = f.meta.exposure.end_time
       >>> obstime_mjd = 0.5*(t1+t2)

    Input Parameters
    ----------------
    telemetry_csvfile: string with the name of the CSV csvfile
    average: argument for whether to average the entire time duration
    obstime_mjd: optional parameter
                 time of observation in modified juliand date (MJD)

    Returns
    -------
    qEST: 4-element array with the quaternion
    SCVel: 3-element array with the spacecraft velocity (VX, VY, VZ) in km/s
    """

    # Read in the telemetry file
    t = ascii.read(telemetry_csvfile, format='csv')

    # Add MJD column to the ascii table
    time_array = t['Primary Time'].data
    new_time_array = []
    for string in time_array:
        new_string = string.replace("/", "-")
        new_time_array.append(new_string)
    time_iso = Time(new_time_array, format='iso', scale='utc')
    time_mjd = time_iso
    time_mjd.format='mjd'
    t['MJD'] = time_mjd.value
    t.sort('time_mjd') # Sort table by time (although, it already should be)

    # Telemetry parameters we're interested in. In principle, we can change this
    # to get average or interplated results for any set of parameters.
    mnemonics = ['SA_ZKFQ1EST' , 'SA_ZKFQ2EST' , 'SA_ZKFQ3EST' , 'SA_ZKFQ4EST' , # Quaternion
                 'SA_ZSUNVLSCX', 'SA_ZSUNVLSCY', 'SA_ZSUNVLSCZ']                 # Spacecraft velocity
    results = []

    if average:
        print("%12s  %11s  %11s" % ("Parameter", "Mean", "Std. Dev."))
        for m in mnemonics:
            mask = t['Telemetry Mnemonic'] == m
            avg = np.mean(t[mask]['EU Value'].data)
            std = np.std(t[mask]['EU Value'].data)
            results.append(avg)
            print("%12s  %11.6f  %11.6f" % (m, avg, std))
        print("")
    else:
        if obstime_mjd < min(time_mjd) or obstime_mjd > max(time_mjd):
            sys.exit('ERROR: Primary Time in telemetry file must bracket observation time')

        for m in mnemonics:
            mask = t['Telemetry Mnemonic'] == m
            xx = t[mask]['MJD'].data
            yy = t[mask]['EU Value'].data
            zz = np.interp(obstime_mjd, xx, yy)
            results.append(zz)

    qEST  = results[:4]
    SCVel = results[4:]

    return qEST, SCVel


def fgstoj_offset(xFGS_meas, yFGS_meas, RA, DEC, PA, qEST, SCVel):

    """Compute the FGS to J frame matrix using the offset method provided by NGAS.
       This requires a pre-loaded FGS-to-J matrix since the algorithm calculates
       the offset resulting from the old matrix vs. new measurements.

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

    Returns
    -------
    FGS1ic_to_J : ndarray
        FGS to J frame alignment matrix, 3x3
    theta : ndarray

    """

    #!!!
    #!!! Hard-coded: Rev. H of ACSK::attDETrHatFGStoJ
    #!!! UPDATE this matrix if necessary. Correct as of Aug 2021.
    #!!!
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

        if i == 0:
            AA = A
            UERR = du_FGS
        else:
            AA = np.vstack((AA, A))
            UERR = np.hstack((UERR, du_FGS))

    if mm == 1:
        theta = np.dot(pinv(AA), UERR)
    else:
        theta, istop = lsqr(AA, UERR)[:2]

    # Finally, apply theta to the old matrix to obtain new matrix
    FGS1ics_to_J = np.dot(expm(-skew(theta)), FGS1ics_to_J_ref.T).T

    return FGS1ics_to_J, theta


def fgstoj_full(xFGS_meas, yFGS_meas, RA, DEC, PA, qEST, SCVel):

    """Compute the FGS to J frame matrix using the full solution.

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

    Returns
    -------
    FGS1ics_to_J : ndarray
        FGS to J frame alignment matrix, 3x3

    """

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

# def fgstoj(xFGS_meas, yFGS_meas, RA, DEC, PA, telemetry_csvfile, obstime_mjd):
#     """Compute the FGS to J frame matrix.
#
#     Description
#     -----------
#     This is a wrapper function dealing with the general case. This function
#     takes the telemetry csv file as input to calculate the average qEST and
#     SCVel at the the of observation (via interpolation). It then determines
#     whether to use the "full" or "offset" solution based on how many
#     input sources there are and returns the FGS1ics_to_J 3x3 array.
#
#     Input Parameters
#     ----------------
#     xFGS_meas : either float or ndarray with size of mm (mm = # of sources)
#         Measured Guide Star centroid x-position in FGS Guider 1 ICS Frame
#         ICS = Ideal Coordinate System (undistorted frame)
#     yFGS_meas : either float or ndarray with size of mm
#         Measured Guide Star centroid y-position in FGS Guider 1 ICS Frame
#     RA : either float or ndarray with size of mm
#         Guide Star Right Ascension, deg
#     DEC : either float or ndarray with size of mm
#         Guide Star Declination, deg
#     PA : either float or ndarray with size of mm
#         Position Angle about Guide Star, deg
#     telemetry_csvfile : name of the telemetry csv file.
#     obstime_mjd : observed date & time in modified julian date (MJD).
#                   this can be extracted from the meta.exposure.mid_time keyword.
#
#     Returns
#     -------
#     FGS1ics_to_J : ndarray
#         FGS to J frame alignment matrix, 3x3
#     """
#
#     qEST, SCVel = read_telemetry(telemetry_csvfile, obstime_mjd)
#
#     mm = len(xFGS_meas)
#
#     if mm < 3:
#         fgs_to_j, theta = fgstoj_offset(xFGS_meas, yFGS_meas, RA, DEC, PA, qEST, SCVel)
#     else:
#         fgs_to_j = fgstoj_full(xFGS_meas, yFGS_meas, RA, DEC, PA, qEST, SCVel)
#
#     return fgs_to_j
