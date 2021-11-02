import numpy as np
import pysiaf
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
np.set_printoptions(precision=6, suppress=True, sign='+')

def verify(FGStoJ_old, FGStoJ_new):

    Ma = FGStoJ_old
    Mb = FGStoJ_new

    # First check whether each input matrix is orthogonal
    #print("*"*100)
    #print("")
    #np.dot(Ma,Ma.T)
    #np.dot(Mb,Mb.T)

    print()
    print("*"*100)
    print("Matrix Verification Report")
    print("*"*100)
    print()
#    print("(+) J2      - Rotation about J3 axis")
#    print("(+) J3      - Rotation about J2 axis")
#    print("(+) J3angle - Rotation about J1 axis")
#    print()

    # Overall rotation b/w the two matrices
    Mab = np.dot(Ma.T,Mb)
    rot_ab = np.rad2deg(np.arccos((np.trace(Mab)-1)/2))

    # Rotational components in each axis
    j1a = np.degrees(np.arctan2(-Ma[2,0],Ma[2,1])) # This is equivalent to "V3IdlYAngle"
    j2a = np.degrees(np.arcsin(Ma[2,2]))           # This is equivalent to "V3Ref"
    j3a = np.degrees(np.arctan2(Ma[1,2],Ma[0,2]))  # This is equivalent to "V2Ref"

    j1b = np.degrees(np.arctan2(-Mb[2,0],Mb[2,1]))
    j2b = np.degrees(np.arcsin(Mb[2,2]))
    j3b = np.degrees(np.arctan2(Mb[1,2],Mb[0,2]))


    print("----------------------")
    print(" Input Matrix 1 (Old) ")
    print("----------------------")
    print(Ma)

    print("")
    print("Matrix above implies that the J-frame are rotated around J1, J2, J3 axes w.r.t. FGS1 ICS of:")
    print("J1: {0:9.6f} deg".format(j1a))
    print("J2: {0:9.6f} deg".format(j2a))
    print("J3: {0:9.6f} deg".format(j3a))
    print("")

    print("----------------------")
    print(" Input Matrix 2 (New) ")
    print("----------------------")
    print(Mb)

    print("")
    print("Matrix above implies that the J-frame are rotated around J1, J2, J3 axes w.r.t. FGS1 ICS of:")
    print("J1: {0:9.6f} deg".format(j1b))
    print("J2: {0:9.6f} deg".format(j2b))
    print("J3: {0:9.6f} deg".format(j3b))
    print("")

    print("--------------------")
    print(" Matrix Comparisons ")
    print("--------------------")
    print("Overall angular separation between Matrices 1 and 2: {0:9.6f} deg".format(rot_ab))
    print("Offsets in each axis:")
    print("J1 offset: {0:9.6f} deg".format(j1b-j1a))
    print("J2 offset: {0:9.6f} deg".format(j2b-j2a))
    print("J3 offset: {0:9.6f} deg".format(j3b-j3a))
    print("**NOTE: offsets in each axis are approximate values, calculated from taking the differences from above.")


def verify_with_unit_vectors(FGStoJ_old, FGStoJ_new):

    # This is adopted from Colin Cox's notebook, but will not likely be used.

    Ma = FGStoJ_old
    Mb = FGStoJ_new

    F = []  # Prepare to make a list of vectors
    F.append(pysiaf.utils.rotations.unit(        0,       0))
    F.append(pysiaf.utils.rotations.unit(-1.0/60.0, -1.0/60))
    F.append(pysiaf.utils.rotations.unit( 1.0/60.0, -1.0/60))
    F.append(pysiaf.utils.rotations.unit( 1.0/60.0,  1.0/60))
    F.append(pysiaf.utils.rotations.unit(-1.0/60.0,  1.0/60))

    J2a = np.zeros(5)
    J3a = np.zeros(5)
    J2b = np.zeros(5)
    J3b = np.zeros(5)
    for v in range(5):
        print('\n', v)
        print('Unit Vector = ',F[v])
        Ja = np.dot(Ma, F[v])
        Jb = np.dot(Mb, F[v])
        print('Ja',Ja)
        print('Jb', Jb)
        scalar = np.dot(Ja, Jb)
        print(scalar)
        angle = np.arccos(scalar)
        print('Angular separation %10.6f degrees' %np.degrees(angle))
        J2a[v] = np.degrees(np.arctan2(Ja[0], Ja[1]))
        J3a[v] = np.degrees(np.arcsin(Ja[2]))
        print('J2J3a %10.6f %10.6f' %(J2a[v], J3a[v]))
        J2b[v] = np.degrees(np.arctan2(Jb[0], Jb[1]))
        J3b[v] = np.degrees(np.arcsin(Jb[2]))
        print('J2J3b %10.6f %10.6f' %(J2b[v], J3b[v]))
        print('Shift %10.6f %10.6f degrees' %(J2b[v]-J2a[v], J3b[v]-J3a[v]))

    J2ac = J2a-J2a.mean()
    J3ac = J3a-J3a.mean()
    J2bc = J2b-J2b.mean()
    J3bc = J3b-J3b.mean()
    print('\nCentralized values')
    print('     J2ac        J3ac        J2bc        J3bc      J2bc-J2ac   J3bc-J3ac')
    for i in range(5):
        print(6*'%12.3e' %(J2ac[i], J3ac[i], J2bc[i], J3bc[i], J2bc[i]-J2ac[i], J3bc[i]-J3ac[i]))
    maxLength = (np.hypot(J2ac-J2bc, J3ac-J3bc)).max()
    print('\nMaximum displacement %8.3e degrees' %maxLength)

    sxx = (J2ac*J2ac).sum()
    syy = (J3ac*J3ac).sum()
    sux = (J2ac*J2bc).sum()
    svx = (J3bc*J2ac).sum()
    suy = (J2bc*J3ac).sum()
    svy = (J3bc*J3ac).sum()

    scale = ((suy-svx)**2 + (sux+svy)**2)/(sxx+syy)
    theta = np.arctan2(suy-svx, sux+svy)
    print('Relative scale %10.6f' %scale)
    print('Rotation %10.6f degrees'%np.degrees(theta))
