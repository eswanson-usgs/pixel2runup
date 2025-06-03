##%% Plot TWL estimates on raw camera image
##% also plots measured runup values on camera image
##% plots cross-shore profile with wl values along it

### IMPORTS ###
import os
import datetime
import time
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
from numpy import random
from scipy.interpolate import interp1d, interp2d, interpn

### FUNCTIONS ###
def coordSys_madbeach(E, N):
    '''
    Transform UTM to a local "madbeach" coordinate system

     Input:
       E, N (float) = UTM coordinates

     Output:
       X, Y (float) = local coordinates

     Usage: x, y= coordSys_madbeach(E,N)

     Written By:
       Jenna Brown, 2017, USGS
           based on Dave Thompson's madbeach_transformUTMToLocal.m
    '''

    #parameters from Madeira Beach calibration (intrinsics, extrinsics) - 01/28/2022
    #angle of rotation, degrees
    theta = 132.246888862048
    #origin Easting
    E0 = 323055.06675
    #origin Northing
    N0 = 3075922.21413775

    X,Y = xyRotate(E, N, theta, xo=E0, yo=N0)
    return X,Y


def xyRotate(x, y, theta, xo=0, yo=0):
    '''
    xyRotate  Rotate data.
    XR, YR = xyRotate(x, y, theta) rotates the coordinates x, y, and theta
    (cartesian decimal degrees).

    XR, YR = xyRotate(x, y, theta, xo, yo) rotates the coordinates around the
    origin xo, yo.

    Inputs:
        x (float) - x coordinate
        y (float) - y coordinate
        theta (float) - angle of rotation
        xo (float) - x origin
        yo (float) - y origin
    Outputs:
        XR (float) - rotated x coordinate
        YR (float) - rotated y coordinate
    '''
    #convert theta to radians
    theta = np.deg2rad(theta)
    
    #rotation matrix. Convert theta from degrees to radian for correct output
    A = np.array([[np.cos(theta), np.sin(theta)], [(-1)*np.sin(theta), np.cos(theta)]])
    
    outX = x - xo
    outY = y - yo
    out = np.array([outX, outY])
    #transpose from shape 2,len(247) to shape len(x),247. e.g. 2,247 --> 247,2
    out = np.transpose(out)
    out = np.matmul(out, A)

    XR = out[0:len(x), 0]
    YR = out[0:len(y), 1]
    return XR,YR


def findUVnDOF(betas, xyz, lcp):
    '''
    [UV] = findUVnDOF(beta, xyz, lcp)
    
    computes the distorted UV coordinates that correspond to a set of xyz
    points where the extrinsic parameters are specified by beta
    There are numerous options for beta, depending on what is known or not 
    known.  xyz is assumed to be an Nx3 list of world coordinates
    beta options
      length 6 (6 dof)    - [xCam yCam zCam azimuth tilt roll]
      length 5 (5 dof)    - [xCam yCam zCam azimuth tilt]
      length 3 (3 dof)    - [azimuth tilt roll]
      length 2 (2 dof)    - [azimuth tilt]
    
    When a beta is less than 6, the missing parameters are assumed to be
    passed as global variables in a structure lcp.lcp, .knownFlags and 
    .knowns where the knowns are listed in their natural order.  
    This structure is passed as an optional third argument in the
    function call normally.  However, if this routine is called as the forward 
    model in nlinfit, only two input arguments are allowed.  In that case,
    the globally-passed version, g2, is used to fill in these variables.
    If beta is a 1x6 vector, the global information is still needed for the
    lcp structure.
    
    NOTE - this returns DISTORTED COORDINATES.  THEY ARE ALSO RETURNED AS A
    COLUMN VECTOR of U(:); V(:) for use with nlinfit!

    Inputs:
        beta (np.array) - numpy array of extrinsic parameters (floats)
        xyz  (list) - list of numpy arrays of xyz coordinates to be converted to UV
        lcp (dict) - dictionary of values used for moving between real world and pixel coordinates
    Outputs:
        UV (np array) - UV coordinates
    '''

    P = lcpBeta2P(lcp, betas)

    onesColumn = np.ones(shape=(1, xyz.shape[1]))
    xyzStack = np.concatenate((xyz, onesColumn),axis=0)
    UV = np.matmul(P, xyzStack)

    #equivalent to repmat(UV(3,:),3,1) in Matlab
    repmat = np.tile(UV[2],(3,1))
    UV = np.divide(UV, repmat)

    U,V = DJIdistort(UV[0], UV[1], lcp)
    UV = np.concatenate((U, V))
    return UV

    
def lcpBeta2P(lcp, betas):
    '''
    Create a P matric from the lcp and beta. Used by findUV6DOF() and findXYZnDOF and to make
    things for pixel geometry toolbox.

    Inputs:
        lcp (dict) - dictionary of values used for moving between real world and pixel coordinates
        betas (np.array) - numpy array of extrinsic parameters (floats)
    Ouputs:
        P (np.array) - 3x3 P matrix
    '''

    #because of way lcp variables are stored in .mat files, have to use weird indexing to get number value. e.g. lcp['fx']['fx'][0][0]
    K = np.array([[lcp['fx']['fx'][0][0], 0, lcp['c0U']['c0U'][0][0]], [0, -1 * lcp['fy']['fy'][0][0], lcp['c0V']['c0V'][0][0]], [0, 0, 1]])

    R = angles2R(betas[3], betas[4], betas[5])
    #equivalent to -betas(1:3)' in matlab
    negativeTranspose = np.atleast_2d(-betas[0:3]).T
    IC = np.concatenate((np.eye(3), negativeTranspose), axis=1)
    KR = np.matmul(K, R)
    P = np.matmul(KR, IC)
    P = np.divide(P, P[2, 3])
    return P


def angles2R(a, t, r):
    '''
    Makes rotation matrix from input azimuth, tilt, and swing (roll)

    Inputs:
        a (float) - azimuth
        t (float) - tilt
        r (float) - roll
    Outputs:
        R (np.array) - 3x3 rotation matrix
    '''
    
    R = np.ndarray(shape=(3,3))
    R[0, 0] = np.cos(a) * np.cos(r) + np.sin(a) * np.cos(t) * np.sin(r)
    R[0, 1] = -np.cos(r) * np.sin(a) + np.sin(r) * np.cos(t) * np.cos(a)
    R[0, 2] = np.sin(r) * np.sin(t)
    R[1, 0] = -np.sin(r) * np.cos(a) + np.cos(r) * np.cos(t) * np.sin(a)
    R[1, 1] = np.sin(r) * np.sin(a) + np.cos(r) * np.cos(t) * np.cos(a)
    R[1, 2] = np.cos(r) * np.sin(t)
    R[2, 0] = np.sin(t) * np.sin(a)
    R[2, 1] = np.sin(t) * np.cos(a)
    R[2, 2] = -np.cos(t)
    return R


def DJIdistort(u, v, lcp):
    '''
    converts from undistorted to distorted pixel locations for a DJI phantom.
    This is based on equations from the Caltech lens distortion manuals.  
    lcp contains all the relevant intrinsic as well as radial and tangential
    distortion coefficients.

    where UV is a matrix of 2 columns [U, V].

    Inputs:
       u (np.array) - undistorted pxiel locations (rows)
       v (np.array) - undistorted pxiel locations (columns)
       lcp (dict) -  dictionary of values used for moving between real world and pixel coordinates
    Outputs:
        ud (np.array) - distorted pixel coordinates (rows)
        vd (np.array) - distorted pixel coordinates (columns)
    '''
    
    #find the range dependent correction factor, fr.
    
    x = (u - lcp['c0U']['c0U'][0][0]) / lcp['fx']['fx'][0][0]
    y = (v - lcp['c0V']['c0V'][0][0]) / lcp['fy']['fy'][0][0]
    r2 = np.multiply(x, x) + np.multiply(y, y)
    
    #need to convert np.object lcp to dict with tolist(), then unpack values r and fr, then squeeze into 1 dimension
    lcp_r = lcp['r'].tolist()
    lcp_r = lcp_r['r']
    lcp_r = lcp_r.squeeze()
    lcp_fr = lcp['fr'].tolist()
    lcp_fr = lcp_fr['fr']
    lcp_fr = lcp_fr.squeeze()
    #can interp1d creates function, which can be called in the same line
    fr = interp1d(lcp_r, lcp_fr, fill_value='extrapolate')(np.sqrt(r2))

    #now do 2d interpolation for dx, dy
    lcp_x = lcp['x'].tolist()
    lcp_x = lcp_x['x']
    lcp_y = lcp['y'].tolist()
    lcp_y = lcp_y['y']
    lcp_dx = lcp['dx'].tolist()
    lcp_dx = lcp_dx['dx']
    lcp_dx = lcp_dx.squeeze()
    lcp_dy = lcp['dy'].tolist()
    lcp_dy = lcp_dy['dy']
    lcp_dy = lcp_dy.squeeze()

    #creates square grid but only want dx and dy to be 1 dimensional, so take last element
    dx = interp2d(lcp_x, lcp_y, lcp_dx)(x, y)[0]
    dy = interp2d(lcp_x, lcp_y, lcp_dy)(x,y)[0]

    x2 = np.multiply(x, fr) + dx
    y2 = np.multiply(y, fr) + dy
    ud = (x2 * lcp['fx']['fx'][0][0]) + lcp['c0U']['c0U'][0][0]      #answer in chip pixel units
    vd = (y2 * lcp['fy']['fy'][0][0]) + lcp['c0V']['c0V'][0][0]
    return ud, vd


def distortUV(U, V, lcp):
    '''
    This function distorts undistorted UV coordinates using distortion
    models from from the Caltech lens distortion manuals.The function also
    suggests whether the UVd coordinate is valid (not having tangential
    distortion values bigger than what is at the corners and being within
    the image).

    Inputs:
       U (np.array) - undistorted pxiel locations (rows)
       V (np.array) - undistorted pxiel locations (columns)
       lcp (dict) -  dictionary of values used for moving between real world and pixel coordinates
    Outputs:
        Ud (np.array) - distorted pixel coordinates (rows)
        Vd (np.array) - distorted pixel coordinates (columns)
    '''

    NU = lcp['NU']['NU'][0][0]
    NV = lcp['NV']['NV'][0][0]
    c0U = lcp['c0U']['c0U'][0][0]
    c0V = lcp['c0V']['c0V'][0][0]
    fx = lcp['fx']['fx'][0][0]
    fy = lcp['fy']['fy'][0][0]
    d1 = lcp['d1']['d1'][0][0]
    d2 = lcp['d2']['d2'][0][0]
    d3 = lcp['d3']['d3'][0][0]
    t1 = lcp['t1']['t1'][0][0]
    t2 = lcp['t2']['t2'][0][0]

    #normalize distances
    x = (U - c0U) / fx
    y = (V - c0V) / fy

    #radial distortion
    r2 = np.multiply(x, x) + np.multiply(y, y)
    fr = 1 + (d1 * r2) + (d2 * np.multiply(r2, r2)) + d3 * (np.multiply(r2, np.multiply(r2, r2)))

    #tangential distortion
    dx = (2 * t1 * np.multiply(x,y)) + t2 * (r2 + (2 * np.multiply(x, x)))
    dy = t1 * (r2 + (2 * np.multiply(y, y))) + (2 * t2 * np.multiply(x, y))

    #apply Correction, answer in chip pixel units
    xd = np.multiply(x, fr) + dx
    yd = np.multiply(y, fr) + dy
    Ud = (xd * fx) + c0U
    Vd = (yd * fy) + c0V

    return Ud, Vd
    
    
    
### MAIN ###
#load image
snapFile = '1665860400.Sat.Oct.15_19_00_00.GMT.2022.madbeach.c1.snap.jpg'
timexFile = '1665860400.Sat.Oct.15_19_00_00.GMT.2022.madbeach.c1.timex.jpg'
snap = plt.imread(snapFile)
timex = plt.imread(timexFile)

#load data
geom_file = 'geomFile_c1.mat'
geom_c1 = scipy.io.loadmat('./matlabcode/' + geom_file)
lcp = {}
lcp['c0U'] = scipy.io.loadmat('./matlabcode/globals/c0U.mat')
lcp['c0V'] = scipy.io.loadmat('./matlabcode/globals/c0V.mat')
lcp['d1'] = scipy.io.loadmat('./matlabcode/globals/d1.mat')
lcp['d2'] = scipy.io.loadmat('./matlabcode/globals/d2.mat')
lcp['d3'] = scipy.io.loadmat('./matlabcode/globals/d3.mat')
lcp['dx'] = np.array(scipy.io.loadmat('./matlabcode/globals/dx.mat'))
lcp['dy'] = np.array(scipy.io.loadmat('./matlabcode/globals/dy.mat'))
lcp['fr'] = np.array(scipy.io.loadmat('./matlabcode/globals/fr.mat'))
lcp['fx'] = scipy.io.loadmat('./matlabcode/globals/fx.mat')
lcp['fy'] = scipy.io.loadmat('./matlabcode/globals/fy.mat')
lcp['NU'] = scipy.io.loadmat('./matlabcode/globals/NU.mat')
lcp['NV'] = scipy.io.loadmat('./matlabcode/globals/NV.mat')
lcp['r']= np.array(scipy.io.loadmat('./matlabcode/globals/r.mat'))
lcp['t1'] = scipy.io.loadmat('./matlabcode/globals/t1.mat')
lcp['t2'] = scipy.io.loadmat('./matlabcode/globals/t2.mat')
lcp['x'] = np.array(scipy.io.loadmat('./matlabcode/globals/x.mat'))
lcp['y'] = np.array(scipy.io.loadmat('./matlabcode/globals/y.mat'))

#use squeeze() to get rid of unnecessary dimensions
betas = geom_c1['betas'].squeeze()

#if 2022 geom file, offset azimuth by rotation angle
if geom_file == 'geomFile_c1.mat':
    betas[3] = betas[3] - np.deg2rad(132.246888862048);

#get datetime elements from filename
filenameElements = snapFile.split('.')
epoch = filenameElements[0]
dayOfWeek = filenameElements[1]
shortMonth = filenameElements[2]
month = time.strptime(shortMonth, '%b').tm_mon
dayHourSecond = filenameElements[3]
day = dayHourSecond.split('_')[0]
hour = dayHourSecond.split('_')[1]
minute = dayHourSecond.split('_')[2]
second = dayHourSecond.split('_')[3]
year = filenameElements[5]
imgDatetimeStr = f'{year}-{month}-{day} {hour}:{minute}:{second}'
imgDatetime = datetime.datetime.strptime(imgDatetimeStr, "%Y-%m-%d %H:%M:%S")

###display image
##plt.imshow(snap)
##plt.show()
##plt.imshow(timex)
##plt.show()

surveysAll = os.listdir('//gs/stpetersburgfl-g/NACCH/Imagery/madbeach/surveys/walking/')
i = 0
for elem in surveysAll:
    #want only most recent folder
    if elem.startswith('202'):
        if i == 0:
            mostRecent = elem
            i = i + 1
        else:
            if int(elem) > int(mostRecent):
                #2nd most recent in case most recent survey hasn't been processed yet
                secondMostRecent = mostRecent
                mostRecent = elem
                
mostRecentDir = os.listdir('//gs/stpetersburgfl-g/NACCH/Imagery/madbeach/surveys/walking/' + mostRecent)
#check that most recent survey folder has .xyz files
hasXYZ = False
for file in mostRecentDir:
    if file.endswith('.xyz'):
        hasXYZ = True
        break

#yes, most recent folder has .xyz files
if hasXYZ:
    mostRecentLine60 = '//gs/stpetersburgfl-g/NACCH/Imagery/madbeach/surveys/walking/' + mostRecent + '/line60.xyz'
else:
    mostRecentLine60 = '//gs/stpetersburgfl-g/NACCH/Imagery/madbeach/surveys/walking/' + secondMostRecent + '/line60.xyz'

###teting 2017 data
##mostRecentLine60 = '//gs/stpetersburgfl-g/NACCH/Imagery/madbeach/surveys/walking/20170217/line59.xyz'

#process line60
profileE = []
profileN = []
profileZ = []
with open(mostRecentLine60, 'r') as xyzFile:
    for line in xyzFile:
        xyz = line.split(' ')
        x = float(xyz[0])
        y = float(xyz[1])
        z = xyz[2]
        #z has newiine character at the end
        z = float(z.replace('\n', ''))
        profileE.append(x)
        profileN.append(y)
        profileZ.append(z)

profileE = np.array(profileE)
profileN = np.array(profileN)
profileZ = np.array(profileZ)
profileX, profileY = coordSys_madbeach(profileE, profileN)

### estimating line 60 position to be -87.7 based on surveys that measured line 60
##profileY = np.ones(shape=profileY.shape) * (-87.7)

#change NaN to 0s
for i, elem in enumerate(profileZ):
    if (elem == np.nan) or (elem == None):
        profileZ[i] = 0

uniqueZ = np.unique(profileZ)
a, b = np.histogram(profileZ, bins=uniqueZ)
mult = np.argwhere(a > 1)
#catch >2 identical values in profileZ
for i in range(0, len(mult)):
    #b is ascending array of unique values. mult is array of bins (arrays) in a with dupllicate values. indices in a match indices in b. So, index mult[i][0] will give index of a value in b that is not unique
    dupValue = b[mult[i][0]]
    dupIndices = np.argwhere(profileZ == dupValue)
    
    for ii in range(0, len(dupIndices)):
        #dupIndices[ii][0] is index in profileZ that contains duplicate values
        index = dupIndices[ii][0]
        #slightly alter value to prevent duplicates
        offset = random.rand() * 0.0000000001
        profileZ[index] = profileZ[index] + offset

uniqueZ = np.unique(profileZ)
a, b = np.histogram(profileZ, bins=uniqueZ)
mult = np.argwhere(a > 1)

lat = 27.796216949206798
lon = -82.796102542950635

all_twl = scipy.io.loadmat('//gs/StPetersburgFL-G/NACCH/Imagery/madbeach/runup/all_TWL_forecast.mat')
Rtime = all_twl['all_twl']['forecastTime']
Rrunup05 = all_twl['all_twl']['runup05']
Rrunup = all_twl['all_twl']['runup']
Rrunup95 = all_twl['all_twl']['runup95']
Rtwl05 = all_twl['all_twl']['twl05']
Rtwl = all_twl['all_twl']['twl']
Rtwl95 = all_twl['all_twl']['twl95']
Rslope = all_twl['all_twl']['slope']

#only save madbeach data
pindex = 0
#find twl data that matches image date. Because of scipy loadmat(), array is 3-layer deep array
for i in range(0, len(Rtime)):
    #convert string to datetime object
    RDatetime = datetime.datetime.strptime(Rtime[i][0][0], "%Y-%m-%d %H:%M:%S.0")

    if (RDatetime < imgDatetime + datetime.timedelta(hours=1)) and (RDatetime > imgDatetime - datetime.timedelta(hours=1)):
        tindex = i

#loop through alongshore locations
num = [x for x in range(-25, -451, -1)]
for i in range(0, len(num)):
    m = Rslope[tindex, pindex]

    R2z = []
    R2z.append(Rrunup05[tindex, pindex][0][0])
    R2z.append(Rrunup[tindex, pindex][0][0])
    R2z.append(Rrunup95[tindex, pindex][0][0])
    R2z = np.array(R2z)

    #interp1d returns a 1d interpolation function (interpolate)
    interpolate = interp1d(profileZ, profileX)

    R2x = []
    interpRX = interpolate(R2z)
    #interpX is array with shape (3, 1, 1). Need to extract individual values and append to R2x
    R2x.append(interpRX[0])
    R2x.append(interpRX[1])
    R2x.append(interpRX[2])
    R2x = np.array(R2x)
    R2y = num[i]*np.ones(len(R2x))

    TWLz = []
    TWLz.append(Rtwl05[tindex,pindex][0][0])
    TWLz.append(Rtwl[tindex,pindex][0][0])
    TWLz.append(Rtwl05[tindex,pindex][0][0])

    TWLx = []
    interpTWLX = interpolate(TWLz)
    TWLx.append(interpTWLX[0])
    TWLx.append(interpTWLX[1])
    TWLx.append(interpTWLX[2])
    TWLy = num[i]*np.ones(len(TWLx))

    ###***world image coordinates***

    #-- profile --
    xyz = np.array([profileX, profileY, profileZ])
    UV = findUVnDOF(betas, xyz, lcp)
    UV = np.round(UV)
    UV = np.reshape(UV, (2, int(len(UV)/2)))
    UV = UV.T
    profileU = UV[:,0]
    profileV = UV[:,1]

    #-- R2 --
    if i == 0:
        R2u = np.ndarray(shape=(3,len(num)))
        R2v = np.ndarray(shape=(3,len(num)))
    xyz = np.array([R2x, R2y, R2z])
    UV = findUVnDOF(betas, xyz, lcp)
    UV = np.round(UV)
    UV = np.reshape(UV, (2, int(len(UV)/2)))
    UV = UV.T
    R2u[:,i] = UV[:,0]
    R2v[:,i] = UV[:,1]

    #-- TWL --
    if i == 0:
        TWLu = np.ndarray(shape=(3,len(num)))
        TWLv = np.ndarray(shape=(3,len(num)))
    xyz = np.array([TWLx, TWLy, TWLz])
    UV = findUVnDOF(betas, xyz, lcp)
    UV = np.round(UV)
    UV = np.reshape(UV, (2, int(len(UV)/2)))
    UV = UV.T
    TWLu[:,i] = UV[:,0]
    TWLv[:,i] = UV[:,1]

## plot forecast ##
snapImg = plt.imread(snapFile)
timexFile = plt.imread(timexFile)
fig, ax = plt.subplots()
ax.imshow(snapImg)
ax.plot(TWLu, TWLv, '-', color='red')
plt.show()

    

