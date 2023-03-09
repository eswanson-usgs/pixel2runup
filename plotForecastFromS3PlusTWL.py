##%% Plot TWL estimates on raw camera image
##% also plots measured runup values on camera image
##% plots cross-shore profile with wl values along it
###Take data from S3 and TWL viewere API

### IMPORTS ###
import os
import datetime
import time
import fsspec
import requests
import json
import csv
import scipy.io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from numpy import random
from tqdm import tqdm
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
    K = np.array([[lcp['fx'], 0, lcp['c0U']], [0, -1 * lcp['fy'], lcp['c0V']], [0, 0, 1]])

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
    
    x = (u - lcp['c0U']) / lcp['fx']
    y = (v - lcp['c0V']) / lcp['fy']
    r2 = np.multiply(x, x) + np.multiply(y, y)
    
    #need to convert np.object lcp to dict with tolist(), then unpack values r and fr, then squeeze into 1 dimension
    lcp_r = lcp['r'].tolist()
    lcp_fr = lcp['fr'].tolist()
    #can interp1d creates function, which can be called in the same line
    fr = interp1d(lcp_r, lcp_fr, fill_value='extrapolate')(np.sqrt(r2))

    #now do 2d interpolation for dx, dy
    lcp_x = lcp['x'].tolist()
    lcp_y = lcp['y'].tolist()
    lcp_dx = lcp['dx'].tolist()
    lcp_dy = lcp['dy'].tolist()

    #creates square grid but only want dx and dy to be 1 dimensional, so take last element
    dx = interp2d(lcp_x, lcp_y, lcp_dx)(x, y)[0]
    dy = interp2d(lcp_x, lcp_y, lcp_dy)(x,y)[0]

    x2 = np.multiply(x, fr) + dx
    y2 = np.multiply(y, fr) + dy
    ud = (x2 * lcp['fx']) + lcp['c0U']      #answer in chip pixel units
    vd = (y2 * lcp['fy']) + lcp['c0V']
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

    NU = lcp['NU']
    NV = lcp['NV']
    c0U = lcp['c0U']
    c0V = lcp['c0V']
    fx = lcp['fx']
    fy = lcp['fy']
    d1 = lcp['d1']
    d2 = lcp['d2']
    d3 = lcp['d3']
    t1 = lcp['t1']
    t2 = lcp['t2']

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


def add_impact(row):
    if row == 0:
        return 'None'
    elif row == 1:
        return 'Collision'
    elif row == 2:
        return 'Overwash'
    elif row == 3:
        return 'Inundation'


def twl_cc_api_download(region_id, init_date, site_id, final_date=None, refresh_id=False,
                        refresh_df=False, pause=0.2, msl_2_navd=-0.096, external=True):
    """
    Download USGS TWL&CC data using the API and store everything in a DataFrame

    region_id: Int with the ID of the region to download from
    init_date: String with the first day of data to download (yyyy-mm-dd)
    site_id: Int with the TWL&CC site ID number
    final_date: String with the last day of data to download (yyyy-mm-dd)
    refresh_id: Bool to force a re-grab of the forecast IDs
    refresh_df: Bool to force a re-grab of the forecast
    pause: Float with the amount of seconds between API calls (Default = 0.2)
    msl_2_navd: Float to convert MSL elevation to NAVD88 (Set to 0 to leave as MSL)
    external: Bool to use the external (True) or internal(?, False) server. Use for testing DDOS failures

    Outputs:
        twl_cc_df (pandas dataframe) - dataframe that holds the resulting forecast(s) from the TWL & CC API
        usedPreviousDay (bool) - flag that relays whether or not a forecast from the previous day was used. This
                                 flag is True when no forecasts are found from the initial query. This is really used
                                 for Eric's querying of a single day (today)

    Written by Michael Itzkin
    Edits by Eric Swanson
    """

    # Set the final date to tomorrow if None
    if final_date is None:
        final_date = datetime.datetime.today() + datetime.timedelta(days=1)
        final_date = final_date.strftime(format='%Y-%m-%d')

    # Build the base request string
    if external:
        API_STR = f'https://coastal.er.usgs.gov/hurricanes/research/twlviewer/api/regions/{region_id}/forecasts'
    else:
        API_STR = f'http://coastal-dmz.er.usgs.gov/hurricanes/research/twlviewer/api/regions/{region_id}/forecasts'

    # Get the total number pages of data
    r = requests.get(f'{API_STR}?date=gte_{init_date}|lt_{final_date}&pageSize=100')
    total_forecast_pages = json.loads(r.headers['X-Pagination'])['totalPages']
    print(f'There are {total_forecast_pages} pages of forecasts to download...')

    #if no forecasts available for init date, use a forecast from the day before
    usedPreviousDay = False
    if total_forecast_pages == 0:

        print(f'No forecasts available for {init_date}-using forecast from previous day')
        usedPreviousDay = True
        previousDay = datetime.datetime.strptime(init_date,'%Y-%m-%d') - datetime.timedelta(days=1)
        init_date = previousDay.strftime(format='%Y-%m-%d')
        r = requests.get(f'{API_STR}?date=gte_{init_date}|lt_{final_date}&pageSize=100')
        total_forecast_pages = json.loads(r.headers['X-Pagination'])['totalPages']
        print(f'There are {total_forecast_pages} pages of forecasts to download...')
        

    # Loop over the pages of data and store the forecast IDs
    forecast_id_file = 'TWLCC Forecast IDs.npy'
    if not os.path.isfile(forecast_id_file) or refresh_id:
        forecast_ids = []
        print('Getting forecast ID numbers:')
        for page in tqdm(range(1, total_forecast_pages + 1)):
            print(page)
            r = requests.get(f'{API_STR}?date=gte_{init_date}|lt_{final_date}&pageSize=100&pageNumber={page}')
            page_data = r.json()
            for item in page_data:
                forecast_ids.append(item['id'])
        np.save(forecast_id_file, np.array(forecast_ids))
    else:
        print('Loading previously collected forecast IDs')
        forecast_ids = np.load(forecast_id_file)
    print('Finished getting forecast ID numbers. Getting TWL&CC forecast records:')

    # Loop over the forecast IDs
    forecast_df_file = 'TWLCC Forecast.csv'
    if not os.path.isfile(forecast_df_file) or refresh_df:
        twl_cc_df_list = []
        for _id in tqdm(forecast_ids):
            try:
                twl_cc_data = requests.get(f'{API_STR}/{_id}/sites/{site_id}/waterLevels')
                json_data = twl_cc_data.json()

                # The ID number changed at some point so if the result is blank try the older ID number
                if len(json_data) == 0:
                    twl_cc_data = requests.get(f'{API_STR}/{_id}/sites/{3195}/waterLevels')
                    json_data = twl_cc_data.json()

                # The ID number changed at some point again so if the result is blank try the older ID number
                if len(json_data) == 0:
                    twl_cc_data = requests.get(f'{API_STR}/{_id}/sites/{307}/waterLevels')
                    json_data = twl_cc_data.json()

                temp_df = pd.DataFrame.from_records(json_data)
                twl_cc_df_list.append(temp_df)
                time.sleep(pause)
            except:
                pass
        twl_cc_df = pd.concat(twl_cc_df_list).reset_index(drop=True)
        twl_cc_df.to_csv(forecast_df_file, index=False)
    else:
        print('Loading previously downloaded forecasts')
        twl_cc_df = pd.read_csv(forecast_df_file, header=0, delimiter=',')
    print('Finished downloading TWL&CC records!')

    # Force data types
    use_cols = twl_cc_df.columns[twl_cc_df.columns != 'dateTime']
    twl_cc_df[use_cols] = twl_cc_df[use_cols].apply(pd.to_numeric, errors='coerce')
    twl_cc_df['Date'] = pd.to_datetime(twl_cc_df['dateTime'], format='%Y-%m-%d %H:%M:%S')

    # Convert MSL to NAVD88
    msl_columns = ['twl', 'twl05', 'twl95', 'setup', 'runup', 'runup05', 'runup95',
                   'tideWindSetup', 'swash', 'incSwash', 'infragSwash', 'hs']
    twl_cc_df[msl_columns] += msl_2_navd

    # Add additional columns
    twl_cc_df['L'] = (9.81 * (twl_cc_df['pp']**2)) / (2 * np.pi)
    twl_cc_df['predictedImpactCode'] = twl_cc_df['predictedImpactCode'].apply(add_impact)

    # Sort the data by date and remove repeated forecasts.
    twl_cc_df = twl_cc_df.sort_values(by='Date').drop_duplicates(subset=['Date'], keep='first').reset_index(drop=True)

    return twl_cc_df, usedPreviousDay

      
    
### MAIN ###
print('start:',datetime.datetime.now())

file_system = fsspec.filesystem('s3', profile='coastcam')

#load image
snapFile = '1665860400.Sat.Oct.15_19_00_00.GMT.2022.madbeach.c1.snap.jpg'
snap = plt.imread(snapFile)

#load data
geom_file_s3 = 's3://cmgp-coastcam/cameras/madeira_beach/forecast/geomFile_c1.mat'
geom_file_local = 'geomFile_c1.mat'
file_system.download(geom_file_s3, geom_file_local)
geom_c1 = scipy.io.loadmat(geom_file_local, squeeze_me=True, struct_as_record=False)

lcp = {}
lcp['c0U'] = geom_c1['meta'].globals.lcp.c0U
lcp['c0V'] = geom_c1['meta'].globals.lcp.c0V
lcp['d1'] = geom_c1['meta'].globals.lcp.d1
lcp['d2'] = geom_c1['meta'].globals.lcp.d2
lcp['d3'] = geom_c1['meta'].globals.lcp.d3
lcp['dx'] = np.array(geom_c1['meta'].globals.lcp.dx)
lcp['dy'] = np.array(geom_c1['meta'].globals.lcp.dy)
lcp['fr'] = np.array(geom_c1['meta'].globals.lcp.fr)
lcp['fx'] = geom_c1['meta'].globals.lcp.fx
lcp['fy'] = geom_c1['meta'].globals.lcp.fy
lcp['NU'] = geom_c1['meta'].globals.lcp.NU
lcp['NV'] = geom_c1['meta'].globals.lcp.NV
lcp['r']= np.array(geom_c1['meta'].globals.lcp.r)
lcp['t1'] = geom_c1['meta'].globals.lcp.t1
lcp['t2'] = geom_c1['meta'].globals.lcp.t2
lcp['x'] = np.array(geom_c1['meta'].globals.lcp.x)
lcp['y'] = np.array(geom_c1['meta'].globals.lcp.y)

#use squeeze() to get rid of unnecessary dimensions
betas = geom_c1['betas'].squeeze()

#if 2022 geom file, offset azimuth by rotation angle
if geom_file_local == 'geomFile_c1.mat':
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
#imgDatetime = datetime.datetime.strptime(imgDatetimeStr, "%Y-%m-%d %H:%M:%S")
######JUST FOR TESTING
imgDatetime = datetime.datetime.now()

surveysAll = file_system.glob('s3://cmgp-coastcam/cameras/madeira_beach/forecast/')
i = 0
for elem in surveysAll:
    
    filename = elem.split('/')[-1]
    filedate = filename[0:8]
    
    #want only most recent survey
    if filename.startswith('202'):
        if i == 0:
            mostRecentDate = filedate
            mostRecentFile = elem
            mostRecentLine60 = filename
            i = i + 1
        else:
            if int(filedate) > int(mostRecentDate):
                mostRecentFile = elem
                mostRecentLine60 = filename
                
file_system.download('s3://' + mostRecentFile, mostRecentLine60)

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

#access TWL data via API
region_id = 4
init_date = datetime.datetime.today().strftime(format='%Y-%m-%d')
site_id = 3667
twl_cc_df, usedPreviousDay = twl_cc_api_download(region_id, init_date, site_id, refresh_id=True, refresh_df=True)

#####code for testing opening csv and reading data without using Pandas dataframe
##with open('TWLCC Forecast.csv', newline='') as csvfile:
##    reader = csv.DictReader(csvfile, delimiter=',')
##    Rtime = []
##    Rrunup05 = []
##    Rrunup = []
##    Rrunup95 = []
##    Rtwl05 = []
##    Rtwl = []
##    Rtwl95 = []
##    for row in reader:
##        Rtime.append(row['dateTime'])
##        Rrunup05.append(float(row['runup05']))
##        Rrunup.append(float(row['runup']))
##        Rrunup95.append(float(row['runup95']))
##        Rtwl05.append(float(row['twl05']))
##        Rtwl.append(float(row['twl']))
##        Rtwl95.append(float(row['twl95']))

Rtime = twl_cc_df['dateTime']
Rrunup05 = twl_cc_df['runup05']
Rrunup = twl_cc_df['runup']
Rrunup95 = twl_cc_df['runup95']
Rtwl05 = twl_cc_df['twl05']
Rtwl = twl_cc_df['twl']
Rtwl95 = twl_cc_df['twl95']

#find twl data that matches image date. Because of scipy loadmat(), array is 3-layer deep array
for i in range(0, len(Rtime)):
    #convert string to datetime object
    RDatetime = datetime.datetime.strptime(Rtime[i], "%Y-%m-%d %H:%M:%S")

    #If using forecast from previous day, look for hour closest (or same as the image's)
    if not usedPreviousDay:
        dateDiff = abs(imgDatetime - RDatetime)
    else:
        previousDay = imgDatetime - datetime.timedelta(days=1)
        dateDiff = abs(previousDay - RDatetime)

    #compare times from dataframe to image time. 
    if i == 0:
        lowestDiff = dateDiff
        forecastTime = Rtime[i]
        tindex = i
    else:
        if dateDiff < lowestDiff:
            lowestDiff = dateDiff
            forecastTime  = Rtime[i]
            tindex = i
    
#loop through alongshore locations
num = [x for x in range(-25, -451, -1)]
for i in range(0, len(num)):

    R2z = []
    R2z.append(Rrunup05[tindex])
    R2z.append(Rrunup[tindex])
    R2z.append(Rrunup95[tindex])
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
    TWLz.append(Rtwl05[tindex])
    TWLz.append(Rtwl[tindex])
    TWLz.append(Rtwl05[tindex])

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
fig, ax = plt.subplots()
ax.imshow(snap)
ax.plot(TWLu, TWLv, '-', color='red')
plt.show()

print('End:', datetime.datetime.now())

    

