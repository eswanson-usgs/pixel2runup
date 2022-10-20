##%% Plot TWL estimates on raw camera image
##% also plots measured runup values on camera image
##% plots cross-shore profile with wl values along it

### IMPORTS ###
import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt


### FUNCTIONS ###


### MAIN ###
#load image data
snapFile = '1665860400.Sat.Oct.15_19_00_00.GMT.2022.madbeach.c1.snap.jpg'
timexFile = '1665860400.Sat.Oct.15_19_00_00.GMT.2022.madbeach.c1.timex.jpg'
snap = plt.imread(snapFile)
timex = plt.imread(timexFile)
geom = scipy.io.loadmat('./matlabcode/geomFile_c1.mat')

#get datetime elements from filename
filenameElements = snapFile.split('.')
epoch = filenameElements[0]
dayOfWeek = filenameElements[1]
shortMonth = filenameElements[2]
dayHourSecond = filenameElements[3]
day = dayHourSecond.split('_')[0]
hour = dayHourSecond.split('_')[1]
second = dayHourSecond.split('_')[2]
year = filenameElements[5]

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
                mostRecent = elem


