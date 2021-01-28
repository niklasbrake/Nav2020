import os
import pyabf
import numpy as np
import mat4py as m4p
import sys

# Get directory containing folders for each experiment
folder = sys.argv[1]
# This is the list of experiment folders
SubSubFolders = np.array(os.listdir(folder))
# Turn it into a full file path
# SubFolders = np.core.defchararray.add(DataDirectory+'\\',SubFolders) 

# # Loop over each experiment folder
for file in SubSubFolders:
    # If the file is an abf file..
    ext = os.path.splitext(file)[-1].lower() 
    if ext == '.abf':
        filename = os.path.splitext(file)[0].lower() 
        savedName = folder + '\\' + filename.split('_')[-1] + '.mat'
        # And the folder doesn't already contain the respective .mat file
        if os.path.isfile(savedName):
            continue
        # Load the abf data and save the voltage, current, and epoch data
        ABFData = pyabf.ABF(folder + '\\' + file)
        V = np.zeros([ABFData.sweepCount,len(ABFData.sweepY)])
        I = np.zeros([ABFData.sweepCount,len(ABFData.sweepY)])
        Epochs = ABFData.sweepEpochs.p1s;
        for i in ABFData.sweepList:
            ABFData.setSweep(i)
            V[i,:] = ABFData.sweepC
            I[i,:] = ABFData.sweepY
        V = V + ABFData.data[1,:ABFData.sweepEpochs.p2s[0]].mean()
        # Data = ABFData.data
        data = {'Voltage':V.tolist(),'Current':I.tolist(),'Epochs':Epochs}
        m4p.savemat(savedName, data)