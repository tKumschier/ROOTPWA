#!/usr/bin/env python
# coding: utf-8

# In[18]:


import uproot
import numpy as np
import os
import math
import matplotlib.pyplot as plt
#import progressbar
import multiprocessing as mp
from tqdm import * # tqdm_notebook as tqdm
import argparse
import time


# In[2]:

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", help="The name in which the data should stored; Format: min_max_BinWidth; Example: 500_600_1")
args = parser.parse_args()

# In[3]:


mainPath = "/nfs/freenas/tuph/e18/project/compass/analysis/tkumschier/ROOTPWA/test/TEST_FOLDER_MASS/"
testFoldersTmp = []
for filename in sorted(os.listdir(mainPath)):
    testFoldersTmp.append(filename)
testFolders = (np.sort(np.array(testFoldersTmp).astype(int))).astype(str)


# # Berechnung der Integralmatrix 
# ## Definiere Funktionen

# In[4]:


def prepareArray(origArray):
    preparedArray = np.zeros(len(origArray), dtype=complex)
    for i in range(0, len(origArray)):
        preparedArray[i] = origArray[i][0]
    return preparedArray


# In[36]:


## Lade Daten in Matrix
def calcErrors(folderI):
    #if bar.currval <= folderI:
    #    bar.update(folderI)
    ampPath = mainPath + testFolders[folderI] + '/amps/'
    onlyfiles = []
    for filename in sorted(os.listdir(ampPath)):
        onlyfiles.append(filename)
    nmbWaves = len(onlyfiles)
    ampFile1 = uproot.open(ampPath + onlyfiles[0])
    nmbEvents = len(ampFile1[ampFile1.keys()[0]]['amplitude'].get('_incohSubAmps._real').array())
    integralMatrix = np.zeros((nmbWaves, nmbWaves), dtype = np.ndarray)
    sd = np.zeros((nmbWaves, nmbWaves), dtype = complex)
    ampMatrix = np.zeros((nmbWaves,nmbEvents), dtype = complex)

    integralMatrix = np.zeros((nmbWaves, nmbWaves), dtype = complex)
    errorMatrixAbsolut = np.zeros((nmbWaves, nmbWaves), dtype = complex)
    errorMatrixRelativ = np.zeros((nmbWaves, nmbWaves), dtype = complex)

    for waveIndexI in range(0, nmbWaves):
        if ampMatrix[waveIndexI].all() == 0:
            ampFile1 = uproot.open(ampPath + onlyfiles[waveIndexI])
            key1 = ampFile1.keys()[0]

            incohSubAmpsReal1 = ampFile1[key1]['amplitude'].get('_incohSubAmps._real').array()
            incohSubAmpsImag1 = ampFile1[key1]['amplitude'].get('_incohSubAmps._imag').array()
            ampMatrix[waveIndexI] = prepareArray(np.array(incohSubAmpsReal1 + 1j * incohSubAmpsImag1))
            #print ampMatrix[waveIndexI]
            #raise ValueError('Exit.')




        for waveIndexJ in range(0, waveIndexI + 1): #nmbWaves):
            if ampMatrix[waveIndexJ].all() == 0:
                ampFile2 = uproot.open(ampPath + onlyfiles[waveIndexJ])
                key2 = ampFile2.keys()[0]

                incohSubAmpsReal2 = ampFile2[key2]['amplitude'].get('_incohSubAmps._real').array()
                incohSubAmpsImag2 = ampFile2[key2]['amplitude'].get('_incohSubAmps._imag').array()
                ampMatrix[waveIndexJ] = prepareArray(np.array(incohSubAmpsReal2 + 1j * incohSubAmpsImag2))

            integral = np.sum(ampMatrix[waveIndexI] * ampMatrix[waveIndexJ].conjugate())/nmbEvents 

            sdReal = sum(pow((ampMatrix[waveIndexI] * ampMatrix[waveIndexJ].conjugate()).real, 2)).real
            sdImag = sum(pow((ampMatrix[waveIndexI] * ampMatrix[waveIndexJ].conjugate()).imag, 2)).real
            #print "Var1:", sdReal
            
            
            sdReal = math.sqrt(sdReal/float(nmbEvents) - pow(integral.real,2))
            sdImag = math.sqrt(sdImag/float(nmbEvents) - pow(integral.imag,2))
            sdComplex = complex(sdReal, sdImag)        
            errorComplex = complex(sdReal/math.sqrt(nmbEvents), sdImag/math.sqrt(nmbEvents))

            sd[waveIndexI][waveIndexJ] = sdComplex
            sd[waveIndexJ][waveIndexI] = sdComplex
   
            #print "%i, %i, %i: %s, %s, %s" %(folderI, waveIndexI, waveIndexJ, integral, sd[waveIndexI][waveIndexJ], errorComplex)
            #print np.cov(amp/nmbEvents)
            #raise ValueError('Exit.')


            integralMatrix[waveIndexI][waveIndexJ] = integral
            integralMatrix[waveIndexJ][waveIndexI] = integral
            errorMatrixAbsolut[waveIndexI][waveIndexJ] = errorComplex
            errorMatrixAbsolut[waveIndexJ][waveIndexI] = errorComplex
            if waveIndexI == waveIndexJ:
                errorMatrixRelativ[waveIndexI][waveIndexJ] = complex(errorComplex.real/math.fabs(integral.real), 0)
            else:
                errorMatrixRelativ[waveIndexI][waveIndexJ] = complex(errorComplex.real/math.fabs(integral.real), errorComplex.imag/math.fabs(integral.imag))
                errorMatrixRelativ[waveIndexJ][waveIndexI] = complex(errorComplex.real/math.fabs(integral.real), errorComplex.imag/math.fabs(integral.imag))

    return (folderI, integralMatrix, errorMatrixAbsolut, errorMatrixRelativ, testFolders[folderI])


# In[34]:


## Lade Daten in eine Variable, weniger Speicher notwendig
def calcErrors2(folderI):
    #if bar.currval <= folderI:
    #bar.update(folderI)
    ampPath = mainPath + testFolders[folderI] + '/amps/'
    onlyfiles = []
    for filename in sorted(os.listdir(ampPath)):
        onlyfiles.append(filename)
    nmbWaves = len(onlyfiles)
    ampFile1 = uproot.open(ampPath + onlyfiles[0])
    nmbEvents = len(ampFile1[ampFile1.keys()[0]]['amplitude'].get('_incohSubAmps._real').array())
    sd = np.zeros((nmbWaves, nmbWaves), dtype = complex)
    
    integralMatrix = np.zeros((nmbWaves, nmbWaves), dtype = complex)
    errorMatrixAbsolut = np.zeros((nmbWaves, nmbWaves), dtype = complex)
    errorMatrixRelativ = np.zeros((nmbWaves, nmbWaves), dtype = complex)
    
    '''
    integralMatrix = np.zeros((nmbWaves, nmbWaves), dtype = np.ndarray)
    errorMatrixAbsolut = np.zeros((nmbWaves, nmbWaves), dtype = np.ndarray)
    errorMatrixRelativ = np.zeros((nmbWaves, nmbWaves), dtype = np.ndarray)
    for i in range(0, nmbWaves):
        for j in range(0, nmbWaves):
            integralMatrix[i][j] = np.zeros(nmbTestFolders, dtype=complex)
            errorMatrixAbsolut[i][j] = np.zeros(nmbTestFolders, dtype=complex)
            errorMatrixRelativ[i][j] = np.zeros(nmbTestFolders, dtype=complex)
    '''

    for waveIndexI in range(0, 30): # nmbWaves):
        ampFile1 = uproot.open(ampPath + onlyfiles[waveIndexI])
        key1 = ampFile1.keys()[0]
        incohSubAmpsReal1 = ampFile1[key1]['amplitude'].get('_incohSubAmps._real').array()
        incohSubAmpsImag1 = ampFile1[key1]['amplitude'].get('_incohSubAmps._imag').array()
        ampMatrixI = prepareArray(np.array(incohSubAmpsReal1 + 1j * incohSubAmpsImag1))
        ampFile1 = None

        for waveIndexJ in range(0, waveIndexI + 1): #nmbWaves):
            if waveIndexI != waveIndexJ:
                ampFile2 = uproot.open(ampPath + onlyfiles[waveIndexJ])
                key2 = ampFile2.keys()[0]
                incohSubAmpsReal2 = ampFile2[key2]['amplitude'].get('_incohSubAmps._real').array()
                incohSubAmpsImag2 = ampFile2[key2]['amplitude'].get('_incohSubAmps._imag').array()
                ampMatrixJ = prepareArray(np.array(incohSubAmpsReal2 + 1j * incohSubAmpsImag2))
                ampFile2 = None
            else:
                ampMatrixJ = ampMatrixI
                
   
            integral = np.sum(ampMatrixI * ampMatrixJ.conjugate())/nmbEvents 
            sdReal = sum(pow((ampMatrixI * ampMatrixJ.conjugate()).real, 2)).real
            sdImag = sum(pow((ampMatrixI * ampMatrixJ.conjugate()).imag, 2)).real
            sdReal = math.sqrt(sdReal/float(nmbEvents) - pow(integral.real,2))
            sdImag = math.sqrt(sdImag/float(nmbEvents) - pow(integral.imag,2))
            sdComplex = complex(sdReal, sdImag)        
            errorComplex = complex(sdReal/math.sqrt(nmbEvents), sdImag/math.sqrt(nmbEvents))
            sd[waveIndexI][waveIndexJ] = sdComplex
            sd[waveIndexJ][waveIndexI] = sdComplex

            #print "%i, %i, %i: %s, %s, %s" %(folderI, waveIndexI, waveIndexJ, integral, sd[waveIndexI][waveIndexJ], errorComplex)


            integralMatrix[waveIndexI][waveIndexJ] = integral
            integralMatrix[waveIndexJ][waveIndexI] = integral
            errorMatrixAbsolut[waveIndexI][waveIndexJ] = errorComplex
            errorMatrixAbsolut[waveIndexJ][waveIndexI] = errorComplex
            if waveIndexI == waveIndexJ:
                errorMatrixRelativ[waveIndexI][waveIndexJ] = complex(errorComplex.real/math.fabs(integral.real), 0)
            else:
                errorMatrixRelativ[waveIndexI][waveIndexJ] = complex(errorComplex.real/math.fabs(integral.real), errorComplex.imag/math.fabs(integral.imag))
                errorMatrixRelativ[waveIndexJ][waveIndexI] = complex(errorComplex.real/math.fabs(integral.real), errorComplex.imag/math.fabs(integral.imag))

    return (folderI, integralMatrix, errorMatrixAbsolut, errorMatrixRelativ, testFolders[folderI])


# ## Berechnung

# In[37]:


nmbTestFolders = len(testFolders)
#nmbTestFolders = 16

ampPathTmp = mainPath + testFolders[0] + '/amps/'
onlyfiles = []
for filename in sorted(os.listdir(ampPathTmp)):
    onlyfiles.append(filename)
nmbWaves = len(onlyfiles)


#bar = progressbar.ProgressBar(maxval=nmbTestFolders, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
#bar.start()
pool = mp.Pool(mp.cpu_count()) # hercules = 8
#results = pool.map(calcErrors2, [folderI for folderI in range(0,nmbTestFolders)])


t0 = time.time()
results = []
with tqdm(total=nmbTestFolders) as pbar:
    #for i, res in tqdm(enumerate(pool.imap_unordered(calcErrors, [folderI for folderI in range(0,nmbTestFolders)]))):
    for i, res in tqdm(enumerate(pool.imap(calcErrors, [folderI for folderI in range(0,nmbTestFolders)]))):
        results.append(res) 
        pbar.update()

pool.close()   
pool.join()
#bar.update(bar.maxval)
t1 = time.time()


#nmbWaves = len([name for name in os.listdir(mainPath + testFolders[0] + '/amps/') if os.path.isfile(os.path.join(mainPath + testFolders[0] + '/amps/', name))])
axisArr = np.zeros(nmbTestFolders, dtype=int)
integralMatrix = np.zeros((nmbWaves, nmbWaves), dtype = np.ndarray)
errorMatrixAbsolut = np.zeros((nmbWaves, nmbWaves), dtype = np.ndarray)
errorMatrixRelativ = np.zeros((nmbWaves, nmbWaves), dtype = np.ndarray)
for i in range(0, nmbWaves):
    for j in range(0, nmbWaves):
        integralMatrix[i][j] = np.zeros(nmbTestFolders, dtype=complex)
        errorMatrixAbsolut[i][j] = np.zeros(nmbTestFolders, dtype=complex)
        errorMatrixRelativ[i][j] = np.zeros(nmbTestFolders, dtype=complex)
        
for resultI in range(0,len(results)):
    folderI = results[resultI][0]
    axisArr[folderI] = results[resultI][4]
    for waveIndexI in range(0, nmbWaves):
        for waveIndexJ in range(0, waveIndexI + 1):
            integralMatrix[waveIndexI][waveIndexJ][folderI] = results[resultI][1][waveIndexI][waveIndexJ]
            integralMatrix[waveIndexJ][waveIndexI][folderI] = results[resultI][1][waveIndexI][waveIndexJ]
            errorMatrixAbsolut[waveIndexI][waveIndexJ][folderI] = results[resultI][2][waveIndexI][waveIndexJ]
            errorMatrixAbsolut[waveIndexJ][waveIndexI][folderI] = results[resultI][2][waveIndexI][waveIndexJ]
            errorMatrixRelativ[waveIndexI][waveIndexJ][folderI] = results[resultI][3][waveIndexI][waveIndexJ]
            errorMatrixRelativ[waveIndexJ][waveIndexI][folderI] = results[resultI][3][waveIndexI][waveIndexJ]

'''        
for resultI in range(0,len(results)):
    folderI = results[resultI][0]
    axisArr[folderI] = results[resultI][4]
    for waveIndexI in range(0, nmbWaves):
        for waveIndexJ in range(0, waveIndexI + 1):
            integralMatrix[waveIndexI][waveIndexJ][folderI] = results[resultI][1][waveIndexI][waveIndexJ][folderI]
            integralMatrix[waveIndexJ][waveIndexI][folderI] = results[resultI][1][waveIndexI][waveIndexJ][folderI]
            errorMatrixAbsolut[waveIndexI][waveIndexJ][folderI] = results[resultI][2][waveIndexI][waveIndexJ][folderI]
            errorMatrixAbsolut[waveIndexJ][waveIndexI][folderI] = results[resultI][2][waveIndexI][waveIndexJ][folderI]
            errorMatrixRelativ[waveIndexI][waveIndexJ][folderI] = results[resultI][3][waveIndexI][waveIndexJ][folderI]
            errorMatrixRelativ[waveIndexJ][waveIndexI][folderI] = results[resultI][3][waveIndexI][waveIndexJ][folderI]
'''
            
print "Finish calculation successfully in %.2f min" %((t1 - t0)/60)


# In[ ]:

if args.filename:
	fileName = 'Fehler_bei_unterschiedlicher_Masse_Auswertung_' +  args.filename + '_Multi.npy'
else:
	fileName = 'Fehler_bei_unterschiedlicher_Masse_Auswertung_1541_1600_1_Multi.npy'
np.save(fileName, np.array([(axisArr, onlyfiles, integralMatrix, errorMatrixAbsolut, errorMatrixRelativ)]))
print "Saving data successfully to file:", fileName
print "Finish"


# ## Plot

# In[ ]:




