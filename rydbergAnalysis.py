#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 09:12:57 2021

@author: 
   ___  _ ____             __          
  / _ \(_) / /__ ___      / /  ___ ___ 
 / // / / / / -_) _ \    / /__/ -_) -_)
/____/_/_/_/\__/_//_/   /____/\__/\__/ 

"""
#import the necessary packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#path to file
path = '/home/dillen/University/Python/Rydberg Constant/'


#create some functions to clean data 
#remove the NaNs from the array
def clearNan(array):
    return [i for i in array if str(i) != 'nan']

#normalise to the central maxima angle and change into radians from degrees
def shiftAndRadians(array):
    return [np.sin((array[0]-i)*(np.pi/180)) for i in array]

def rydbergX(n):
    return ((1/n**2)-(1/4))





#skip all the unnecesary data
rydbergData = pd.read_csv(path+'rydberg.csv',header=0,usecols=[3,6,9,12])

red80,blue80,red300,blue300 = rydbergData.T.values



varNames = [red80,blue80,red300,blue300]

for i in range(4):
    varNames[i] = clearNan(varNames[i])
    varNames[i] = shiftAndRadians(varNames[i])
    
red80,blue80,red300,blue300 = varNames
dslits = [1/78.8,1/300]
wavelengths = []
errors = []
form = ['red','blue','orange','lightblue']


counter = 0
for var in varNames:
    #set the x values
    m = np.arange(0,len(var))
    
    #error in the reading of the vernier 
    errorInReading = [1.05e-5 for x in m]
    #curve fit
    fit,cov = np.polyfit(m,var,1 ,w=[1/x for x in errorInReading],cov=True)      
    linear = np.poly1d(fit)
    print('The gradient is %f with an uncertainty of %f'%(fit[0],np.sqrt(cov[0][0])))
    
    #find the lambda
    #lambda = d*sin(θ)/m which is d*gradient
    dInd = int(np.floor(counter/2))
    d = dslits[dInd]
    wavelength = d*fit[0]
    error = d*np.sqrt(cov[0][0])
    errors.append(error)
    print('The wavelength of this light is %.2f ± %f nm'%(wavelength*10**6,error*10**6))
    wavelengths.append(wavelength)
    
    #plot the graph
    plt.plot(m,linear(m),form[counter])
    plt.errorbar(m,var,yerr=errorInReading,fmt='xg',barsabove=True,capsize=3)
    plt.legend(['Curve fit with gradient %.2f'%fit[0],'Measured data'])
    plt.xlabel('Amount of peaks from central maximum')
    plt.title('Graph of wavelength %d nm through a %.1f slits/mm grating'%(wavelength*1e6,1/d))
    plt.ylabel('sin(θ)')
    
    counter += 1

    
    plt.grid()
    plt.show()
    
    
#calculate the rydberg constant
#first average out the two unique wavelengths from both slits
red = ((wavelengths[0]+wavelengths[2])/2)*1e-3
blue = ((wavelengths[1]+wavelengths[3])/2)*1e-3
redErrors = (np.sqrt(errors[0]**2+errors[2]**2)/2)*1e-3
blueErrors = (np.sqrt(errors[1]**2+errors[3]**2)/2)*1e-3

wavelengths = [1/red,1/blue]
totalErrors = [redErrors,blueErrors]
#now I'm going to take another assumption and assume red wavelength comes from deexcitations from n = 3 and blue, n = 4
xVals = [rydbergX(3),rydbergX(4)]
print('average red wavelength is %f ± %f nm'%(red*1e9,redErrors*1e9))
print('average blue wavelength is %f ± %f nm'%(blue*1e9,blueErrors*1e9))

rydberg = -(wavelengths[1]-wavelengths[0])/(xVals[1]-xVals[0])
#this error analysis looks pretty complicated but really it is adding the errors of the error in each wavelength in quadrature 
#i.e R = 1/lambdaRed*constant - 1/lambdaBlue*constant 
#so first calculate the red error then the blue and finally add in quadrature 
rydbergError = np.sqrt((redErrors/((xVals[1]-xVals[0])*(red**2)))**2+(blueErrors/((xVals[1]-xVals[0])*(blue**2)))**2)



plt.errorbar(xVals,wavelengths,yerr=[redErrors,blueErrors],fmt='|r--',capsize=5,label='Gradient is %f'%(-rydberg))
plt.xlabel(r'$ \frac{1}{n^2} - \frac{1}{2^2}$')
plt.ylabel(r'1/λ ($m^{-1}$)')
plt.title('Linear plot to find Rydberg constant')
plt.legend(loc='upper center')
plt.grid()
plt.show()
print('Rydberg constant is %d ± %f'%(rydberg,rydbergError))






