#! /usr/bin/ipython

#########################################################################################################################
#	cloudObject is a python script designed to convert a data file to a cloud object for use with SAMAC software	#
#	The data file is a parameter											#
#	The script is executed with the command "ipython cloudObject.py"						#
#	Created June 2016												#
#	Copyright 2016 Alex Butland											#
#########################################################################################################################


import numpy as np
import csv
import Cloud  #SAMAC object
import pickle

f=open('BasicAircraftData.dat','r') 					#user imports data into python
lines=f.readlines()
f.close()

for i in range(len(lines)):         					#removing end of line remarks
	lines[i]=lines[i].strip()

titles=list(map(str,lines[0].split(',')))				#getting titles and units from file and putting them into lists
units=list(map(str,lines[1].split(',')))

df=list()

df=np.loadtxt('BasicAircraftData.dat',skiprows=2,delimiter=',')		#loading data as floats into a numpy array, skip 2 rows of header, commas as delimiters
data=np.ma.array(df)							#formatting the data as a masked array

times={									#creating an empty times of maneuvers box as a dictionary
"abovecloud": np.zeros((0,2)),
"belowcloud": np.zeros((0,2)),
"cloud": np.zeros((0,2)),
"horicloud": np.zeros((0,2)),
"verticloud": np.zeros((0,2))
}

sizedist=list()								#create an empty size distribution box

ExampleCloud=Cloud.Cloud(data.transpose(),titles,units,times,sizedist)		#creating a cloud object with the data we have

ExampleCloud.times["cloud"]=np.array([[ExampleCloud.data[0][0], ExampleCloud.data[0][1945]]])		#defining the start and end times of the cloud object

pickle.dump(ExampleCloud,open("ExampleCloud.p","wb"))				#saving the object using pickle

ExCloud=pickle.load(open("/usr/local/SAMAC/tags/0.9.2/ExampleCloud.p","rb"))	#cloud object


