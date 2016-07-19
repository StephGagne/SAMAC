#! /usr/bin/ipython

"""
   NAME:
     LoadDataExample.py

   PURPOSE:
     LoadDataExample is a python script designed to convert a data file to a cloud object for use with SAMAC software

   CALLS:
     Nothing

   MODIFICATIONS:
     Alex Butland   - 160712:  Modification

   USAGE:
     /usr/local/SAMAC/tags/0.9.2/LoadDataExample.py

   COPYRIGHT NOTICE:
     Copyright 2016 Alex Butland                                                                                  
"""


import numpy as np
import csv
import sys
#SAMAC object
import Cloud  
import pickle

infile = sys.argv[-1]

#Opening the file and storing its data. 
f = open(infile,'r') 
lines = f.readlines()
f.close()

#Removing end of line remarks.
for i in range(len(lines)):         					
    lines[i] = lines[i].strip()

#Getting titles and units from data and putting them into lists.
titles = list(map(str,lines[0].split(',')))				
units = list(map(str,lines[1].split(',')))

df = list()

#Loading data as floats into a numpy array, skipping 2 rows of titles and units.
#Commas used as delimiters between data.
df = np.loadtxt(infile,skiprows = 2,delimiter = ',')	

#Formatting the data as a masked array.	
data = np.ma.array(df)							

#Creating an empty times of maneuvers box as a dictionary.
times = {									
"abovecloud": np.zeros((0,2)),
"belowcloud": np.zeros((0,2)),
"cloud": np.zeros((0,2)),
"horicloud": np.zeros((0,2)),
"verticloud": np.zeros((0,2))
}

#Create an empty size distribution box.
sizedist = list()								

#Creating a cloud object with the data we have.
ExampleCloud = Cloud.Cloud(data.transpose(),titles,units,times,sizedist)		

#Defining the start and end times of the cloud object.
ExampleCloud.times["cloud"] = np.array([[ExampleCloud.data[0][0], ExampleCloud.data[0][1945]]])		

#Saving the object using pickle.
pickle.dump(ExampleCloud,open("ExampleCloud.p","wb"))				

#Cloud object.
ExCloud = pickle.load(open("/usr/local/SAMAC/tags/0.9.2/ExampleCloud.p","rb"))	
