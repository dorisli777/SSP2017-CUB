from math import *
from vpython import *

def centerOfMass(fileName):
    
    myFile = open(fileName,'r')
    myArr = []
    nbrRows = 0
    for line in myFile:
        myArr.append(line)
        nbrRows += 1
    for row in range (0,nbrRows):
        myArr[row] = myArr[row].strip() #takes off the /n
        myArr[row] = myArr[row].split() #splits the lines into an array
    for i in range (0,nbrRows):
        for j in range (0,nbrRows):
            myArr[i][j] = float(myArr[i][j]) #from string to float
    
    origin = input("Enter the center pixel (row, col): ")
    n = input("Enter the dimension of pixel: ")

    #calculation for center of mass
    xCm = 0
    yCm = 0
    sumMass = 0
    #loops through each point in the NxN matrix 
    for i in range (origin[1]- n/2, origin[1] + n/2 +1):
        for j in range(origin[0]- n/2, origin[0] + n/2 +1):
            xCm += (j - origin[1]) * myArr[i][j]
            yCm += (origin[0] - i) * myArr[i][j]
            sumMass += myArr[i][j]
    
    xCm2 = xCm/sumMass
    yCm2 = yCm/sumMass
    
    print "Xcm: ", xCm2
    print "Ycm: ", yCm2
    
    #calculation for uncertainty
    uncertainX = 0
    uncertainY = 0
    #loops through each point in the NxN matrix 
    for i in range(origin[1] - n/2, origin[1] + n/2 +1):
        for j in range(origin[0] - n/2, origin[0] + n/2 +1):
            uncertainX += ((j - origin[1] - yCm2)**2)*myArr[i][j]
            uncertainY += ((origin[1] - i - xCm2)**2)*myArr[i][j]
    uncertainX = 1/float(sumMass)*sqrt(uncertainX)
    uncertainY = 1/float(sumMass)*sqrt(uncertainY)
    
    print "Uncertainty of x-coordinate: ", uncertainX
    print "uncertainty of y-coordinate: ", uncertainY

    #visualization
    for i in range(0,nbrRows):
        for j in range(nbrRows):
            ball = sphere(pos=vector(j,n-1-i,0), radius=myArr[i][j]*0.001, color=color.red)
    
    #finding centroid
    x = origin[1] + xCm2
    y = origin[0] + yCm2
    ball = sphere(pos=vector(x,y,0), radius=0.05, color=color.green)
            
centerOfMass("datafile.txt")

