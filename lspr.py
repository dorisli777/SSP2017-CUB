# This function finds the RA and DEC of an asteroid given the known coordinates of a set of stars in the image.

from math import *

#import det3x3 and other solvers
def det3x3(x):
    a = float(x[0][0])*(x[1][1]*x[2][2] - x[1][2]*x[2][1])
    b = float(x[0][1])*(x[1][0]*x[2][2] - x[1][2]*x[2][0])
    c = float(x[0][2])*(x[1][0]*x[2][1] - x[1][1]*x[2][0])
    return a - b + c

def solve3eq(coeff, const):
    a = [[const[0],coeff[0][1],coeff[0][2]],[const[1],coeff[1][1],coeff[1][2]],[const[2],coeff[2][1],coeff[2][2]]]
    b = [[coeff[0][0],const[0],coeff[0][2]],[coeff[1][0],const[1],coeff[1][2]],[coeff[2][0],const[2],coeff[2][2]]]
    c = [[coeff[0][0],coeff[0][1],const[0]],[coeff[1][0],coeff[1][1],const[1]],[coeff[2][0],coeff[2][1],const[2]]]
    return [float(det3x3(a)) / det3x3(coeff), float(det3x3(b)) / det3x3(coeff), float(det3x3(c)) / det3x3(coeff)]

#setting up an array
myFile = open("lspr.txt",'r')
myArr = []
nbrRows = 0
for line in myFile:
    myArr.append(line)
    nbrRows += 1
for row in range (0,nbrRows):
    myArr[row] = myArr[row].strip() #takes off the /n
    myArr[row] = myArr[row].split() #splits the lines into an array

for i in range (0,6):
    for j in range (0,8):
        myArr[i][j] = float(myArr[i][j]) #from string to float

#part2 (decimalizes ra and dec)
x = []
y = []
ra = []
dec = []

for i in range (0,6):
    x.append(myArr[i][0])
    y.append(myArr[i][1])
    if (myArr[i][2] >= 0):
        ra.append(myArr[i][2]*15. + myArr[i][3]/4. + myArr[i][4]/240.)
    else: #takes into account negative ra
        ra.append(myArr[i][2]*15. - myArr[i][3]/4. - myArr[i][4]/240.)
    if (myArr[i][5] >= 0):
        dec.append(myArr[i][5] + myArr[i][6]/60. + myArr[i][7]/3600.)
    else: #takes into account negative dec
        dec.append(myArr[i][5] - myArr[i][6]/60. - myArr[i][7]/3600.)

#part3 (best plate coefficients)
N = 6.
sumx = 0
sumy = 0
sumxy = 0
sumx2 = 0
sumy2 = 0
sumra = 0
sumrax = 0
sumray = 0
sumdec = 0
sumdecx = 0
sumdecy = 0

#setting up variables
for i in range (0,6):
    sumx += x[i]
    sumy += y[i]
    sumxy += x[i]*y[i]
    sumx2 += x[i]*x[i]
    sumy2 += y[i]*y[i]
    sumra += ra[i]
    sumrax += ra[i]*x[i]
    sumray += ra[i]*y[i]
    sumdec += dec[i]
    sumdecx += dec[i]*x[i]
    sumdecy += dec[i]*y[i]

#set up for cramer's rule
coeff = [[N,sumx,sumy],[sumx,sumx2,sumxy],[sumy,sumxy,sumy2]]
ra1 = [sumra,sumrax,sumray]
dec1 = [sumdec,sumdecx,sumdecy]

#using solver frm the beginning
rA = solve3eq(coeff, ra1)
dEC = solve3eq(coeff, dec1)

#assigning coefficient values (finishes cramer's rule)
b1 = rA[0]
b2 = dEC[0]
a11 = rA[1]
a21 = dEC[1]
a12 = rA[2]
a22 = dEC[2]

print "Plate Coefficients: ", b1,b2,a11,a12,a21,a22

#finding residuals
rRA = []
rDEC = []
for i in range (0,6):
    rRA.append((ra[i] - (b1 + a11*x[i] + a12*y[i]))*3600.)
    rDEC.append((dec[i] - (b2 + a21*x[i] + a22*y[i]))*3600.)

for i in range(0,6):
    print "Star", i+1, ":", rRA[i],rDEC[i]

#Begin standard deviation of ra and dec
sumrRA = 0
sumrDEC = 0

#plug into the SD formula
for i in range(0,6):
    sumrRA += rRA[i]*rRA[i]
    sumrDEC += rDEC[i]*rDEC[i]

sdRA = sqrt(1/(N-3.)*sumrRA)
sdDEC = sqrt(1/(N-3.)*sumrDEC)

print "SD of RA: ", sdRA, "SD of DEC: ", sdDEC

#finding asteroid ra and dec
asterx = 57.373320
astery = 22.486560
asterRA = b1 + a11*asterx + a12*astery
asterDEC = b2 + a21*asterx + a22*astery

print "Asteroid RA(degrees): ", asterRA, "+-", sdRA
print "Asteroid DEC(degrees): ", asterDEC, "+-", sdDEC

#converting from decimal degrees to ra in hr, min, sec
h = asterRA
hr = h/15.
minn = abs(hr-int(hr))*60.
sec = abs((hr-int(hr))*60. - int(minn))*60.
print "RA = hr:", int(hr), "min:", int(minn), "sec:", sec

#converting from dec deg to dec in deg, min, sec
g = asterDEC
deg = int(g)
minn1 = abs(int((g-deg)*60.))
sec1 = abs((g-deg)*60. - int((g-deg)*60.))*60.
print "DEC = deg:", deg, "min:",minn1, "sec:", sec1

