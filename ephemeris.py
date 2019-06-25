# This function generates an ephemeris for an asteroid given the inputs as directed below 
# Last modified August 2017 at the Summer Science Program in Boulder, CO

from math import *

#import magnitude finder 
def magnitude_x(x):
    y = 0
    for i in x:
        y = y + i**2 
    return sqrt(y) 

#import angle checker
def angle_quad(sine, cosine):
    angle1 = float(acos(cosine))
    if sine >= 0 and cosine >= 0:
        angle = angle1
    if sine < 0 and cosine >= 0:
        angle = (2*pi) - angle1
    if sine >= 0 and cosine < 0:
        angle = angle1
    if sine < 0 and cosine < 0:
        angle = pi - angle1
    return angle

def ephemeris(e, a, i1, bOmega1, sOmega1, m01, t):
    #change values into radians
    i = i1 * (pi/180)
    bOmega = bOmega1 * (pi/180)
    sOmega = sOmega1 * (pi/180)
    m0 = m01 * (pi/180)

    #find mean anomaly (M)
    k = 0.01720209894
    tE = 2458000.5
    M = (k*sqrt(1/a**3)*(t - tE)) + m0
    print "Mean Anomaly (M): ", M

    #find eccentric anomaly (E, vertical from orbit if circumscribe)
    E = M
    for i0 in range(5):
        E = E - ((E-e*sin(E)-M)/(1-e*cos(E)))
    print "Eccentric Anomaly (E): ", E

    #find cartesian coordinates (orbital plane)
    xCar = a*(cos(E)-e)
    yCar = a*sqrt(1-e**2)*sin(E)
    zCar = 0
    rCar = [xCar, yCar, zCar]
    print "Cartesian coordinates: ", rCar

    #matrix transformation (argument of perihelion)
    xEclw = (xCar*cos(sOmega)) - (yCar*sin(sOmega))
    yEclw = (xCar*sin(sOmega)) + (yCar*cos(sOmega))
    zEclw = zCar
    rEclw = [xEclw, yEclw, zEclw]
    print "Small Omega: ", rEclw
    
    #matrix transformation (inclination)
    xEcli = xEclw
    yEcli = yEclw*cos(i) - zEclw*sin(i)
    zEcli = yEclw*sin(i) + zEclw*cos(i)
    rEcli = [xEcli, yEcli, zEcli]
    print "Inclination: ", rEcli
    
    #matrix transformation (longitude of asceding node)
    xEcl = xEcli*cos(bOmega) - yEcli*sin(bOmega)
    yEcl = xEcli*sin(bOmega) + yEcli*cos(bOmega)
    zEcl = zEcli
    rEcl = [xEcl, yEcl, zEcl]
    print "Big Omega: ", rEcl
    
    #matrix transformation (epsilon)
    epsilon = 0.4090926277
    xEq = xEcl
    yEq = yEcl*cos(epsilon) - zEcl*sin(epsilon)
    zEq = yEcl*sin(epsilon) + zEcl*cos(epsilon)
    rEq = [xEq, yEq, zEq]
    print "Epsilon: ", rEq
    
    #find range of vector (p) using Sun vector
    R = [-0.404120107345586, 0.855666323896357, 0.370934219454331]
    p = [rEq[0] + R[0], rEq[1] + R[1], rEq[2] + R[2]]
    print "Range of vector: ", p
    
    #find p hat
    pMag = float(magnitude_x(p))
    pHat = [p[0]/pMag, p[1]/pMag, p[2]/pMag]
    print "Rho Hat: ", pHat
    
    #find declination and ra 
    dec = asin(pHat[2])
    decDeg = dec * (180/pi)
    cosRa = pHat[0]/cos(dec)
    sinRa = pHat[1]/cos(dec)
    ra = angle_quad(sinRa, cosRa) * (180/pi)
    print "Declination: ", dec, "RA: ", ra
    
    #convert dec into deg, min, sec
    decDeg1 = int(decDeg)
    decMin = int((decDeg - decDeg1) * 60)
    decSec = (60*(decDeg - decDeg1) - decMin)*60
    print "DecDeg: ", decDeg1, "DecMin: ", decMin, "DecSec: ", decSec
    
    #convert ra into hr, min, sec
    raHr = int(ra/15)
    raMin = int(60 * (ra/15 - raHr))
    raSec = (60 * (ra/15 - raHr) - raMin) * 60
    print "RA Hr: ", raHr, "RA Min: ", raMin, "RA Sec: ", raSec
    
# a test input
ephemeris(0.7654035254, 2.236572976, 4.959459269, 338.8496637, 168.1622182, 43.510113400153, 2457957.7190240743)
