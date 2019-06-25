# This function gives the six orbital elements of an asteroid given a data file containing RA and DEC.  

from math import *

#magnitude function
def mag(x):
    y = 0
    for i in x:
        y = y + i**2 
    return sqrt(y) 

#dot product function
def dot(x,y):
    n = 0
    i = 0
    while i <= (len(x)-1):
        n = x[i]*y[i] + n 
        i = i + 1
    return n

#cross product function
def cross(x,y):
    g = x[1]*y[2] - x[2]*y[1]
    h = x[0]*y[2] - x[2]*y[0]
    i = x[0]*y[1] - x[1]*y[0]
    return g,-h,i

#angle checker
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

#find orbital elements function
def orbitalElements(r, rDot):
    r0 = sqrt(r[0]**2 + r[1]**2 + r[2]**2)
    v02 = dot(rDot,rDot)
    a = 1/(2/r0 - v02)
    print "Semi-Major Axis (a): ", a
    
    rDot2 = (mag(cross(r,rDot)))**2
    e = sqrt(1-(rDot2/a))
    print "Eccentricity (e): ", e
    
    hVec = cross(r,rDot)
    i = atan(sqrt(hVec[0]**2+hVec[1]**2)/hVec[2]) * (180/pi)
    print "Inclination (i): ", i
    
    bigOmega = acos(-hVec[1]/(mag(hVec)*sin(i*pi/180))) * (180/pi)
    print "Longitude of Ascending Node (big omega): ", bigOmega
    
    cosv0 = (1/e)*((a*(1-e**2))/r0 - 1)
    sinv0 = a*(1-e**2) / (e*mag(hVec))*dot(r,rDot)/r0

    #define u (angle between ascending node and asteroid)
    cosu0 = (r[0]*cos(bigOmega*(pi/180)) + r[1]*sin(bigOmega*pi/180))/r0
    sinu01 = r[2]/(r0*sin(i*pi/180))
    #check the sign of u0
    u0 = angle_quad(sinu01,cosu0)
    v0 = angle_quad(sinv0,cosv0)

    #find argument of perihelion (angle between ascending node and perihelion)
    w = (u0 - v0) % (2 * pi)
    print "Argument of Perihelion (small omega): ", w * (180/pi)
    
    E = acos(1/e*(1-mag(r)/a))
    M = E - e*sin(E)
    print "Mean Anomaly (M): ", M*(180/pi)
    
################################################
def orbitdet(fileName):

    myFile = open(fileName,'r') #open file
    myArr = []
    num_rows = 0
    #create array from file
    for line in myFile:
        myArr.append(line)
        num_rows +=1
    for row in range(0,num_rows):
        myArr[row] = myArr[row].strip() #takes off the /n
        myArr[row] = myArr[row].split() #splits the lines into an array
    for i in range(0,3):
        for j in range(0,13):
            myArr[i][j] = float(myArr[i][j]) #change all the string elements in the array to floats
    
    ra = []
    dec = []
    #changing ra and dec to decimal degrees
    for i in range(0,3):
        if myArr[i][0] > 0:
            ra.append(radians((myArr[i][0]*15.) + (myArr[i][1]/4.) + (myArr[i][2]/240.)))
        elif myArr[i][0] < 0: 
            ra.append(radians((myArr[i][0]*15.) - (myArr[i][1]/4.) - (myArr[i][2]/240.))) 
        if myArr[i][3] > 0:
            dec.append(radians((myArr[i][3] + (myArr[i][4]/60.) + (myArr[i][5]/3600.)))) 
        elif myArr[i][3] < 0:
            dec.append((radians(myArr[i][3] - (myArr[i][4]/60.) - (myArr[i][5]/3600.))))
    
    observation1 = myArr[0] #seperate lists for each observation
    observation2 = myArr[1]
    observation3 = myArr[2]
    
    time1 = observation1[6]  #in JD
    time2 = observation2[6]
    time3 = observation3[6]
    
    RA1 = ra[0] #decimal hours
    RA2 = ra[1] 
    RA3 = ra[2]

    dec1 = dec[0] #decimal degrees
    dec2 = dec[1] 
    dec3 = dec[2]

    R1 = [observation1[7],observation1[8],observation1[9]] #sun vectors
    R2 = [observation2[7],observation2[8],observation2[9]]
    R3 = [observation3[7],observation3[8],observation3[9]]
    
    cR1 = [observation1[7],observation1[8],observation1[9]] #sun vectors
    cR2 = [observation2[7],observation2[8],observation2[9]]
    cR3 = [observation3[7],observation3[8],observation3[9]]
    
    vR1 = [observation1[10], observation1[11], observation1[12]] #sun velocity vectors
    vR2 = [observation2[10], observation2[11], observation2[12]]
    vR3 = [observation3[10], observation3[11], observation3[12]]
    
    ##find rho hat vector### (in AU)
    p_hat1 = [cos(dec1)*cos(RA1), cos(dec1)*sin(RA1), sin(dec1)] 
    p_hat2 = [cos(dec2)*cos(RA2), cos(dec2)*sin(RA2), sin(dec2)]
    p_hat3 = [cos(dec3)*cos(RA3), cos(dec3)*sin(RA3), sin(dec3)]

    ##find tau##
    k = 0.01720209895
    tau1 = k*(time1 - time2)
    tau2 = k*(time3 - time1)
    tau3 = k*(time3 - time2)
    
    ##find initial constants a1 and a3##
    a1_ini = (time3 - time2)/(time3 - time1)
    a3_ini = -(time1 - time2)/(time3 - time1) 
    
    #hairy triple vector products
    D1 = dot(cross(p_hat1,p_hat2), p_hat3)
    D2 = dot(cross(p_hat2,p_hat3), p_hat1)
    D01 = dot(cross(R1, p_hat2), p_hat3)
    D02 = dot(cross(R2, p_hat2), p_hat3)
    D03 = dot(cross(R3, p_hat2), p_hat3)
    D11 = dot(cross(p_hat1, R1), p_hat3)
    D12 = dot(cross(p_hat1, R2), p_hat3)
    D13 = dot(cross(p_hat1, R3), p_hat3)
    D21 = dot(cross(p_hat2, R1), p_hat1)
    D22 = dot(cross(p_hat2, R2), p_hat1)
    D23 = dot(cross(p_hat2, R3), p_hat1)
    
    ##find rho vector magnitudes##
    pmag1_ini= (a1_ini*D01 - D02 + a3_ini*D03) / (a1_ini*D2)
    pmag2_ini = -(a1_ini*D11 - D12 + a3_ini*D13) / D2
    pmag3_ini = (a1_ini*D21 - D22 + a3_ini*D23) / (a3_ini*D1)
    
    ##find rho vectors##
    p1_ini = [pmag1_ini*p_hat1[0], pmag1_ini*p_hat1[1], pmag1_ini*p_hat1[2]]
    p2_ini = [pmag2_ini*p_hat2[0], pmag2_ini*p_hat2[1], pmag2_ini*p_hat2[2]]
    p3_ini = [pmag3_ini*p_hat3[0], pmag3_ini*p_hat3[1], pmag3_ini*p_hat3[2]]
    
    ######LIGHT CORRECTION#####
    c = 173.1446 #in AU per day
    time1new = time1 - pmag1_ini/c
    time2new = time2 - pmag2_ini/c
    time3new = time3 - pmag3_ini/c
    
    #another way to write pmag/c (change in time)
    x1 = time1new - time1
    x2 = time2new - time2
    x3 = time3new - time3
    
    #correct for sun vector
    cR1 = [(R1[0] + x1*vR1[0]), R1[1] + x1*vR1[1], R1[2] + x1*vR1[2]]
    cR2 = [(R2[0] + x2*vR2[0]), R2[1] + x2*vR2[1], R2[2] + x2*vR2[2]]
    cR3 = [(R3[0] + x3*vR3[0]), R3[1] + x3*vR3[1], R3[2] + x3*vR3[2]]
    
    #correct tau values
    tau1 = k*(time1new - time2new)
    tau2 = k*(time3new - time1new)
    tau3 = k*(time3new - time2new)
    
    ##find position vectors (r = p - R) ###
    r1_ini = [p1_ini[0] - cR1[0], p1_ini[1] - cR1[1], p1_ini[2] - cR1[2]]
    r2_ini = [p2_ini[0] - cR2[0], p2_ini[1] - cR2[1], p2_ini[2] - cR2[2]]  #r0
    r3_ini = [p3_ini[0] - cR3[0], p3_ini[1] - cR3[1], p3_ini[2] - cR3[2]]
    magr2_i = mag(r2_ini)
    
    #find r dot vector 
    rdot_ini = [(r3_ini[0] - r1_ini[0])/tau2, (r3_ini[1] - r1_ini[1])/tau2, (r3_ini[2] - r1_ini[2])/tau2,]
    
    ##f and g initial guesses##
    f1_ini = 1 - (tau1**2)/(2*magr2_i**3) + dot(rdot_ini, r2_ini)*(tau1**3)/(2*magr2_i**5) 
    g1_ini = tau1 - ((tau1**3)/(6*(magr2_i**3)))  
    f3_ini = 1 - (tau3**2)/(2*magr2_i**3) + dot(rdot_ini, r2_ini)*(tau3**3)/(2 * magr2_i**5)
    g3_ini = tau3 - ((tau3**3)/(6*(magr2_i**3)))  
    
    #find new values of constants 
    a1 = g3_ini/(f1_ini*g3_ini - f3_ini*g1_ini)
    a3 = -g1_ini/(f1_ini*g3_ini - f3_ini*g1_ini)
    
    #find new magnitudes with new a constants 
    pmag1 = (a1*D01 - D02 + a3*D03)/(a1*D2)
    pmag2 = -(a1*D11 - D12 + a3*D13)/D2
    pmag3 = (a1*D21 - D22 + a3*D23)/(a3*D1)

    #new rho vectors 
    p1 = [pmag1*p_hat1[0], pmag1*p_hat1[1], pmag1*p_hat1[2]]
    p2 = [pmag2*p_hat2[0], pmag2*p_hat2[1], pmag2*p_hat2[2]]
    p3 = [pmag3*p_hat3[0], pmag3*p_hat3[1], pmag3*p_hat3[2]]
    
    #######LIGHT CORRRECTION#######
    time1new = time1 - pmag1/c
    time2new = time2 - pmag2/c
    time3new = time3 - pmag3/c
    
    #another way to write pmag/c (time difference)
    x1 = time1new - time1
    x2 = time2new - time2
    x3 = time3new - time3
    
    #correct for sun vector
    cR1 = [(R1[0] + x1*vR1[0]), R1[1] + x1*vR1[1], R1[2] + x1*vR1[2]]
    cR2 = [(R2[0] + x2*vR2[0]), R2[1] + x2*vR2[1], R2[2] + x2*vR2[2]]
    cR3 = [(R3[0] + x3*vR3[0]), R3[1] + x3*vR3[1], R3[2] + x3*vR3[2]]
    
    #new tau values
    tau1 = k*(time1new - time2new)
    tau2 = k*(time3new - time1new)
    tau3 = k*(time3new - time2new)
    
    ##find position vectors ##
    r1 = [p1[0] - cR1[0], p1[1] - cR1[1], p1[2] - cR1[2]]
    r2 = [p2[0] - cR2[0], p2[1] - cR2[1], p2[2] - cR2[2]]  #r0
    r3 = [p3[0] - cR3[0], p3[1] - cR3[1], p3[2] - cR3[2]]
    magr2 = mag(r2)
    
    #find r dot vectors
    rdot = [(r3[0] - r1[0])/tau2, (r3[1] - r1[1])/tau2, (r3[2] - r1[2])/tau2]
     
    #another f and g series to make values of r and r dot more precise 
    f1 = 1 - (tau1**2)/(2*magr2**3) + dot(rdot, r2)*(tau1**3)/(2*magr2**5) #f is dimensionless
    g1 = tau1 - ((tau1**3)/(6 * (magr2**3)))  #g has units of time
    f3 = 1 - (tau3**2)/(2*magr2**3) + dot(rdot, r2)*(tau3**3)/(2*magr2**5)
    g3 = tau3 - ((tau3**3)/(6*(magr2**3)))  
    
    #find r dot vector
    c1 = f3/(g1*f3 - g3*f1)
    c2 = f1/(g1*f3 - g3*f1)
    rdot = [c1*r1[0] - c2*r3[0], c1*r1[1] - c2*r3[1], c1*r1[2] - c2*r3[2]] 

    ##iteration loop for r and rdot vectors##
    r1f = [] #set empty final lists for r position vector
    r3f = []
    count = 0 #set a counter for number of iterations
    while abs(magr2 - mag(r2_ini)) != 0:
        count += 1
        r2_ini = r2 #sets the new value to the "old" value
        
        #f and g series 
        f1 = 1 - (tau1**2)/(2*magr2**3) + dot(rdot, r2)*(tau1**3)/(2*magr2**5) #f is dimensionless
        g1 = tau1 - ((tau1**3)/(6*(magr2**3)))  #g has units of time
        f3 = 1 - (tau3**2)/(2*magr2**3) + dot(rdot, r2)*(tau3**3)/(2*magr2**5)
        g3 = tau3 - ((tau3**3)/(6*(magr2**3))) 
        
        #new a constants 
        a1 = g3/(f1*g3 - f3*g1)
        a3 = -g1/(f1*g3 - f3*g1)
        
        #new magnitude vectors 
        pmag1 = (a1*D01 - D02 + a3*D03)/(a1*D2)
        pmag2 = -(a1*D11 - D12 + a3*D13)/D2
        pmag3 = (a1*D21 - D22 + a3*D23)/(a3*D1)
        
        #########rho vectors########
        p1 = [pmag1*p_hat1[0], pmag1*p_hat1[1], pmag1*p_hat1[2]]
        p2 = [pmag2*p_hat2[0], pmag2*p_hat2[1], pmag2*p_hat2[2]]
        p3 = [pmag3*p_hat3[0], pmag3*p_hat3[1], pmag3*p_hat3[2]]
        
        time1new = time1 - pmag1/c
        time2new = time2 - pmag2/c
        time3new = time3 - pmag3/c
        
        #another way to write pmag/c (change in time)
        x1 = time1new - time1
        x2 = time2new - time2
        x3 = time3new - time3

        #correct for sun vector
        cR1 = [(R1[0] + x1*vR1[0]), R1[1] + x1*vR1[1], R1[2] + x1*vR1[2]]
        cR2 = [(R2[0] + x2*vR2[0]), R2[1] + x2*vR2[1], R2[2] + x2*vR2[2]]
        cR3 = [(R3[0] + x3*vR3[0]), R3[1] + x3*vR3[1], R3[2] + x3*vR3[2]]

        tau1 = k*(time1new - time2new)
        tau2 = k*(time3new - time1new)
        tau3 = k*(time3new - time2new)
        #####END OF LIGHT #######

        #find position vectors 
        r1f = [p1[0] - cR1[0], p1[1] - cR1[1], p1[2] - cR1[2]]
        r2 = [p2[0] - cR2[0], p2[1] - cR2[1], p2[2] - cR2[2]] 
        r3f = [p3[0] - cR3[0], p3[1] - cR3[1], p3[2] - cR3[2]]
        magr2 = mag(r2)
        
        #find r dot vector
        c1 = f3/(g1*f3 - g3*f1)
        c2 = f1/(g1*f3 - g3*f1)
        rdot = [c1*r1f[0] - c2*r3f[0], c1*r1f[1] - c2*r3f[1], c1*r1f[2] - c2*r3f[2]]
        
    print "r2 Vector: ", r2
    print "rdot Vector: ", rdot
    print "Number of iterations: ", count
    
    #convert from equatorial to ecliptic
    epsilon = 0.4090926277
    #rotation matrix for r position vector 
    xEcl = r2[0]
    yEcl = r2[1]*cos(-epsilon) - r2[2]*sin(-epsilon)
    zEcl = r2[1]*sin(-epsilon) + r2[2]*cos(-epsilon)
    rEcl = [xEcl, yEcl, zEcl]
    
    #rotation matrix for r dot (velocity) vector 
    xEcl1 = rdot[0]
    yEcl1 = rdot[1]*cos(-epsilon) - rdot[2]*sin(-epsilon)
    zEcl1 = rdot[1]*sin(-epsilon) + rdot[2]*cos(-epsilon)
    rdotEcl = [xEcl1, yEcl1, zEcl1]
    print rEcl, rdotEcl
    
    #call orbital elements function
    orbitalElements(rEcl, rdotEcl)

#own asteroid 
orbitdet("asterorbit.txt")

#test data
#orbitdet("orbitdet.txt")
