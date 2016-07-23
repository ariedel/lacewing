import numpy as np
import time
import kinematics
#from matplotlib import pyplot as plt

#  Functions in this module:
# fitellipse: Fits a 3D ellipsoid by performing linear fits to successive planes of the 3D array, derotating the ellipsoid each time. 
# rotate: Constructs the rotation matrix from fitellipse using the input Tait-Bryant angles.
# fitellipse2d: A modified version of fitellipse for two dimensions.  Useful if you want to make plots of specific planes (UV,UW,VW, for example) because the pices of an ellipse are not interchangeable.
# triaxial: Takes the dictionary produced by fitellipse, returns an object that can be plotted by mplot3D's surface plotter.


# This function fits a 9 parameter ellipse by making linear fits on an axis-by-axis basis.
# The input values are numpy arrays of points to fit.
# This function is intended to be run in a monte carlo simulation.
def fitellipse(x,y,z):

    meanx = np.mean(x)
    meany = np.mean(y)
    meanz = np.mean(z)
    
    vec0 = np.asarray([x,y,z])

    # 1: determine XY rotation angle
    coeff = np.polyfit(vec0[0],vec0[1],1)
    anglexy = np.arctan2(coeff[0],1)
    
    cosc = np.cos(anglexy)
    sinc = np.sin(anglexy)
    rotz_xy = np.array([ [cosc,sinc,0],
                            [-sinc,cosc,0],
                            [0,0,1] ])
    
    vec1 = np.dot(rotz_xy,vec0)
    
    # 2: determine XZ 
    
    m = 1
    b = 1
    coeff = np.polyfit(vec1[0],vec1[2],1)
    
    anglexz = np.arctan2(coeff[0],1)
    
    cosb = np.cos(anglexz) # sazzle frazzle
    sinb = np.sin(anglexz)
    roty_xz = np.array([ [cosb,0, sinb],
                            [0,1,0],
                            [-sinb,0, cosb] ])
    
    vec2 = np.dot(roty_xz,vec1)
    
    # 3: determine YZ
    
    m = 1
    b = 1
    coeff = np.polyfit(vec2[1],vec2[2],1)
    
    angleyz = np.arctan2(coeff[0],1)
    cosa = np.cos(angleyz)
    sina = np.sin(angleyz)
    rotx_yz = np.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])
    
    vec3 = np.dot(rotx_yz,vec2)
    
    #now rotate and translate
    # yz rotates around the x axis, xz rotates around the y axis, and xy rotates around the z axis

    rotmatrix = np.dot(rotx_yz,np.dot(roty_xz,rotz_xy))

    # sanity check: If this was done correctly, vec3 should be identical to vec.
    vec = np.dot(rotmatrix,vec0)

    #ax1.hist(vec[0])
    #ax2.hist(vec[1])
    #ax3.hist(vec[2])


    # now obtain the standard deviations
    stda = np.std(vec[0],ddof=1)
    stdb = np.std(vec[1],ddof=1)
    stdc = np.std(vec[2],ddof=1)

    return {'x':meanx,'y':meany,'z':meanz,'a':stda,'b':stdb,'c':stdc,'xy':anglexy,'xz':anglexz,'yz':angleyz}#,ax1,ax2,ax3

def fitellipse2d(x,y):

    ##figellipse = pyplot.figure()
    ##axa = figellipse.add_subplot(221)
    ##axb = figellipse.add_subplot(222)
    ##axc = figellipse.add_subplot(223)
    
    meanx = np.average(x)
    meany = np.average(y)
    
    vec0 = np.asarray([x,y])
    ##axa.scatter(vec0[0],vec0[1],c='k',s=3)
    ##axb.scatter(vec0[0],vec0[2],c='k',s=3)
    ##axc.scatter(vec0[1],vec0[2],c='k',s=3)
    

    # 1: determine XY rotation angle
    m = 1
    b = 1
    coeff = np.polyfit(vec0[0],vec0[1],1)
    anglexy = np.arctan2(coeff[0],1)
    cosc = np.cos(anglexy)
    sinc = np.sin(anglexy)
    rot = np.array([ [cosc,sinc],
                        [-sinc,cosc] ])

    vec = np.dot(rot,vec0)
    stda = np.std(vec[0],ddof=1)
    stdb = np.std(vec[1],ddof=1) 

    return {'x':meanx,'y':meany,'a':stda,'b':stdb,'xy':anglexy}

# code adapted from http://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib
def triaxial(object):

    u = np.linspace(0,2*np.pi,200)
    v = np.linspace(0,np.pi,200)
    
    # generate the triaxial ellipsoid
    plota = object['a'] * np.outer(np.cos(u),np.sin(v))
    plotb = object['b'] * np.outer(np.sin(u),np.sin(v))
    plotc = object['c'] * np.outer(np.ones_like(u),np.cos(v))
    
    #now rotate and translate
    # yz rotates around the x axis, xz rotates around the y axis, and xy rotates around the z axis
    # we have figured out the rotation, so now we need to UNrotate
    cosa = np.cos(-1*object['yz'])
    sina = np.sin(-1*object['yz'])
    cosb = np.cos(-1*object['xz']) # sazzle frazzle
    sinb = np.sin(-1*object['xz'])
    cosc = np.cos(-1*object['xy'])
    sinc = np.sin(-1*object['xy'])
    
    rotz_xy = np.array([ [cosc,sinc,0],
                            [-sinc,cosc,0],
                            [0,0,1] ])
    roty_xz = np.array([ [cosb,0, sinb],
                            [0,1,0],
                            [-sinb,0, cosb] ])
    rotx_yz = np.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])    

    rotmatrix = np.dot(rotx_yz,np.dot(roty_xz,rotz_xy))
    
    vec = np.asarray([plota.reshape([-1]),plotb.reshape([-1]),plotc.reshape([-1])])
    #   print rotmatrix.shape,vec.shape
    triaxial = np.dot(rotmatrix,vec)
    
    triaxialx = triaxial[0].reshape([200,-1]) + object['x']
    triaxialy = triaxial[1].reshape([200,-1]) + object['y']
    triaxialz = triaxial[2].reshape([200,-1]) + object['z']
    
    #print triaxialx
    
    return triaxialx,triaxialy,triaxialz

def rotate(xy,xz,yz):

    # This routine produces a rotation matrix from the angles determined by fitellipse
    cosa = np.cos(yz)
    sina = np.sin(yz)
    cosb = np.cos(xz) # sazzle frazzle
    sinb = np.sin(xz)
    cosc = np.cos(xy)
    sinc = np.sin(xy)

    rotx_yz = np.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])
   
    roty_xz = np.array([ [cosb, 0, sinb],
                            [0,1,0],
                            [-sinb, 0, cosb]])

    rotz_xy = np.array([ [cosc, sinc, 0],
                            [-sinc, cosc, 0], 
                            [0, 0, 1] ])

    # The transpose of this matrix can be used to reverse the rotation, 
    # for instance, to generate points within the ellipse.
    rotmatrix = np.dot(rotx_yz,np.dot(roty_xz,rotz_xy))

    return rotmatrix


def fitellipse_gagne(x,y,z):

    meanx = np.average(x)
    meany = np.average(y)
    meanz = np.average(z)
    
    vec0 = np.asarray([x,y,z])

    # 3: determine YZ
    
    m = 1
    b = 1
    coeff = np.polyfit(vec0[1],vec0[2],1)
    
    angleyz = np.arctan2(coeff[0],1)
    cosa = np.cos(angleyz)
    sina = np.sin(angleyz)
    rotx_yz = np.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])
    
    vec1 = np.dot(rotx_yz,vec0)
    
    # 2: determine XZ
    m = 1
    b = 1
    coeff = np.polyfit(vec1[0],vec1[2],1)
    
    anglexz = np.arctan2(coeff[0],1)
    
    cosb = np.cos(anglexz) # sazzle frazzle
    sinb = np.sin(anglexz)
    roty_xz = np.array([ [cosb,0, -sinb],
                            [0,1,0],
                            [sinb,0, cosb] ])
    
    vec2 = np.dot(roty_xz,vec1)
    
    # 1: determine XY rotation angle
    m = 1
    b = 1
    coeff = np.polyfit(vec2[0],vec2[1],1)
    anglexy = np.arctan2(coeff[0],1)
    
    cosc = np.cos(anglexy)
    sinc = np.sin(anglexy)
    rotz_xy = np.array([ [cosc,sinc,0],
                            [-sinc,cosc,0],
                            [0,0,1] ])
    
    vec3 = np.dot(rotz_xy,vec2)
    
    #now rotate and translate
    # yz rotates around the x axis, xz rotates around the y axis, and xy rotates around the z axis

    rotmatrix = np.dot(rotz_xy,np.dot(roty_xz,rotx_yz))

    # sanity check: If this was done correctly, vec3 should be identical to vec.
    vec = np.dot(rotmatrix,vec0)

    # now obtain the standard deviations
    stda = np.std(vec[0],ddof=1)
    stdb = np.std(vec[1],ddof=1)
    stdc = np.std(vec[2],ddof=1)

    return {'x':meanx,'y':meany,'z':meanz,'a':stda,'b':stdb,'c':stdc,'xy':anglexy,'xz':anglexz,'yz':angleyz}

def rotate_gagne(xy,xz,yz):

    # This routine produces a rotation matrix from the angles determined by fitellipse
    cosa = np.cos(yz)
    sina = np.sin(yz)
    cosb = np.cos(xz) # sazzle frazzle
    sinb = np.sin(xz)
    cosc = np.cos(xy)
    sinc = np.sin(xy)

    rotx_yz = np.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])
   
    roty_xz = np.array([ [cosb, 0, -sinb],
                            [0,1,0],
                            [sinb, 0, cosb]])

    rotz_xy = np.array([ [cosc, sinc, 0],
                            [-sinc, cosc, 0], 
                            [0, 0, 1] ])

    # The transpose of this matrix can be used to reverse the rotation, 
    # for instance, to generate points within the ellipse.
    rotmatrix = np.dot(rotz_xy,np.dot(roty_xz,rotx_yz))

    return rotmatrix

def rotationtester():
    # This tests whether both versions of the rotation matrix do the same thing
    tra = np.random.rand()*360 + np.random.randn(10000) * 0.00004
    tdec = np.random.rand()*180-90 + np.random.randn(10000) * 0.00003
    tdist = np.random.rand()*300 + np.random.randn(10000)* 1.1
    tpmra =  np.random.rand()*2-1 + np.random.randn(10000) * 0.00320
    tpmdec = np.random.rand()*2-1 + np.random.randn(10000) * 0.00320
    trv = np.random.rand()*100-50 + np.random.randn(10000) * 0.74

    u,v,w,x,y,z = kinematics.gal_uvwxyz(ra=tra,dec=tdec,distance=tdist,pmra=tpmra,pmdec=tpmdec,vrad=trv)
    
    print "Regular"
    obj = fitellipse(x,y,z)
    rotmatrix = rotate(obj['xy'],obj['xz'],obj['yz'])
    print rotmatrix

    print "Gagne"
    obj2 = fitellipse_gagne(x,y,z)
    rotmatrix2 = rotate_gagne(obj2['xy'],obj2['xz'],obj2['yz'])
    print rotmatrix2
    
    # The rotation matrices are not the same, but if you rotate with one
    #  and derotate with the other, it's zero.

    rotated = np.dot(rotmatrix,[x,y,z])
    
    
    derotate = np.dot(rotmatrix2.transpose(),rotated)

    print derotate
    print rotated
    print [x,y,z]
    print [x,y,z]-derotate
