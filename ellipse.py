import numpy
import time
import kinematics

#  Functions in this module:
# fitellipse: Fits a 3D ellipsoid by performing linear fits to successive planes of the 3D array, derotating the ellipsoid each time. 
# rotate: Constructs the rotation matrix from fitellipse using the input Tait-Bryant angles.
# fitellipse2d: A modified version of fitellipse for two dimensions.  Useful if you want to make plots of specific planes (UV,UW,VW, for example) because the pices of an ellipse are not interchangeable.
# triaxial: Takes the dictionary produced by fitellipse, returns an object that can be plotted by mplot3D's surface plotter.


# This function fits a 9 parameter ellipse by making linear fits on an axis-by-axis basis.
# The input values are numpy arrays of points to fit.
# This function is intended to be run in a monte carlo simulation.
def fitellipse(x,y,z):

    meanx = numpy.average(x)
    meany = numpy.average(y)
    meanz = numpy.average(z)    outfile.write('\n')
    
    vec0 = numpy.asarray([x,y,z])

    # 1: determine XY rotation angle
    coeff = numpy.polyfit(vec0[0],vec0[1],1)
    anglexy = numpy.arctan2(coeff[0],1)
    
    cosc = numpy.cos(anglexy)
    sinc = numpy.sin(anglexy)
    rotz_xy = numpy.array([ [cosc,sinc,0],
                            [-sinc,cosc,0],
                            [0,0,1] ])
    
    vec1 = numpy.dot(rotz_xy,vec0)
    
    # 2: determine XZ 
    
    m = 1
    b = 1
    coeff = numpy.polyfit(vec1[0],vec1[2],1)
    
    anglexz = numpy.arctan2(coeff[0],1)
    
    cosb = numpy.cos(anglexz) # sazzle frazzle
    sinb = numpy.sin(anglexz)
    roty_xz = numpy.array([ [cosb,0, sinb],
                            [0,1,0],
                            [-sinb,0, cosb] ])
    
    vec2 = numpy.dot(roty_xz,vec1)
    
    # 3: determine YZ
    
    m = 1
    b = 1
    coeff = numpy.polyfit(vec2[1],vec2[2],1)
    
    angleyz = numpy.arctan2(coeff[0],1)
    cosa = numpy.cos(angleyz)
    sina = numpy.sin(angleyz)
    rotx_yz = numpy.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])
    
    vec3 = numpy.dot(rotx_yz,vec2)
    
    #now rotate and translate
    # yz rotates around the x axis, xz rotates around the y axis, and xy rotates around the z axis

    rotmatrix = numpy.dot(rotx_yz,numpy.dot(roty_xz,rotz_xy))

    # sanity check: If this was done correctly, vec3 should be identical to vec.
    vec = numpy.dot(rotmatrix,vec0)

    # now obtain the standard deviations
    stda = numpy.std(vec[0],ddof=1)
    stdb = numpy.std(vec[1],ddof=1)
    stdc = numpy.std(vec[2],ddof=1)

    return {'x':meanx,'y':meany,'z':meanz,'a':stda,'b':stdb,'c':stdc,'xy':anglexy,'xz':anglexz,'yz':angleyz}

def fitellipse2d(x,y):

    ##figellipse = pyplot.figure()
    ##axa = figellipse.add_subplot(221)
    ##axb = figellipse.add_subplot(222)
    ##axc = figellipse.add_subplot(223)
    
    meanx = numpy.average(x)
    meany = numpy.average(y)
    
    vec0 = numpy.asarray([x,y])
    ##axa.scatter(vec0[0],vec0[1],c='k',s=3)
    ##axb.scatter(vec0[0],vec0[2],c='k',s=3)
    ##axc.scatter(vec0[1],vec0[2],c='k',s=3)
    

    # 1: determine XY rotation angle
    m = 1
    b = 1
    coeff = numpy.polyfit(vec0[0],vec0[1],1)
    anglexy = numpy.arctan2(coeff[0],1)
    cosc = numpy.cos(anglexy)
    sinc = numpy.sin(anglexy)
    rot = numpy.array([ [cosc,sinc],
                        [-sinc,cosc] ])

    vec = numpy.dot(rot,vec0)
    stda = numpy.std(vec[0],ddof=1)
    stdb = numpy.std(vec[1],ddof=1) 

    return {'x':meanx,'y':meany,'a':stda,'b':stdb,'xy':anglexy}

# code adapted from http://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib
def triaxial(object):

    u = numpy.linspace(0,2*numpy.pi,200)
    v = numpy.linspace(0,numpy.pi,200)
    
    # generate the triaxial ellipsoid
    plota = object['a'] * numpy.outer(numpy.cos(u),numpy.sin(v))
    plotb = object['b'] * numpy.outer(numpy.sin(u),numpy.sin(v))
    plotc = object['c'] * numpy.outer(numpy.ones_like(u),numpy.cos(v))
    
    #now rotate and translate
    # yz rotates around the x axis, xz rotates around the y axis, and xy rotates around the z axis
    # we have figured out the rotation, so now we need to UNrotate
    cosa = numpy.cos(-1*object['yz'])
    sina = numpy.sin(-1*object['yz'])
    cosb = numpy.cos(-1*object['xz']) # sazzle frazzle
    sinb = numpy.sin(-1*object['xz'])
    cosc = numpy.cos(-1*object['xy'])
    sinc = numpy.sin(-1*object['xy'])
    
    rotz_xy = numpy.array([ [cosc,sinc,0],
                            [-sinc,cosc,0],
                            [0,0,1] ])
    roty_xz = numpy.array([ [cosb,0, sinb],
                            [0,1,0],
                            [-sinb,0, cosb] ])
    rotx_yz = numpy.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])    

    rotmatrix = numpy.dot(rotx_yz,numpy.dot(roty_xz,rotz_xy))
    
    vec = numpy.asarray([plota.reshape([-1]),plotb.reshape([-1]),plotc.reshape([-1])])
    #   print rotmatrix.shape,vec.shape
    triaxial = numpy.dot(rotmatrix,vec)
    
    triaxialx = triaxial[0].reshape([200,-1]) + object['x']
    triaxialy = triaxial[1].reshape([200,-1]) + object['y']
    triaxialz = triaxial[2].reshape([200,-1]) + object['z']
    
    #print triaxialx
    
    return triaxialx,triaxialy,triaxialz

def rotate(xy,xz,yz):

    # This routine produces a rotation matrix from the angles determined by fitellipse
    cosa = numpy.cos(yz)
    sina = numpy.sin(yz)
    cosb = numpy.cos(xz) # sazzle frazzle
    sinb = numpy.sin(xz)
    cosc = numpy.cos(xy)
    sinc = numpy.sin(xy)

    rotx_yz = numpy.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])
   
    roty_xz = numpy.array([ [cosb, 0, sinb],
                            [0,1,0],
                            [-sinb, 0, cosb]])

    rotz_xy = numpy.array([ [cosc, sinc, 0],
                            [-sinc, cosc, 0], 
                            [0, 0, 1] ])

    # The transpose of this matrix can be used to reverse the rotation, 
    # for instance, to generate points within the ellipse.
    rotmatrix = numpy.dot(rotx_yz,numpy.dot(roty_xz,rotz_xy))

    return rotmatrix


def fitellipse_gagne(x,y,z):

    meanx = numpy.average(x)
    meany = numpy.average(y)
    meanz = numpy.average(z)
    
    vec0 = numpy.asarray([x,y,z])

    # 3: determine YZ
    
    m = 1
    b = 1
    coeff = numpy.polyfit(vec0[1],vec0[2],1)
    
    angleyz = numpy.arctan2(coeff[0],1)
    cosa = numpy.cos(angleyz)
    sina = numpy.sin(angleyz)
    rotx_yz = numpy.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])
    
    vec1 = numpy.dot(rotx_yz,vec0)
    
    # 2: determine XZ
    m = 1
    b = 1
    coeff = numpy.polyfit(vec1[0],vec1[2],1)
    
    anglexz = numpy.arctan2(coeff[0],1)
    
    cosb = numpy.cos(anglexz) # sazzle frazzle
    sinb = numpy.sin(anglexz)
    roty_xz = numpy.array([ [cosb,0, -sinb],
                            [0,1,0],
                            [sinb,0, cosb] ])
    
    vec2 = numpy.dot(roty_xz,vec1)
    
    # 1: determine XY rotation angle
    m = 1
    b = 1
    coeff = numpy.polyfit(vec2[0],vec2[1],1)
    anglexy = numpy.arctan2(coeff[0],1)
    
    cosc = numpy.cos(anglexy)
    sinc = numpy.sin(anglexy)
    rotz_xy = numpy.array([ [cosc,sinc,0],
                            [-sinc,cosc,0],
                            [0,0,1] ])
    
    vec3 = numpy.dot(rotz_xy,vec2)
    
    #now rotate and translate
    # yz rotates around the x axis, xz rotates around the y axis, and xy rotates around the z axis

    rotmatrix = numpy.dot(rotz_xy,numpy.dot(roty_xz,rotx_yz))

    # sanity check: If this was done correctly, vec3 should be identical to vec.
    vec = numpy.dot(rotmatrix,vec0)

    # now obtain the standard deviations
    stda = numpy.std(vec[0],ddof=1)
    stdb = numpy.std(vec[1],ddof=1)
    stdc = numpy.std(vec[2],ddof=1)

    return {'x':meanx,'y':meany,'z':meanz,'a':stda,'b':stdb,'c':stdc,'xy':anglexy,'xz':anglexz,'yz':angleyz}

def rotate_gagne(xy,xz,yz):

    # This routine produces a rotation matrix from the angles determined by fitellipse
    cosa = numpy.cos(yz)
    sina = numpy.sin(yz)
    cosb = numpy.cos(xz) # sazzle frazzle
    sinb = numpy.sin(xz)
    cosc = numpy.cos(xy)
    sinc = numpy.sin(xy)

    rotx_yz = numpy.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])
   
    roty_xz = numpy.array([ [cosb, 0, -sinb],
                            [0,1,0],
                            [sinb, 0, cosb]])

    rotz_xy = numpy.array([ [cosc, sinc, 0],
                            [-sinc, cosc, 0], 
                            [0, 0, 1] ])

    # The transpose of this matrix can be used to reverse the rotation, 
    # for instance, to generate points within the ellipse.
    rotmatrix = numpy.dot(rotz_xy,numpy.dot(roty_xz,rotx_yz))

    return rotmatrix

def rotationtester():
    # This tests whether both versions of the rotation matrix do the same thing
    tra = numpy.random.rand()*360 + numpy.random.randn(10000) * 0.00004
    tdec = numpy.random.rand()*180-90 + numpy.random.randn(10000) * 0.00003
    tdist = numpy.random.rand()*300 + numpy.random.randn(10000)* 1.1
    tpmra =  numpy.random.rand()*2-1 + numpy.random.randn(10000) * 0.00320
    tpmdec = numpy.random.rand()*2-1 + numpy.random.randn(10000) * 0.00320
    trv = numpy.random.rand()*100-50 + numpy.random.randn(10000) * 0.74

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

    rotated = numpy.dot(rotmatrix,[x,y,z])
    
    
    derotate = numpy.dot(rotmatrix2.transpose(),rotated)

    print derotate
    print rotated
    print [x,y,z]
    print [x,y,z]-derotate
