import numpy as np
import scipy
from astropy.io import ascii
from sys import argv
import time

import kinematics
import lacewing
import ellipse

#########################################################
#########################################################
### MAIN ROUTINE
### ARR 2015-10-16
### 1.3: Now compatible with LACEwING v1.3, and uses the 
###      lacewing.lacewing() and lacewing_moving_group_loader()
###      functions
#########################################################
#########################################################

iterations = int(argv[1])
number = argv[2]

moving_groups = lacewing.moving_group_loader()

outfile = []
for i in xrange(len(moving_groups)):
    outfile.append(open('{0:}{1:}'.format(moving_groups[i].name.replace(' ','_'),number),'wb'))

weightednumber = []
for j in xrange(len(moving_groups)):
    weightednumber.append(moving_groups[j].weightednumber)

for i in xrange(iterations):
    # Now begin the random number generators
    time20 = time.time()
    # 1. Decide which type of star this is going to be.
    selector = np.random.rand()
    # Choose the first type greater than 'selector'
    startype = np.where(np.asarray(weightednumber) > selector)[0][0]
    mgp = moving_groups[startype].name.replace(' ','_')
    
    time21 = time.time()

    # 2. Now generate a star within the dispersion of the group.
    tu = np.random.randn()*float(moving_groups[startype].A)
    tv = np.random.randn()*float(moving_groups[startype].B)
    tw = np.random.randn()*float(moving_groups[startype].C)
    if moving_groups[startype].uniform == 0: # uniform random distribution of positions
        tx = ((np.random.rand()*2)-1)*float(moving_groups[startype].D)
        ty = ((np.random.rand()*2)-1)*float(moving_groups[startype].E)
        tz = ((np.random.rand()*2)-1)*float(moving_groups[startype].F)
        while ((tx/float(moving_groups[startype].D))**2 + (ty/float(moving_groups[startype].E))**2 + (tz/float(moving_groups[startype].F)**2)) > 1:
            tx = ((np.random.rand()*2)-1)*float(moving_groups[startype].D)
            ty = ((np.random.rand()*2)-1)*float(moving_groups[startype].E)
            tz = ((np.random.rand()*2)-1)*float(moving_groups[startype].F)
        # A quick test shows that an ellipse fit to a uniform distribution
        # has axes smaller than the initial ellipse by a factor of sqrt(3)
        tx = tx * np.sqrt(3)
        ty = ty * np.sqrt(3)
        tz = tz * np.sqrt(3)

    if moving_groups[startype].uniform == 1: # for clusters: use gaussians
        tx = np.random.randn()*float(moving_groups[startype].D)
        ty = np.random.randn()*float(moving_groups[startype].E)
        tz = np.random.randn()*float(moving_groups[startype].F)

    if moving_groups[startype].uniform == 2: # exponential disk dropoff (for the field stars)
        tx = ((np.random.rand()*2)-1)*float(moving_groups[startype].D)
        ty = ((np.random.rand()*2)-1)*float(moving_groups[startype].E)
        tz = np.random.exponential(scale=300)*float(moving_groups[startype].F) # problem: This generates points from 0 to infinity with a scale height of 300.
        while ((tx/float(moving_groups[startype].D))**2 + (ty/float(moving_groups[startype].E))**2 + (tz/float(moving_groups[startype].F)**2)) > 1:
            tx = ((np.random.rand()*2)-1)*float(moving_groups[startype].D)
            ty = ((np.random.rand()*2)-1)*float(moving_groups[startype].E)
            tz = np.random.exponential(scale=300)*float(moving_groups[startype].F) # problem: This generates points from 0 to infinity with a scale height of 300.
        if np.random.rand() > 0.5: # mirror half the points below the galactic plane
            tz = -1*tz
         
    time22 = time.time()
  
    # we need to rotate the points back to the correct ellipse
    rotmatrix = ellipse.rotate(moving_groups[startype].UV,moving_groups[startype].UW,moving_groups[startype].VW)
    
    tu,tv,tw = np.dot(rotmatrix.transpose(),[tu,tv,tw])

    rotmatrix = ellipse.rotate(moving_groups[startype].XY,moving_groups[startype].XZ,moving_groups[startype].YZ)
    
    tx,ty,tz = np.dot(rotmatrix.transpose(),[tx,ty,tz])
   
    # Now we need to translate these distributed points to where they belong in space.
    tu = tu + moving_groups[startype].U
    tv = tv + moving_groups[startype].V
    tw = tw + moving_groups[startype].W
    tx = tx + moving_groups[startype].X
    ty = ty + moving_groups[startype].Y
    tz = tz + moving_groups[startype].Z
    
    # now we have a properly generated star; convert it to observables and finish setting up the "observation"
   
    ra,dec,dist,pmra,pmdec,rv = kinematics.gal_rdp(tu,tv,tw,tx,ty,tz)
    plx = 1/dist
    era = abs(np.random.randn()*0.1/3600.)
    edec = abs(np.random.randn()*0.1/3600.)
    eplx = abs(np.random.randn()*0.0005)
    epmra = abs(np.random.randn()*0.01)
    epmdec = abs(np.random.randn()*0.01)
    erv = abs(np.random.randn()*1.0)

    # real observations will not be exactly the right value, with an error. Let's add in observational errors to the values.

    ra = ra + np.random.randn()*era/np.cos(dec)
    dec = dec + np.random.randn()*edec
    plx = plx + np.random.randn()*eplx
    pmra = pmra + np.random.randn()*epmra
    pmdec = pmdec + np.random.randn()*epmdec

    # We've got our star set up.  Now run it through the convergence code
    out = lacewing.lacewing(moving_groups,iterate=1000,ra=ra,era=era,dec=dec,edec=edec,pmra=pmra,epmra=epmra,pmdec=pmdec,epmdec=epmdec,plx=plx,eplx=eplx,rv=rv,erv=erv)
    for k in xrange(len(out)):
        # these are the possibilities
        sig_all = [out[k]['pmsig'],out[k]['distsig'],out[k]['rvsig'],out[k]['possig']]
        sig_pm = [out[k]['pmsig'],out[k]['posksig']]
        sig_dist = [out[k]['possig']]
        sig_rv = [out[k]['rvsig']]
        sig_pmdist = [out[k]['pmsig'],out[k]['distsig'],out[k]['possig']]
        sig_pmrv = [out[k]['pmsig'],out[k]['rvsig'],out[k]['posksig']]
        sig_distrv = [out[k]['rvsig'],out[k]['possig']]
        
        sigma_all = lacewing.weightadd(sig_all)
        sigma_pm = lacewing.weightadd(sig_pm)
        sigma_dist = lacewing.weightadd(sig_dist)
        sigma_rv = lacewing.weightadd(sig_rv)
        sigma_pmdist = lacewing.weightadd(sig_pmdist)
        sigma_pmrv = lacewing.weightadd(sig_pmrv)
        sigma_distrv = lacewing.weightadd(sig_distrv)

        time210 = time.time()      
        outfile[k].write('{0:},{1:8.3f},{2:8.3f},{3:8.3f},{4:8.3f},{5:8.3f},{6:8.3f},{7:8.3f},{8:},\n'.format(out[k]['group'],sigma_all,sigma_pm,sigma_dist,sigma_rv,sigma_pmdist,sigma_pmrv,sigma_distrv,mgp))

for i in outfile:
    i.close()
