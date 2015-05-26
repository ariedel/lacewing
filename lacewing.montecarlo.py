pimport numpy
from numpy import linalg
import scipy
from astropy.io import ascii
from sys import argv
import time

import kinematics
import converge
import ellipse

###################################
###################################
### MAIN ROUTINE
###################################
###################################

class Mgp:
    def __init__(self,mgp):
        self.name = mgp["Name"]
        self.U = numpy.float(mgp["U"])
        self.A = numpy.float(mgp["A"])
        self.V = numpy.float(mgp["V"])
        self.B = numpy.float(mgp["B"])
        self.W = numpy.float(mgp["W"])
        self.C = numpy.float(mgp["C"])
        self.UV = numpy.float(mgp["UV"])
        self.UW = numpy.float(mgp["UW"])
        self.VW = numpy.float(mgp["VW"])
        self.X = numpy.float(mgp["X"])
        self.D = numpy.float(mgp["D"])
        self.Y = numpy.float(mgp["Y"])
        self.E = numpy.float(mgp["E"])
        self.Z = numpy.float(mgp["Z"])
        self.F = numpy.float(mgp["F"])
        self.XY = numpy.float(mgp["XY"])
        self.XZ = numpy.float(mgp["XZ"])
        self.YZ = numpy.float(mgp["YZ"])

        self.color = [mgp["Red"],mgp["Green"],mgp["Blue"]]

        self.coeff_all = [ [mgp['field_all_A'],mgp['field_all_S']],[mgp['field_pm_A'],mgp['field_pm_S']],[mgp['field_dist_A'],mgp['field_dist_S']],[mgp['field_rv_A'],mgp['field_rv_S']],[mgp['field_pmdist_A'],mgp['field_pmdist_S']],[mgp['field_pmrv_A'],mgp['field_pmrv_S']],[mgp['field_distrv_A'],mgp['field_distrv_S']]]
        self.coeff_young = [ [mgp['young_all_A'],mgp['young_all_S']],[mgp['young_pm_A'],mgp['young_pm_S']],[mgp['young_dist_A'],mgp['young_dist_S']],[mgp['young_rv_A'],mgp['young_rv_S']],[mgp['young_pmdist_A'],mgp['young_pmdist_S']],[mgp['young_pmrv_A'],mgp['young_pmrv_S']],[mgp['young_distrv_A'],mgp['young_distrv_S']]]
        for x in range(len(self.coeff_all)):
            if kinematics.isnumber(self.coeff_all[x][0]):
                self.coeff_all[x][0] = numpy.float(self.coeff_all[x][0])
            if kinematics.isnumber(self.coeff_all[x][1]):
                self.coeff_all[x][1] = numpy.float(self.coeff_all[x][1])
        for x in range(len(self.coeff_young)):
            if kinematics.isnumber(self.coeff_young[x][0]):
                self.coeff_young[x][0] = numpy.float(self.coeff_young[x][0])
            if kinematics.isnumber(self.coeff_young[x][1]):
                self.coeff_young[x][1] = numpy.float(self.coeff_young[x][1])

# read in ellipse parameters for associations

balance = argv[1]

file = open('Moving_Group_all.csv','rb')
readtable = ascii.get_reader(Reader=ascii.Basic)
readtable.header.splitter.delimiter = ','
readtable.data.splitter.delimiter = ','
readtable.header.start_line = 0
readtable.data.start_line = 1
groups = readtable.read(file)
file.close()

moving_groups = []

for i in range(len(groups)):
    moving_groups.append(Mgp(groups[i]))

outfile = []
file = argv[2]

outfile.append(open('epscha{0:}.csv'.format(file),'wb'))
outfile.append(open('etacha{0:}.csv'.format(file),'wb'))
outfile.append(open('twhya{0:}.csv'.format(file),'wb'))
outfile.append(open('betapic{0:}.csv'.format(file),'wb'))
outfile.append(open('octans{0:}.csv'.format(file),'wb'))
#outfile.append(open('chanear{0:}.csv'.format(file),'wb'))
outfile.append(open('tuchor{0:}.csv'.format(file),'wb'))
outfile.append(open('columba{0:}.csv'.format(file),'wb'))
#outfile.append(open('carina{0:}.csv'.format(file),'wb'))
outfile.append(open('argus{0:}.csv'.format(file),'wb'))
#outfile.append(open('carinanear{0:}.csv'.format(file),'wb'))
#outfile.append(open('carinavela{0:}.csv'.format(file),'wb'))
#outfile.append(open('ic2391{0:}.csv'.format(file),'wb'))
#outfile.append(open('b4{0:}.csv'.format(file),'wb'))
outfile.append(open('abdor{0:}.csv'.format(file),'wb'))
outfile.append(open('plei{0:}.csv'.format(file),'wb'))
outfile.append(open('herlyr{0:}.csv'.format(file),'wb'))
#outfile.append(open('castor{0:}.csv'.format(file),'wb'))
outfile.append(open('comaber{0:}.csv'.format(file),'wb'))
outfile.append(open('uma{0:}.csv'.format(file),'wb'))
outfile.append(open('hya{0:}.csv'.format(file),'wb'))
outfile.append(open('field{0:}.csv'.format(file),'wb'))

#for e in range(len(outfile)):
#   outfile[e].write('idx,Name,RA,eRA,DEC,eDEC,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv,mgp,d1,sig_all,sig_padist,sig_parv,sig_distrv,d3,d4,sig_pa,d5,d6,d7,sig_dist,d8,d9,d10,sig_rv,d11,MGP1,Note\n')

# set up iterations loop
iterations = 2000000

lineno = 0
for i in xrange(iterations):
    # Now begin the random number generators
    time20 = time.time()
    # 1. Decide which type of star this is going to be.
    selector = numpy.random.rand()
    # Choose the first type greater than 'selector'
    startype = numpy.where(numpy.asarray(groups['Weightednumber']) > selector)[0][0]
    mgp1 = groups[startype]['Name'].replace(' ','_')
    name = '{0:8d}_{1:}'.format(i,mgp1)
    
    time21 = time.time()

    # 2. Now generate a star within the boundaries of that group
    tu = numpy.random.randn()*float(groups[startype]['A'])
    tv = numpy.random.randn()*float(groups[startype]['B'])
    tw = numpy.random.randn()*float(groups[startype]['C'])
    if groups[startype]['uniform'] == 0: # uniform random distribution of positions
        tx = ((numpy.random.rand()*2)-1)*float(groups[startype]['D'])
        ty = ((numpy.random.rand()*2)-1)*float(groups[startype]['E'])
        tz = ((numpy.random.rand()*2)-1)*float(groups[startype]['F'])
        while ((tx/float(groups[startype]['D']))**2 + (ty/float(groups[startype]['E']))**2 + (tz/float(groups[startype]['F'])**2)) > 1:
            tx = ((numpy.random.rand()*2)-1)*float(groups[startype]['D'])
            ty = ((numpy.random.rand()*2)-1)*float(groups[startype]['E'])
            tz = ((numpy.random.rand()*2)-1)*float(groups[startype]['F'])
        # A quick test shows that an ellipse fit to a uniform distribution
        # has axes smaller than the initial ellipse by a factor of sqrt(3)
        tx = tx * numpy.sqrt(3)
        ty = ty * numpy.sqrt(3)
        tz = tz * numpy.sqrt(3)

    if groups[startype]['uniform'] == 1: # for clusters: use gaussians
        tx = numpy.random.randn()*float(groups[startype]['D'])
        ty = numpy.random.randn()*float(groups[startype]['E'])
        tz = numpy.random.randn()*float(groups[startype]['F'])

    if groups[startype]['uniform'] == 2: # exponential disk dropoff (for the field stars)
        tx = ((numpy.random.rand()*2)-1)*float(groups[startype]['D'])
        ty = ((numpy.random.rand()*2)-1)*float(groups[startype]['E'])
        tz = numpy.random.exponential(scale=300)*float(groups[startype]['F']) # problem: This generates points from 0 to infinity with a scale height of 300.
        while ((tx/float(groups[startype]['D']))**2 + (ty/float(groups[startype]['E']))**2 + (tz/float(groups[startype]['F'])**2)) > 1:
            tx = ((numpy.random.rand()*2)-1)*float(groups[startype]['D'])
            ty = ((numpy.random.rand()*2)-1)*float(groups[startype]['E'])
            tz = numpy.random.exponential(scale=300)*float(groups[startype]['F']) # problem: This generates points from 0 to infinity with a scale height of 300.
        if numpy.random.rand() > 0.5: # mirror half the points below the galactic plane
            tz = -1*tz
         
    time22 = time.time()
  
    # we need to rotate the points back to the correct ellipse
    #now rotate and translate
    # yz rotates around the x axis, xz rotates around the y axis, and xy rotates around the z axis
    # we have figured out the rotation, so now we need to UNrotate
    cosa = numpy.cos(-1*groups[startype]['VW'])
    sina = numpy.sin(-1*groups[startype]['VW'])
    cosb = numpy.cos(-1*groups[startype]['UW']) # sazzle frazzle
    sinb = numpy.sin(-1*groups[startype]['UW'])
    cosc = numpy.cos(-1*groups[startype]['UV'])
    sinc = numpy.sin(-1*groups[startype]['UV'])
    
    rotz_xy = numpy.array([ [cosc,sinc,0],
                           [-sinc,cosc,0],
                           [0,0,1] ])
    roty_xz = numpy.array([ [cosb,0, sinb],
                           [0,1,0],
                           [-sinb,0, cosb] ])
    rotx_yz = numpy.array([ [1, 0, 0],
                           [0, cosa, sina],
                           [0, -sina, cosa] ])
   
    rotmatrix = numpy.dot(rotz_xy,numpy.dot(roty_xz,rotx_yz))
    
    tu,tv,tw = numpy.dot(rotmatrix,[tu,tv,tw])

    cosa = numpy.cos(-1*groups[startype]['YZ'])
    sina = numpy.sin(-1*groups[startype]['YZ'])
    cosb = numpy.cos(-1*groups[startype]['XZ']) # sazzle frazzle
    sinb = numpy.sin(-1*groups[startype]['XZ'])
    cosc = numpy.cos(-1*groups[startype]['XY'])
    sinc = numpy.sin(-1*groups[startype]['XY'])
    
    rotz_xy = numpy.array([ [cosc,sinc,0],
                            [-sinc,cosc,0],
                            [0,0,1] ])
    roty_xz = numpy.array([ [cosb,0, sinb],
                            [0,1,0],
                            [-sinb,0, cosb] ])
    rotx_yz = numpy.array([ [1, 0, 0],
                            [0, cosa, sina],
                            [0, -sina, cosa] ])    
    
    rotmatrix = numpy.dot(rotz_xy,numpy.dot(roty_xz,rotx_yz))
    
    tx,ty,tz = numpy.dot(rotmatrix,[tx,ty,tz])
   
    tu = tu + groups[startype]['U']+numpy.random.randn()*0.02
    tv = tv + groups[startype]['V']+numpy.random.randn()*0.02
    tw = tw + groups[startype]['W']+numpy.random.randn()*0.02
    tx = tx + groups[startype]['X']+numpy.random.randn()*0.02
    ty = ty + groups[startype]['Y']+numpy.random.randn()*0.02
    tz = tz + groups[startype]['Z']+numpy.random.randn()*0.02
    
    # now we have a properly generated star; convert it to observables and finish setting up the "observation"
   
    ra,dec,dist,pmra,pmdec,rv = kinematics.gal_rdp(tu,tv,tw,tx,ty,tz)
    plx = 1/dist
    era = abs(numpy.random.randn()*0.1/3600.)
    edec = abs(numpy.random.randn()*0.1/3600.)
    eplx = abs(numpy.random.randn()*dist/20.)/(dist**2)
    epmra = abs(numpy.random.randn()*0.0005)
    epmdec = abs(numpy.random.randn()*0.0005)
    erv = abs(numpy.random.randn()*1.0)
    #print ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv

    pmexists = 1
    distexists = 1
    rvexists = 1

    # We've got our star set up.  Now run it through the convergence code

    for a in xrange(len(moving_groups)):
       exp_pmra,exp_epmra,exp_pmdec,exp_epmdec,exp_rv,exp_erv = converge.converge(moving_groups[a],ra,era,dec,edec,1000)
       #print moving_groups[a].name,exp_pmra,exp_epmra,exp_pmdec,exp_epmdec,exp_rv,exp_erv
       mgpname = moving_groups[a].name.replace(' ','_')
       
       exp_pm,exp_epm,exp_pa,exp_epa = astrometry.pmjoin(exp_pmra,exp_epmra,exp_pmdec,exp_epmdec)
       pm_new = 0
       epm_new = 0
       cosa = 0
       sina = 0
       exp_dist = 0
       exp_edist = 0
       exp_pm_new = 0
       exp_epm_new = 0
       pm_sig = 0
       rv_sig = 0
       dist_sig = 0
       pos_sig = 0
       sig = 0

       # There are now 7 cases to test for:
       # 0: Nothing. We can't do anything without input.
       # 1. Only PM: the most basic case, check the agreement between the expected and actual vector. We can also run the spatial test.
       # 2. Only DIST: We can run the spatial test, but that's basically it.
       # 3. Only RV: we DO predict an RV, so this case WILL let us do something. This is the only case that won't allow the spatial test.
       # 4. PM + DIST: We can predict both the proper motion, and the proper motion with a distance, and do a spatial test.
       # 5. PM + RV : Now with two options.
       # 6. DIST + RV: RV test with spatial test.
       # 7. PM + DIST + RV: All four measures come into play here.

       # There are only four cases that really matter: PM, PM+DIST, RV, and the spatial test.
       sig_pm = []
       weight_pm=[]
       sig_dist = []
       weight_dist=[]
       sig_rv=[]
       weight_rv=[]
       sig_pmdist = []
       weight_pmdist = []
       sig_pmrv = []
       weight_pmrv = []
       sig_distrv = []
       weight_distrv = []
       sig_all = []
       weight_all = []

       pm_string = 'PM:,s=,"", '
       dist_string = 'Dist:,s=,"", '
       kinemag_string = 'Kinemag:,s=,"",  '
       absmag_string = 'Absmag:,s=,"", '
       rv_string = 'RV:,s=,"", '
       pos_string = 'XYZ: s=,"",'

       # step 1: compare proper motions. We're going to rotate the star's proper motion 
       #  vectors (pmRA, pmDEC) by the angle of the association's proper motion, to 
       #  transform the motion into pm= and pm_|_.
       if pmexists == 1:
           # 1a. Calculate the parallel and perpendicular components of the proper motion
           #   First: compute a rotation matrix to transform pmRA,pmDEC into pm= and pm_|_
           cosa = numpy.cos(exp_pa*numpy.pi/180.0)
           sina = numpy.sin(exp_pa*numpy.pi/180.0)
           rotmatrix = [ [ cosa, -sina ] , [sina, cosa] ]
           
           #   Second: Apply this rotation matrix to the measured and expected proper motion vectors.
           pm_new = numpy.dot(rotmatrix,[pmra,pmdec])
           epm_new = numpy.dot(rotmatrix,[epmra,epmdec])
           exp_pm_new = numpy.dot(rotmatrix,[exp_pmra,exp_pmdec]) # sanity check! Should be [0,pm] if correctly done
           exp_epm_new = numpy.dot(rotmatrix,[exp_epmra,exp_epmdec])
           # The 0th element is the perpendicular component, the 1st is parallel

           # 1b. Compute a kinematic distance by comparing the magnitude of the expected pm vector
           #   (which was calculated for a distance of 10 parsecs) to the measured parallel one
           exp_dist = exp_pm_new[1]/pm_new[1] * 10.0
           exp_edist = numpy.sqrt((exp_epm_new[1]/exp_pm_new[1])**2 + (epm_new[1]/pm_new[1])**2) * exp_dist

           # 1c. Now we use that kinematic distance to scale the error out. Why? If the proper motion was
           #  in completely the wrong direction, but tiny, it would come out as a match when it shouldn't. 
           #  This scales the error appropriately to stop that from happening.
           exp_epm_new = exp_epm_new * (10/exp_edist)
            
           # 1d. Compute the proper motion sigma based on the perpendicular component (with scaled errors)
           pm_sig = numpy.sqrt(pm_new[0]**2)/(numpy.sqrt(epm_new[0]**2 + exp_epm_new[0]**2))

           # 1e. this should weed out anything where the pm vector is substantial and flipped (and
           #   distance was negative). I initially thought this should go under the
           #   "distance" section, but I can verify it entirely based on proper
           #   motions even when a distance isn't known, so it should go here.
           if ((abs(exp_epm_new[1]) < abs(exp_pm_new[1])) & (exp_dist < 0)):
               pm_sig = pm_sig + 1000.
               # I might also try this relative to a km/sec velocity dispersion, pm_new[0]*exp_dist*4.74

           # 1f. Now normalize the sigma and add the weight to the array of weights (if desired)
           # 2014.0202: PA_SIG values tend to be very, very small for good results.  Normalize them appropriately (on a per-association basis)
           if balance == 'weighted':
               pm_sig = pm_sig / pm_norm
               weight_all.append(pm_weight)
               weight_pmdist.append(pm_weight)
               weight_pmrv.append(pm_weight)
               
           sig_all.append(pm_sig)
           sig_pmdist.append(pm_sig)
           sig_pm.append(pm_sig)
           sig_pmrv.append(pm_sig)
           pm_string = 'PM:,s=,{0: 8.3f},  '.format(pm_sig)

           pospm_sig,pos_esig = converge.spatial(moving_groups[a],ra,era,dec,edec,1/exp_dist,exp_edist/(exp_dist)**2,1000)

           # 1g. normalize the sigma and add the weight to the array of weights (if desired)
           if balance == 'weighted':
               pospm_sig = pospm_sig / pospm_norm
               weight_pm.append(pospm_weight)
               weight_pmrv.append(pospm_weight)

           sig_pm.append(pospm_sig)
           sig_pmrv.append(pospm_sig)
           pospm_string = 'XYZpm: s=,{0: 8.3f},'.format(pospm_sig)


       # step 2: Distances exist. Without a proper motion, we can run the spatial test but nothing else.
       if distexists == 1:
           dist_string = 'Dist:,s=,"",  '
           posdist_sig,pos_esig = converge.spatial(moving_groups[a],ra,era,dec,edec,plx,eplx,1000)
           if balance == 'weighted':
               posdist_sig = posdist_sig / pos_norm
               weight_dist.append(posdist_weight)
               weight_pmdist.append(posdist_weight)
               weight_distrv.append(posdist_weight)
               weight_all.append(posdist_weight)

           sig_dist.append(posdist_sig)
           sig_pmdist.append(posdist_sig)
           sig_distrv.append(posdist_sig)
           sig_all.append(posdist_sig)
           posdist_string = 'XYZdist: s=,{0: 8.3f},'.format(posdist_sig)


       # step 3: Distances AND proper motions exist. Now we can compute a kinematic distance error.
       if (pmexists == 1) & (distexists == 1):

           # 3a. use the kinematic distance calculated in step 1c, which should have already run if pmexists==1.
           dist_sig = numpy.abs(1/plx - exp_dist) / numpy.sqrt((eplx/plx**2)**2 + exp_edist**2)

           # 3b. normalize the sigma and add the weight to the array of weights (if desired)
           # 2014.0202: DIST_SIG values aren't quite as small, but still need to be normalized.
           if balance == 'weighted':
               dist_sig = dist_sig / dist_norm
               weight_all.append(dist_weight)
               weight_pmdist.append(dist_weight)
               weight_distrv.append(dist_weight)
         
           sig_all.append(dist_sig)
           sig_pmdist.append(dist_sig)
           sig_distrv.append(dist_sig)
           dist_string = 'Dist:,s=,{0: 8.3f}, '.format(dist_sig)



       # step 4: RVs exist. This is straightforward: We have an estimated RV from the convergence.
       if (rvexists == 1):
           rv_sig = numpy.abs(rv - exp_rv) / numpy.sqrt(erv**2 + exp_erv**2)

           # 4b. normalize the sigma and add the weight to the array of weights (if desired)
           if balance == 'weighted':
               rv_sig = rv_sig / rv_norm
               weight_rv.append(rv_weight)
               weight_pmrv.append(rv_weight)
               weight_distrv.append(rv_weight)
               weight_all.append(rv_weight)

           sig_all.append(rv_sig)
           sig_rv.append(rv_sig)
           sig_pmrv.append(rv_sig)
           sig_distrv.append(rv_sig)
           rv_string = 'RV:,s=,{0: 8.3f}, '.format(rv_sig)


       sigma_all = converge.weightadd(sig_all)
       sigma_pm = converge.weightadd(sig_pm)
       sigma_dist = converge.weightadd(sig_dist)
       sigma_rv = converge.weightadd(sig_rv)
       sigma_pmdist = converge.weightadd(sig_pmdist)
       sigma_pmrv = converge.weightadd(sig_pmrv)
       sigma_distrv = converge.weightadd(sig_distrv)

       #print '{14:},"{0:16}",{1:09.5f},{2:+08.5f},{3:+6.4f},{4:+6.4f},I,{5:+6.2f},{6:12},SIG=,{7: 8.2f},{8:},{9:},{10:},{11:},{12:},{13:},.'.format(name,ra,dec,pmra,pmdec,imag,assocname,sig_all,sig_pm,sig_rv,sig_pmdist,sig_pmrv,sig_distrv,pm_string,dist_string,rv_string,mgp1,lineno)
       time210 = time.time()      
       outfile[a].write('{0:},"{1:16}",{2:09.5f},{3:09.6f},{4:+08.5f},{5:09.6f},{6:09.5f},{7:09.5f},{8:+6.4f},{9:6.4f},{10:+6.4f},{11:6.4f},{12:+6.4f},{13:+6.4f},{14:},SIG=,{15: 8.3f},{16: 8.3f},{17: 8.3f},{18: 8.3f},{19:},{20:},{21:},{22:},{23:},{24:}.\n'.format(lineno,name,ra,era,dec,edec,1/plx,eplx/(plx**2),pmra,epmra,pmdec,epmdec,rv,erv,mgpname,sigma_all,sigma_pmdist,sigma_pmrv,sigma_distrv,pm_string,dist_string,rv_string,pospm_string,posdist_string,mgp1))

       lineno = lineno+1


for i in outfile:
    i.close()
