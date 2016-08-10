import numpy as np
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u 
from sys import argv
import kinematics
import astrometry
from scipy.special import erfc
import ellipse

###############################################
###############################################
### MAIN ROUTINE
###############################################
###############################################

###############################################
### LACEwING 1.5
### ARR 2016-07-05
### 1.0: Everything
### 1.1: Removed dependence on GROUP columns in input data table
###      Fixed dependency on astrometry module
###      Default output should be a summary mode now
### 1.2: Fixed input format
###      All astrometry routines now point to astrometry module
###      Number of groups loaded is now controlled by input .csv file
### 1.3: The LACEwING algorithm is now a separate function from the
###      moving group loader and the csv file reader function.
### 1.4: CSV loader now optionally understands the AAS-journal-style
###      headers in the catalog
### 1.5: Now handles the case of multiple "field" populations correctly
###      and also different numbers/orders of groups
###############################################

######################################################
## Load the moving groups into a moving group class ##
######################################################

def moving_group_loader():

    class Mgp:
        def __init__(self,mgp):
            self.name = mgp["Name"]
            self.U = np.float(mgp["U"])
            self.A = np.float(mgp["A"])
            self.V = np.float(mgp["V"])
            self.B = np.float(mgp["B"])
            self.W = np.float(mgp["W"])
            self.C = np.float(mgp["C"])
            self.UV = np.float(mgp["UV"])
            self.UW = np.float(mgp["UW"])
            self.VW = np.float(mgp["VW"])
            self.X = np.float(mgp["X"])
            self.D = np.float(mgp["D"])
            self.Y = np.float(mgp["Y"])
            self.E = np.float(mgp["E"])
            self.Z = np.float(mgp["Z"])
            self.F = np.float(mgp["F"])
            self.XY = np.float(mgp["XY"])
            self.XZ = np.float(mgp["XZ"])
            self.YZ = np.float(mgp["YZ"])

            self.A2 = np.float(mgp["A2"])
            self.B2 = np.float(mgp["B2"])
            self.C2 = np.float(mgp["C2"])
            self.UV2 = np.float(mgp["UV2"])
            self.UW2 = np.float(mgp["UW2"])
            self.VW2 = np.float(mgp["VW2"])
            self.D2 = np.float(mgp["D2"])
            self.E2 = np.float(mgp["E2"])
            self.F2 = np.float(mgp["F2"])
            self.XY2 = np.float(mgp["XY2"])
            self.XZ2 = np.float(mgp["XZ2"])
            self.YZ2 = np.float(mgp["YZ2"])

            self.color = [mgp["Red"],mgp["Green"],mgp["Blue"]]
            self.weightednumber = np.float(mgp["Weightednumber"])
            self.uniform = np.int(mgp["uniform"])

            self.coeff_all = [ [mgp['field_all_M'],mgp['field_all_S']],[mgp['field_pm_M'],mgp['field_pm_S']],[mgp['field_dist_M'],mgp['field_dist_S']],[mgp['field_rv_M'],mgp['field_rv_S']],[mgp['field_pmdist_M'],mgp['field_pmdist_S']],[mgp['field_pmrv_M'],mgp['field_pmrv_S']],[mgp['field_distrv_M'],mgp['field_distrv_S']]]
            self.coeff_young = [ [mgp['young_all_M'],mgp['young_all_S']],[mgp['young_pm_M'],mgp['young_pm_S']],[mgp['young_dist_M'],mgp['young_dist_S']],[mgp['young_rv_M'],mgp['young_rv_S']],[mgp['young_pmdist_M'],mgp['young_pmdist_S']],[mgp['young_pmrv_M'],mgp['young_pmrv_S']],[mgp['young_distrv_M'],mgp['young_distrv_S']]]
            for x in range(len(self.coeff_all)):
                if astrometry.isnumber(self.coeff_all[x][0]):
                    self.coeff_all[x][0] = np.float(self.coeff_all[x][0])
                if astrometry.isnumber(self.coeff_all[x][1]):
                    self.coeff_all[x][1] = np.float(self.coeff_all[x][1])
            for x in range(len(self.coeff_young)):
                if astrometry.isnumber(self.coeff_young[x][0]):
                    self.coeff_young[x][0] = np.float(self.coeff_young[x][0])
                if astrometry.isnumber(self.coeff_young[x][1]):
                    self.coeff_young[x][1] = np.float(self.coeff_young[x][1])

    # read in ellipse parameters for associations

    file = open('Moving_Group_all.csv','rb')
    readtable = ascii.get_reader(Reader=ascii.Basic)
    readtable.header.splitter.delimiter = ','
    readtable.data.splitter.delimiter = ','
    readtable.header.start_line = 0
    readtable.data.start_line = 1
    groups = readtable.read(file)
    file.close()

    moving_groups = []
    # AR 2015.0326: Number of moving groups is now controlled by the input file. (Last ones should be "field")
    for i in range(len(groups)):
        moving_groups.append(Mgp(groups[i]))

    return moving_groups

# Calculate the proper motion vector at the star's RA and DEC and a distance of 10 pc
def converge(assoc,ra,era,dec,edec,n_init):

   ra_tmp = ra + np.random.randn(n_init) * era / np.cos(dec*180/np.pi)
   dec_tmp = dec + np.random.randn(n_init) * edec

   # generate random points within the moving group distribution
   a = assoc.A * np.random.randn(n_init)
   b = assoc.B * np.random.randn(n_init)
   c = assoc.C * np.random.randn(n_init)
   # rotate them to account for the skewed distribution
   rotmatrix = ellipse.rotate(assoc.UV,assoc.UW,assoc.VW)
   rotmatrix = rotmatrix.transpose()
   vec = np.dot(rotmatrix,np.asarray([a,b,c]))

   u = vec[0]+assoc.U
   v = vec[1]+assoc.V
   w = vec[2]+assoc.W

   # conversion between galactic and equatorial coordinates
   k = 4.74047     #Equivalent of 1 A.U/yr in km/s  
   a_g = np.array([[-0.0548755604, -0.8734370902,-0.4838350155],
                [+0.4941094279, -0.4448296300, 0.7469822445],
                [-0.8676661490, -0.1980763734, +0.4559837762]]) # rotation matrix for Galactic-->J2000
   radeg = 180/np.pi

   cosd = np.cos(dec_tmp/radeg)
   sind = np.sin(dec_tmp/radeg)
   cosa = np.cos(ra_tmp/radeg)
   sina = np.sin(ra_tmp/radeg)
    
   a_c = np.array([ [cosa*cosd,-sina,-cosa*sind],
                       [sina*cosd,cosa,-sina*sind],
                       [sind,0,cosd] ]) # rotation matrix for cartesian to spherical

   b = np.dot(a_g,a_c)

   vec1 = (b[0,0] * u + b[1,0] * v + b[2,0] * w)
   vec2 = (b[0,1] * u + b[1,1] * v + b[2,1] * w)
   vec3 = (b[0,2] * u + b[1,2] * v + b[2,2] * w)

   vrad = vec1
   pmra = vec2 / (k * 10.0)
   pmdec = vec3 / (k * 10.0)
   
   exp_rv = np.mean(vrad)
   exp_erv = np.std(vrad,ddof=1)
   exp_pmra = np.mean(pmra)
   exp_epmra = np.std(pmra,ddof=1)
   exp_pmdec = np.mean(pmdec)
   exp_epmdec = np.std(pmdec,ddof=1)

   return exp_pmra,exp_epmra,exp_pmdec,exp_epmdec,exp_rv,exp_erv


# given the ra, dec, and parallax (or expected distance), check the spatial agreement.
def spatial(mgp,ra,era,dec,edec,pi,epi,n_init):
   
   # first, get the XYZ position
   tra = ra + (np.random.randn(n_init)*era)*np.cos(dec*np.pi/180.)
   tdec = dec + np.random.randn(n_init)*edec
   tpi = pi + np.random.randn(n_init)*epi
   u,v,w,x,y,z = kinematics.gal_uvwxyz(ra=tra,dec=tdec,plx=tpi,pmra=np.zeros_like(tra),pmdec=np.zeros_like(tra),vrad=np.zeros_like(tra))

   # second: We now need to transform the XYZ position into the coordinate space of the randomly oriented ellipse.
   #  step 1: get the position of the target relative to the moving group
   star = [x-mgp.X,y-mgp.Y,z-mgp.Z]
   dist = np.sqrt(np.mean(star[0])**2+np.mean(star[1])**2+np.mean(star[2])**2)

   #  step 2: rotate that position into the moving group's coordinate system
   rotmatrix = ellipse.rotate(mgp.XY,mgp.XZ,mgp.YZ)
   rotmatrix = rotmatrix.transpose()

   pos = np.dot(rotmatrix,star)
   # pos now contains the position of the star relative to the center of the moving group, rotated to the frame of the moving group.
   
   # third (and finally): Get a sigma agreement.
   sigmas = np.sqrt((pos[0]/mgp.D)**2+(pos[1]/mgp.E)**2+(pos[2]/mgp.F)**2)
   sigma = np.mean(sigmas)
   esigma = np.std(sigmas,ddof=1)

   return sigma,esigma,dist

def weightadd(sigma,weight):
   sig = 0
   for g in range(len(sigma)):
      sig = sig + weight[g]**2*sigma[g]**2
   sig = np.sqrt(sig)/np.sum(weight)

   return sig

def gaus(x,a,sigma):
    return a*np.exp(-(x)**2/(2*sigma**2))

def gauscdf(x,m,sigma):
    return 0.5*(erfc((x-m)/(sigma*np.sqrt(2))))


#################################
## The Core LACEwING algorithm ##
#################################

def lacewing(moving_groups,young=None,iterate=None,ra=None,era=None,dec=None,edec=None,pmra=None,epmra=None,pmdec=None,epmdec=None,plx=None,eplx=None,rv=None,erv=None):
    if ra is None or dec is None:
        raise Exception('ERROR - You need at least an RA and a DEC for this to work')
    if pmra is None or pmdec is None:
        pmexists = 0
    else:
        pmexists = 1
    if plx is None:
        distexists = 0
    else:
        distexists = 1
    if rv is None:
        rvexists = 0
    else:
        rvexists = 1
    if iterate is None:
        iterate = 10000
    #print iterate

    output = []

    for i in range(len(moving_groups)):
        mgp = moving_groups[i]
        # this routine computes the expected pm (at 10 parsecs) and RV given an RA and DEC.
        exp_pmra,exp_epmra,exp_pmdec,exp_epmdec,exp_rv,exp_erv = converge(mgp,ra,era,dec,edec,iterate)
        # AR 2015.0526: This is now an astrometry routine only.
        exp_pm,exp_epm,exp_pa,exp_epa = astrometry.pmjoin(exp_pmra,exp_epmra,exp_pmdec,exp_epmdec)
        pm_new = None
        epm_new = None
        cosa = None
        sina = None
        exp_dist = 0
        exp_edist = 0
        pm_sig = None
        rv_sig = None
        dist_sig = None
        pos_sig = None
        pos_ksig = None
        pos_esig = None
        pos_eksig = None
        pos_sep = None
        pos_ksep = None
        sig = None

        sigma = []
        weight = []

        # step 1: compare proper motions. We're going to rotate the star's proper motion 
        #  vectors (pmRA, pmDEC) by the angle of the association's proper motion, to 
        #  transform the motion into pm= and pm_|_.
        if pmexists == 1:
            # 1a. Calculate the parallel and perpendicular components of the proper motion
            #   First: compute a rotation matrix to transform pmRA,pmDEC into pm= and pm_|_
            cosa = np.cos(exp_pa*np.pi/180.0)
            sina = np.sin(exp_pa*np.pi/180.0)
            rotmatrix = [ [ cosa, -sina ] , [sina, cosa] ]
          
            #   Second: Apply this rotation matrix to the measured and expected proper 
            #    motion vectors.
            pm_new = np.dot(rotmatrix,[pmra,pmdec])
            epm_new = np.dot(rotmatrix,[epmra,epmdec])
            exp_pm_new = np.dot(rotmatrix,[exp_pmra,exp_pmdec]) # should be [0,pm] if correctly done
            exp_epm_new = np.dot(rotmatrix,[exp_epmra,exp_epmdec])
            # The 0th element is the perpendicular component, the 1st is parallel
            
            # 1b. Compute a kinematic distance by comparing the magnitude of the 
            #   expected pm vector (which was calculated for a distance of 10 
            #   parsecs) to the measured parallel one
            exp_dist = exp_pm_new[1]/pm_new[1] * 10.0
            exp_edist = np.sqrt((exp_epm_new[1]/exp_pm_new[1])**2 + (epm_new[1]/pm_new[1])**2) * exp_dist

            # 1c. Now we use that kinematic distance to scale the error out. 
            #   Why? If the proper motion was in completely the wrong direction, 
            #   but tiny, it would come out as a match when it shouldn't. 
            #   This scales the error appropriately to stop that from happening.
            exp_epm_new = exp_epm_new * (10/exp_edist)

            # 1d. Compute the proper motion sigma based on the perpendicular 
            #   component (with scaled errors)
            pm_sig = np.sqrt(pm_new[0]**2)/(np.sqrt(epm_new[0]**2 + exp_epm_new[0]**2))

            # 1e. this should weed out anything where the pm vector was flipped (and
            #   distance was negative). I initially thought this should go under the
            #   "distance" section, but I can verify it entirely based on proper
            #   motions even when a distance isn't known, so it should go here.
            if ((abs(exp_epm_new[1]) < abs(exp_pm_new[1])) & (exp_dist < 0)):
                pm_sig = pm_sig + 1000.
                # I might also try this relative to a km/sec velocity dispersion, 
                #  pm_new[0]*exp_dist*4.74

            sigma.append(pm_sig)
            weight.append(1)

            pos_ksig,pos_eksig,pos_ksep = spatial(mgp,ra,era,dec,edec,1/exp_dist,exp_edist/(exp_dist**2),iterate)

        if (pmexists == 1) & (distexists == 0): # calculate the positional uncertainty for distance, but only use it if the distance does not exist. 
            sigma.append(pos_ksig)
            weight.append(1)

        # step 2: Distances exist. We can do the real test for spatial position.
        if (distexists == 1):
            pos_sig,pos_esig,pos_sep = spatial(mgp,ra,era,dec,edec,plx,eplx,iterate)
            sigma.append(pos_sig)
            weight.append(1)

        # step 3: Distances AND proper motions exist. Now we can compute a kinematic distance error.
        if (pmexists == 1) & (distexists == 1):
            # 3a. use the kinematic distance calculated in step 1c, which should have already run if pmexists==1.
            dist_sig = np.abs(1/plx - exp_dist) / np.sqrt((eplx/plx**2)**2 + exp_edist**2)
            sigma.append(dist_sig)
            weight.append(1)


        # step 4: RVs exist. This is straightforward: We have an estimated RV from the convergence.
        if (rvexists == 1):
            rv_sig = np.abs(rv - exp_rv) / np.sqrt(erv**2 + exp_erv**2)
            sigma.append(rv_sig)
            weight.append(1)

        sig = weightadd(sigma,weight)
                        
        #print sig
        if young == "young":
            coeff = mgp.coeff_young
        else:
            coeff = mgp.coeff_all

        #print coeff[0][0],astrometry.isnumber(coeff[0][0])
        if pmexists == 1:
            if distexists == 1:
                if rvexists == 1:
                    percent = gauscdf(sig,coeff[0][0],coeff[0][1])
                else:
                    percent = gauscdf(sig,coeff[4][0],coeff[4][1])
            else:
                if rvexists == 1:
                    percent = gauscdf(sig,coeff[5][0],coeff[5][1])
                else:
                    percent = gauscdf(sig,coeff[1][0],coeff[1][1])
        else:
            if distexists == 1:
                if rvexists == 1:
                    percent = gauscdf(sig,coeff[6][0],coeff[6][1])
                else:
                    percent = gauscdf(sig,coeff[2][0],coeff[2][1])
            else:
                if rvexists == 1:
                    percent = gauscdf(sig,coeff[3][0],coeff[3][1])
                else:
                    percent = 0
         
        percent = percent * 100

        output.append({'group':mgp.name, 'gof':sig, 'probability': percent, 'pmsig': pm_sig, 'kin_pmra':exp_pmra,'kin_epmra':exp_epmra,'kin_pmdec':exp_pmdec,'kin_epmdec':exp_epmdec,'distsig':dist_sig,'kin_dist':exp_dist, 'kin_edist':exp_edist, 'rvsig':rv_sig,'kin_rv':exp_rv,'kin_erv':exp_erv,'possig':pos_sig,'pos_esig':pos_esig,'pos_sep':pos_sep, 'posksig':pos_ksig,'pos_eksig':pos_eksig,'pos_ksep':pos_ksep})

    return output

#################################
## Default LACEwING csv loader ##
#################################

def csv_loader(infilename):

    file = open(infilename,'rb')
    readtable = ascii.get_reader(Reader=ascii.Basic)
    readtable.header.splitter.delimiter = ','
    readtable.data.splitter.delimiter = ','
    readtable.header.start_line = 0
    readtable.data.start_line = 1
    star = readtable.read(file)
    file.close()
    
    name = []
    coord = []
    era = []
    edec = []
    pmra = []
    epmra = []
    pmdec = []
    epmdec = []
    plx = []
    eplx = []
    rv = []
    erv = []
    note = []
    
    for i in np.arange(0,len(star)):
       
        name.append(star[i]['Name'])
        try:
            tcoord = SkyCoord(ra=star[i]['RA'],dec=star[i]['DEC'], unit=(u.degree,u.degree))
        except (ValueError,IndexError,KeyError):
            try:
                tcoord = SkyCoord(ra=star[i]['RAdeg'],dec=star[i]['DEdeg'], unit=(u.degree,u.degree))
            except (ValueError,IndexError,KeyError):
                tcoord = SkyCoord('{0:} {1:} {2:} {3:}{4:} {5:} {6:}'.format(star[i]['RAh'],star[i]['RAm'],star[i]['RAs'],star[i]['DE-'],star[i]['DEd'],star[i]['DEm'],star[i]['DEs']), unit=(u.hourangle, u.deg))

        coord.append(tcoord)

        try:
            tera = np.float(star[i]['eRA'])/1000./3600.
            tedec = np.float(star[i]['eDEC'])/1000./3600.
        except (ValueError,IndexError,KeyError):
            try:
                tera = np.float(star[i]['e_RAdeg'])/1000./3600.
                tedec = np.float(star[i]['e_DEdeg'])/1000./3600.
            except (ValueError,IndexError,KeyError):
                tera = 1.0/3600.
                tedec = 1.0/3600.
        era.append(tera)
        edec.append(tedec)
    
        try:
            tpmra = np.float(star[i]['pmRA'])/1000.
            tepmra = np.float(star[i]['epmRA'])/1000.
            tpmdec = np.float(star[i]['pmDEC'])/1000.
            tepmdec = np.float(star[i]['epmDEC'])/1000.
        except (ValueError,IndexError,KeyError):
            try:
                tpmra = np.float(star[i]['pmRA'])/1000.
                tepmra = np.float(star[i]['e_pmRA'])/1000.
                tpmdec = np.float(star[i]['pmDE'])/1000.
                tepmdec = np.float(star[i]['e_pmDE'])/1000.
            except (ValueError,IndexError,KeyError):
                try:
                    tpmra = np.float(star[i]['pmRA'])/1000.
                    tepmra = 0.01
                    tpmdec = np.float(star[i]['pmDEC'])/1000.
                    tepmdec = 0.01
                except (ValueError,IndexError,KeyError):
                    tpmra = None
                    tepmra = None
                    tpmdec = None
                    tepmdec = None
        pmra.append(tpmra)
        epmra.append(tepmra)
        pmdec.append(tpmdec)
        epmdec.append(tepmdec)
    
        try:
            tplx = np.float(star[i]['pi'])/1000.
            teplx = np.float(star[i]['epi'])/1000.
        except (ValueError, IndexError,KeyError):
            try:
                tplx = np.float(star[i]['plx'])/1000.
                teplx = np.float(star[i]['e_plx'])/1000.
            except (ValueError, IndexError,KeyError):
                tplx = None
                teplx = None
        plx.append(tplx)
        eplx.append(teplx)

        try:
            trv = np.float(star[i]['rv'])
            terv = np.float(star[i]['erv'])
        except (ValueError, IndexError,KeyError):
            try:
                trv = np.float(star[i]['HRV'])
                terv = np.float(star[i]['e_HRV'])
            except (ValueError, IndexError,KeyError):
                trv = None
                terv = None
        rv.append(trv)
        erv.append(terv)

        try:
            tnote = star[i]['Note']
        except(ValueError, IndexError,KeyError):
            tnote = None
        note.append(tnote)

    return name,coord,era,edec,pmra,epmra,pmdec,epmdec,rv,erv,plx,eplx,note

if __name__ == "__main__":
    try:
        infilename = argv[1]
    except IndexError: 
        print 'syntax: python lacewing.py inputfile [young] [outfilename] [verbose] [g.o.f.]'
    try:
        young = argv[2]
    except IndexError:
        young = 'field'
    try:
        outfilename = argv[3]
    except IndexError:
        outfilename = 'lacewing_output.csv'
    try:
        verbose = argv[4]
    except IndexError:
        verbose = 'quiet'
    try:
        percentages = argv[5]
    except IndexError:
        percentages = 'percentage'

    lineno = 1

    name,coord,era,edec,pmra,epmra,pmdec,epmdec,rv,erv,plx,eplx,note = csv_loader(infilename)

    moving_groups = moving_group_loader()

    moving_groups = [moving_groups[x] for x in range(len(moving_groups)) if moving_groups[x].name != "Field"]

    outfile = open(outfilename,'wb')
    if verbose == "verbose":
        # verbose output is one line per moving group per entry (1 line per group per star) in CSV format
        outfile.write('lineno,Name,RA,DEC,Group,d1,Probability,d2,sig_pm,kin_pmra,kin_epmra,kin_pmdec,kin_epmdec,pmra,epmra,pmdec,epmdec,d3,sig_dist,kin_dist,kin_edist,dist,edist,d4,sig_rv,kin_rv,kin_erv,rv,erv,d5,sig_pos,sep,Note,.\n')
    else:
        # regular output is a one line per entry summary, also in .csv form
        outfile.write('Name,Note,Group,Probability,Predicted Dist,Predicted Dist uncertainty,Predicted RV,Predicted RV uncertainty,')
        for i in range(len(moving_groups)):
            outfile.write('{0:},'.format(moving_groups[i].name))
        outfile.write('\n')
    


    for i in range(len(coord)):
        out = lacewing(moving_groups, young=young, ra=coord[i].ra.degree,era=era[i],dec=coord[i].dec.degree,edec=edec[i],pmra=pmra[i],epmra=epmra[i],pmdec=pmdec[i],epmdec=epmdec[i],rv=rv[i],erv=erv[i],plx=plx[i],eplx=eplx[i])

        # if verbose output has been selected, print our strings to a file here
        if verbose == "verbose":
            for j in range(len(out)):
                outfile.write('{0:},{1:},{2:09.5f},{3:+08.5f},{4:12},'.format(lineno,name[i],coord[i].ra.degree,coord[i].dec.degree,out[j]['group']))
                if percentages == "percentage":
                    outfile.write('PROB=,{0: 8.0f},'.format(out[j]['probability']))
                else:
                    outfile.write('GOF=,{0: 6.2f},'.format(out[j]['gof']))
                if pmra[i] != None:
                    outfile.write('PM=,{0: 6.2f}, {1:+6.2f},{2:6.2f}, {3:+6.2f},{4:6.2f}, {5:+6.2f},{6:6.2f}, {7:+6.2f},{8:6.2f},'.format(out[j]['pmsig'],out[j]['kin_pmra']*1000.,out[j]['kin_epmra']*1000.,out[j]['kin_pmdec']*1000.,out[j]['kin_epmdec']*1000.,pmra[i]*1000.,epmra[i]*1000.,pmdec[i]*1000.,epmdec[i]*1000.))
                else:
                    outfile.write('PM=,, {0:+6.2f},{1:6.2f},{2:+6.2f},{3:6.2f},,,,,'.format(out[j]['kin_pmra']*1000.,out[j]['kin_epmra']*1000.,out[j]['kin_pmdec']*1000.,out[j]['kin_epmdec']*1000.))                    
                if (plx[i] != None) & (pmra[i] != None):
                    outfile.write('DIST=,{0: 6.2f}, {1: 8.2f},{2:7.2f},{3:+8.2f},{4:8.2f},'.format(out[j]['distsig'],out[j]['kin_dist'],out[j]['kin_edist'],1/plx[i],eplx[i]/(plx[i]**2)))
                elif (plx[i] == None) & (pmra[i] != None):
                    outfile.write('DIST=,, {0: 8.2f},{1:7.2f},,,'.format(out[j]['kin_dist'],out[j]['kin_edist']))
                elif (plx[i] != None) & (pmra[i] == None):
                    outfile.write('DIST=,,,, {0: 8.2f},{1:7.2f},'.format(1/plx[i],eplx[i]/(plx[i]**2)))
                elif (plx[i] == None) & (pmra[i] == None):
                    outfile.write('DIST=,,,,,,')
                if rv[i] != None:
                    outfile.write('RV=,{0: 6.2f}, {1:+8.2f},{2:8.2f},{3:+8.2f},{4:8.2f},'.format(out[j]['rvsig'],out[j]['kin_rv'],out[j]['kin_erv'],rv[i],erv[i]))
                else:
                    outfile.write('RV=,, {0:+8.2f},{1:8.2f},,,'.format(out[j]['kin_rv'],out[j]['kin_erv']))
                if (plx[i] != None):
                    outfile.write('POS=,{0: 6.2f},{1: 8.2f},'.format(out[j]['possig'],out[j]['pos_sep']))
                else:
                    if (plx[i] == None) & (pmra[i] != None):
                        outfile.write('KPOS=,{0: 6.2f},{1: 8.2f},'.format(out[j]['posksig'],out[j]['pos_ksep']))
                    else:
                        outfile.write('POS=,,,')
                outfile.write('{0:},\n'.format(note[i]))
                lineno = lineno+1

        else:
            # We've finished processing one star. If regular output was selected, create and print the output string.
            probs = [out[j]['probability'] for j in range(len(out))]
            order = np.argsort(probs)[::-1]
            # if even the best match isn't above the threshold of consideration:
            if out[order[0]]['probability'] < 20:
                outfile.write('{0:},{1:},(None),,,,,,'.format(name[i],note[i]))
            else:
                outfile.write('{0:},{1:},{2:},{3: 5.0f},{4: 8.2f},{5:7.2f},{6:+6.2f},{7: 5.2f},'.format(name[i],note[i],out[order[0]]['group'],out[order[0]]['probability'],out[order[0]]['kin_dist'],out[order[0]]['kin_edist'],out[order[0]]['kin_rv'],out[order[0]]['kin_erv']))
            for k in range(len(out)):
                outfile.write('{0: 5.0f},'.format(out[k]['probability']))
            outfile.write('\n')

    outfile.close()

