import numpy
import scipy
from astropy.io import ascii
from sys import argv
import kinematics
import converge
import astrometry

###############################################
###############################################
### MAIN ROUTINE
###############################################
###############################################

###############################################
### LACEwING 1.2
### ARR 2015-05-26
### 1.0: Everything
### 1.1: Removed dependence on GROUP columns in input data table
###      Fixed dependency on astrometry module
###      Default output should be a summary mode now
### 1.2: Fixed input format
###      All astrometry routines now point to astrometry module
###      Number of groups loaded is now controlled by input .csv file
###############################################

try:
    infilename = argv[1]
except IndexError: 
    print 'syntax: python lacewing.py inputfile young=young outname=outfile verbose=verbose rawgofs=percentage'
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

        self.coeff_all = [ [mgp['field_all_M'],mgp['field_all_S']],[mgp['field_pm_M'],mgp['field_pm_S']],[mgp['field_dist_M'],mgp['field_dist_S']],[mgp['field_rv_M'],mgp['field_rv_S']],[mgp['field_pmdist_M'],mgp['field_pmdist_S']],[mgp['field_pmrv_M'],mgp['field_pmrv_S']],[mgp['field_distrv_M'],mgp['field_distrv_S']]]
        self.coeff_young = [ [mgp['young_all_M'],mgp['young_all_S']],[mgp['young_pm_M'],mgp['young_pm_S']],[mgp['young_dist_M'],mgp['young_dist_S']],[mgp['young_rv_M'],mgp['young_rv_S']],[mgp['young_pmdist_M'],mgp['young_pmdist_S']],[mgp['young_pmrv_M'],mgp['young_pmrv_S']],[mgp['young_distrv_M'],mgp['young_distrv_S']]]
        for x in range(len(self.coeff_all)):
            if astrometry.isnumber(self.coeff_all[x][0]):
                self.coeff_all[x][0] = numpy.float(self.coeff_all[x][0])
            if astrometry.isnumber(self.coeff_all[x][1]):
                self.coeff_all[x][1] = numpy.float(self.coeff_all[x][1])
        for x in range(len(self.coeff_young)):
            if astrometry.isnumber(self.coeff_young[x][0]):
                self.coeff_young[x][0] = numpy.float(self.coeff_young[x][0])
            if astrometry.isnumber(self.coeff_young[x][1]):
                self.coeff_young[x][1] = numpy.float(self.coeff_young[x][1])

# read in ellipse parameters for associations

filename = argv[1]
file = open('Moving_Group_all.csv','rb')
readtable = ascii.get_reader(Reader=ascii.Basic)
readtable.header.splitter.delimiter = ','
readtable.data.splitter.delimiter = ','
readtable.header.start_line = 0
readtable.data.start_line = 1
groups = readtable.read(file)
file.close()

moving_groups = []

# AR 2015.0326: Number of moving groups is now controlled by the input file. (Last one should be "field" and not included)
for i in range(len(groups)-1):
    moving_groups.append(Mgp(groups[i]))

# AR 2015.0326: Input file format is now a .csv file with a one line header.
file = open(infilename,'rb')
readtable = ascii.get_reader(Reader=ascii.Basic)
readtable.header.splitter.delimiter = ','
readtable.data.splitter.delimiter = ','
readtable.header.start_line = 0
readtable.data.start_line = 1
star = readtable.read(file)
file.close()

outfile = open(outfilename,'wb')
if verbose == "verbose":
    # verbose output is one line per moving group per entry (14 lines per star) in CSV format
    outfile.write('lineno,Name,RA,DEC,pmRA,pmDEC,MGP,d1,sig,d2,sig_pm,pm_|_,exp_pm_|_,d3,sig_dist,dist,exp_dist,exp_edist,d4,sig_rv,rv,exp_rv,exp_erv,d5,sig_pos,mgp1,mgp2,mgp3,.\n')
else:
    # regular output is a one line per entry summary, also in .csv form
    outfile.write('Name,Group,Percent,Predicted Dist,Predicted Dist uncertainty,Predicted RV,Predicted RV uncertainty,eps_Cha,eta_Cha,TW_Hya,beta_Pic,Octans,Tuc-Hor,Columba,Argus,AB_Dor,Pleiades,Her-Lyr,Coma_Ber,Ursa_Major,Hyades\n')

lineno = 1

for i in numpy.arange(0,len(star)):
   
    name = star[i]['Name']
    try:
        ra = float(star[i]['RA'])
        dec = float(star[i]['DEC'])
    except (ValueError,IndexError,KeyError):
        try:
            ra = astrometry.ten((float(star[i]['Rah']),float(star[i]['Ram']),float(star[i]['Ras'])))*15.
            dec = astrometry.ten((numpy.abs(float(star[i]['DECd'])),float(star[i]['DECm']),float(star[i]['DECs'])))
            if star[i]['DECf'].strip() == '-':
                dec = dec * -1.0
        except ValueError:
            ra = astrometry.ten((float(star[i]['Rah']),float(star[i]['Ram']),float(star[i]['Ras'])))*15.
            dec = astrometry.ten((numpy.abs(float(star[i]['DECd'])),float(star[i]['DECm'])))
            if star[i]['DECf'] == '-':
                dec = dec * -1.0
    try:
        era = numpy.float(star[i]['eRA'])/1000./3600.
        edec = numpy.float(star[i]['eDEC'])/1000./3600.
    except (ValueError,IndexError,KeyError):
        era = 1.0/3600.
        edec = 1.0/3600.

    pmexists = 0
    try:
        pmra = numpy.float(star[i]['pmRA'])/1000.
        epmra = numpy.float(star[i]['epmRA'])/1000.
        pmdec = numpy.float(star[i]['pmDEC'])/1000.
        epmdec = numpy.float(star[i]['epmDEC'])/1000.
        pmexists = 1
    except ValueError:
      try:
          pmra = numpy.float(star[i]['pmRA'])/1000.
          epmra = 0.01
          pmdec = numpy.float(star[i]['pmDEC'])/1000.
          epmdec = 0.01
          pmexists = 1
      except ValueError:
          pmra = 999.99
          epmra = 999.99
          pmdec = -999.99
          epmdec = 999.99
          pa = 999.99
          pm = 999.99

    distexists=0
    rvexists=0
    try:
        plx = numpy.float(star[i]['pi'])/1000.
        eplx = numpy.float(star[i]['epi'])/1000.
        distexists = 1
    except (ValueError, IndexError):
        distexists=0
    try:
        rv = numpy.float(star[i]['rv'])
        erv = numpy.float(star[i]['erv'])
        rvexists = 1
    except (ValueError, IndexError):
        rvexists=0

#    print "{8:}  {0:09.5f} {1:7.5f}  {2:+09.5f} {3:7.5f}  {4:+9.5f} {5:8.5f}  {6:+9.5f} {7:8.5f}".format(ra,era,dec,edec,pmra,epmra,pmdec,epmdec,name)

    #print "Running Convergence"

    # these lists will only be needed to save output for the regular summary.
    if verbose != "verbose":
        matchgroup = []
        matchsig = []
        matchdist = []
        matchedist = []
        matchrv = []
        matcherv = []

    for i in range(len(moving_groups)):
        mgp = moving_groups[i]
        # this routine computes the expected pm (at 10 parsecs) and RV given an RA and DEC.
        exp_pmra,exp_epmra,exp_pmdec,exp_epmdec,exp_rv,exp_erv = converge.converge(mgp,ra,era,dec,edec,100000)
        # AR 2015.0526: This is now an astrometry routine only.
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

        sigma = []

        pm_string = 'PM: s=,"",   "",    {0:+6.2f}, '.format(0) # ok, cheating, but the expected pm_|_ is 0 by definition.
        dist_string = 'Dist: s=, "",    "",    "",    "", '
        rv_string = 'RV: s=,"",    "",    {0:+8.2f}, {1: 7.2f},'.format(exp_rv,exp_erv)
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
          
            #   Second: Apply this rotation matrix to the measured and expected proper 
            #    motion vectors.
            pm_new = numpy.dot(rotmatrix,[pmra,pmdec])
            epm_new = numpy.dot(rotmatrix,[epmra,epmdec])
            exp_pm_new = numpy.dot(rotmatrix,[exp_pmra,exp_pmdec]) # should be [0,pm] if correctly done
            exp_epm_new = numpy.dot(rotmatrix,[exp_epmra,exp_epmdec])
            # The 0th element is the perpendicular component, the 1st is parallel
            
            # 1b. Compute a kinematic distance by comparing the magnitude of the 
            #   expected pm vector (which was calculated for a distance of 10 
            #   parsecs) to the measured parallel one
            exp_dist = exp_pm_new[1]/pm_new[1] * 10.0
            exp_edist = numpy.sqrt((exp_epm_new[1]/exp_pm_new[1])**2 + (epm_new[1]/pm_new[1])**2) * exp_dist

            # 1c. Now we use that kinematic distance to scale the error out. 
            #   Why? If the proper motion was in completely the wrong direction, 
            #   but tiny, it would come out as a match when it shouldn't. 
            #   This scales the error appropriately to stop that from happening.
            exp_epm_new = exp_epm_new * (10/exp_edist)

            # 1d. Compute the proper motion sigma based on the perpendicular 
            #   component (with scaled errors)
            pm_sig = numpy.sqrt(pm_new[0]**2)/(numpy.sqrt(epm_new[0]**2 + exp_epm_new[0]**2))

            # 1e. this should weed out anything where the pm vector was flipped (and
            #   distance was negative). I initially thought this should go under the
            #   "distance" section, but I can verify it entirely based on proper
            #   motions even when a distance isn't known, so it should go here.
            if ((abs(exp_epm_new[1]) < abs(exp_pm_new[1])) & (exp_dist < 0)):
                pm_sig = pm_sig + 1000.
                # I might also try this relative to a km/sec velocity dispersion, 
                #  pm_new[0]*exp_dist*4.74

            sigma.append(pm_sig)
            pm_string = 'PM: s=,{0: 8.2f}, {1:+6.2f},{2:+6.2f}, '.format(pm_sig,pm_new[0],exp_pm_new[0])

            dist_string = 'Dist: s=,"",    "",    {0:+8.2f}, {1: 7.2f},'.format(exp_dist,exp_edist)
            converge_plx = 1/exp_dist
            converge_eplx = exp_edist/(exp_dist**2)

        # step 2: Distances exist. Without a proper motion, we can't do anything.
        if distexists == 1:
            dist_string = 'Dist: s=,"",   {0:+8.2f},    "",    "",'.format(1/plx)
            converge_plx = plx
            converge_eplx = eplx

        # step 3: Distances AND proper motions exist. Now we can compute a kinematic distance error.
        if (pmexists == 1) & (distexists == 1):

            # 3a. use the kinematic distance calculated in step 1c, which should have already run if pmexists==1.
            dist_sig = numpy.abs(1/plx - exp_dist) / numpy.sqrt((eplx/plx**2)**2 + exp_edist**2)

            sigma.append(dist_sig)
            dist_string = 'Dist: s=,{0: 8.2f}, {1:+8.2f},{2:+8.2f},{3: 7.2f},'.format(dist_sig,1/plx,exp_dist,exp_edist)


        # step 4: RVs exist. This is straightforward: We have an estimated RV from the convergence.
        if (rvexists == 1):
            rv_sig = numpy.abs(rv - exp_rv) / numpy.sqrt(erv**2 + exp_erv**2)

            sigma.append(rv_sig)
            rv_string = 'RV: s=,{0: 8.2f},   {1:+8.2f},{2:+8.2f},{3: 7.2f},'.format(rv_sig,rv,exp_rv,exp_erv)

        # step 5: A test for spatial position. I feed in RA, DEC, and the distance (either trig
        #   parallax or kinematic distance) and converge.spatial() tells me how close it is to
        #   the ellipse describing the moving group.
        if (pmexists ==1) | (distexists == 1):
            pos_sig,pos_esig = converge.spatial(mgp,ra,era,dec,edec,converge_plx,converge_eplx,10000)

            sigma.append(pos_sig)
            pos_string = 'XYZ: s=,{0: 8.2f},'.format(pos_sig)


        sig = converge.weightadd(sigma)
                        
        #print sig
        if young == "young":
            coeff = mgp.coeff_young
        else:
            coeff = mgp.coeff_all

        #print coeff[0][0],astrometry.isnumber(coeff[0][0])
        percent = 0
        if percentages == 'percentage':
            if pmexists == 1:
                if distexists == 1:
                    if rvexists == 1:
                        if numpy.logical_not(astrometry.isnumber(coeff[0][0])):
                            if sig < coeff[0][1]:
                                percent = 1
                        else:
                            percent = converge.gauscdf(sig,coeff[0][0],coeff[0][1])
                    else:
                        if numpy.logical_not(astrometry.isnumber(coeff[4][0])):
                            if sig < coeff[4][1]:
                                percent = 1
                        else:
                            percent = converge.gauscdf(sig,coeff[4][0],coeff[4][1])
                else:
                    if rvexists == 1:
                        if numpy.logical_not(astrometry.isnumber(coeff[5][0])):
                            if sig < coeff[5][1]:
                                percent = 1
                        else:
                            percent = converge.gauscdf(sig,coeff[5][0],coeff[5][1])
                    else:
                        if numpy.logical_not(astrometry.isnumber(coeff[1][0])):
                            if sig < coeff[1][1]:
                                percent = 1
                        else:
                            percent = converge.gauscdf(sig,coeff[1][0],coeff[1][1])
            else:
                if distexists == 1:
                    if rvexists == 1:
                        if numpy.logical_not(astrometry.isnumber(coeff[6][0])):
                            if sig < coeff[6][1]:
                                percent = 1
                        else:
                            percent = converge.gauscdf(sig,coeff[6][0],coeff[6][1])
                    else:
                        if numpy.logical_not(astrometry.isnumber(coeff[2][0])):
                            if sig < coeff[2][1]:
                                percent = 1
                        else:
                            percent = converge.gauscdf(sig,coeff[2][0],coeff[2][1])
                else:
                    if rvexists == 1:
                        if numpy.logical_not(astrometry.isnumber(coeff[3][0])):
                            if sig < coeff[3][1]:
                                percent = 1
                        else:
                            percent = converge.gauscdf(sig,coeff[3][0],coeff[3][1])
         
            percent = percent * 100
        else:
            # if percentage is set to anything other than "percentage", output the raw goodness-of-fit.
            percent = sig
      

        # remove spaces from name - helps with programming later on
        name = name.replace(' ', '_')
        mgpname = mgp.name.replace(' ', '_')

        #print '{11:},"{0:16}",{1:09.5f},{2:+08.5f},{3:+6.4f},{4:+6.4f},{5:12},SIG=,{6: 8.2f},{7:} {8:} {9:} {10:}'.format(name,ra,dec,pmra,pmdec,mgpname,percent,pm_string,dist_string,rv_string,pos_string,lineno)
        # if verbose output has been selected, print our strings to a file here
        if verbose == "verbose":
            outfile.write('{11:},"{0:16}",{1:09.5f},{2:+08.5f},{3:+6.4f},{4:+6.4f},{5:12},SIG=,{6: 8.2f},{7:} {8:} {9:} {10:}, .\n'.format(name,ra,dec,pmra,pmdec,mgpname,percent,pm_string,dist_string,rv_string,pos_string,lineno))
        else:
            # if regular output was selected, save the info for each moving group for consideration once they're all done.
            matchgroup.append(mgpname)
            matchsig.append(percent)
            matchdist.append(exp_dist)
            matchedist.append(exp_edist)
            matchrv.append(exp_rv)
            matcherv.append(exp_erv)

        lineno = lineno+1
    # We've finished processing one star. If regular output was selected, create and print the output string.
    if verbose != "verbose":
        order = numpy.argsort(matchsig)[::-1]
        # if even the best match isn't above the threshold of consideration:
        if matchsig[order[0]] < 20:
            outfile.write('{0:},(None),    ,    ,    ,    ,    ,'.format(name))
        else:
            outfile.write('{0:},{1:},{2: 5.2f},{3: 6.2f},{4: 5.2f},{5:+6.2f},{6: 5.2f},'.format(name,matchgroup[order[0]],matchsig[order[0]],matchdist[order[0]],matchedist[order[0]],matchrv[order[0]],matcherv[order[0]]))
        for j in range(len(moving_groups)):
            outfile.write('{0: 5.2f},'.format(matchsig[j]))
        outfile.write('\n')

outfile.close()
