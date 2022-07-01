import os
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units
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
### LACEwING 2.0
### ARR 2019-12-23
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
### 2.0: Complete refactor as installable Python module
###############################################

# conversion between galactic and equatorial coordinates
K = 4.74047     #Equivalent of 1 A.U/yr in km/s  
R_GAL = np.array([[-0.0548755604, -0.8734370902,-0.4838350155],
            [+0.4941094279, -0.4448296300, 0.7469822445],
            [-0.8676661490, -0.1980763734, +0.4559837762]]) # rotation matrix for Galactic-->J2000
dataformat = [
    {'name': 'name', 'quant': 1, 'default': 'name', 'forcalc': None, 'keys': ['name', 'source', 'id', 'target']},
    {'name': 'ra', 'quant': units.deg, 'default': "has_coord", 'forcalc': units.deg, 'keys': ['raj2000', 'radeg', 'ra', 'ra_icrs', 'alpha', 'r.a.']},
    {'name': 'dec', 'quant': units.deg, 'default': "has_coord", 'forcalc': units.deg, 'keys': ['dej2000', 'dec', 'dedeg', 'decdeg', 'decj2000', 'de_icrs', 'delta', 'decl.']},
    {'name': 'era', 'quant': units.mas, 'default': 1/3600*units.deg, 'forcalc': units.deg, 'keys': ['e_raj2000', 'era', "e_ra", "e_radeg", "e_ra_icrs", "e_alpha", "ura", "u_ra", "u_radeg", "u_raj2000", "u_alpha"]},
    {'name': 'edec', 'quant': units.mas, 'default': 1/3600*units.deg, 'forcalc': units.deg, 'keys': ['e_dej2000', 'edec', "e_dec", "e_dedeg", "e_decdeg", "e_dec_icrs", "e_decj2000", "e_de_icrs", "e_delta", "udec", "u_dec", "u_dedeg", "u_dej2000", "u_de_icrs", "u_decdeg", "u_decj2000", "u_dec_icrs", "u_delta"]},
    {'name': 'pmra', 'quant': units.mas/units.yr, 'default': "has_pm", 'forcalc': units.arcsec/units.yr, 'keys': ['pmra', "pm_ra", "pm_ra_cos_dec", "pm_alpha", "mu_a", "mu_alpha", "mu_alpha_cos_delta", "mu_ra"]},
    {'name': 'pmdec', 'quant': units.mas/units.yr, 'default': "has_pm", 'forcalc': units.arcsec/units.yr, 'keys': ['pmdec', "pm_dec", "pmde", "pm_de", "pm_delta", "mu_d", "mu_delta", "mu_dec", "mu_de"]},
    {'name': 'epmra', 'quant': units.mas/units.yr, 'default': 0.001*units.arcsec/units.yr, 'forcalc': units.arcsec/units.yr, 'keys': ['epmra', "e_pmra", "e_pm_ra", "e_pm_ra_cos_dec", "e_pm_alpha", "e_mu_a", "e_mu_a_cos_d", "e_mu_alpha", "e_mu_alpha_cos_delta", "e_mu_ra"]},
    {'name': 'epmdec', 'quant': units.mas/units.yr, 'default': 0.001*units.arcsec/units.yr, 'forcalc': units.arcsec/units.yr, 'keys': ['epmdec', "e_pmdec", "e_pm_dec", "e_pmde", "e_pm_de", "e_pm_delta", "e_mu_d", "e_mu_delta", "e_mu_dec", "e_mu_de", "u_pmdec", "u_pm_dec", "u_pm_de", "u_pm_delta", "u_mu_d", "u_mu_delta", "u_mu_dec", "u_mu_de"]},
    {'name': 'plx', 'quant': units.mas, 'default': "has_plx", 'forcalc': units.arcsec, 'keys': ['pi', 'plx']},
    {'name': 'eplx', 'quant': units.mas, 'default': 0.001*units.arcsec, 'forcalc': units.arcsec, 'keys': ['epi', 'eplx', 'e_pi', 'e_plx', 'u_plx', 'u_plx']},
    {'name': 'rv', 'quant': units.km/units.s, 'default': "has_rv", 'forcalc': units.km/units.s, 'keys': ['hrv', 'rv', 'r.v.', 'vrad']},
    {'name': 'erv', 'quant': units.km/units.s, 'default': 1*units.km/units.s, 'forcalc': units.km/units.s, 'keys': ['e_hrv', 'erv', 'e_rv', 'e_r.v.', 'e_vrad', 'u_hrv', 'urv', 'u_rv', 'u_r.v.', 'u_vrad']},
    {'name': 'note', 'quant': 1, 'default': "note", 'forcalc': None, 'keys': ['note', "comment"]}
]



class MovingGroup:
    def __init__(self, definition, n_iter):
        self.n_iter = n_iter

        self.definition = definition
        # everything in the definition dictionary should become a class attribute
        self.name = self._load(["Name"])
        self.color = tuple(self._load(["Red", "Green", "Blue"]))
        self.param3d = self._load(["U", "A", "V", "B", "W", "C", "UV", "UW", "VW", "X", "D", "Y", "E", "Z", "F", "XY", "XZ", "YZ"])
        self.param2d = self._load(["U", "A2", "V", "B2", "W", "C2", "UV2", "UW2", "VW2", "X", "D2", "Y", "E2", "Z", "F2", "XY2", "XZ2", "YZ2"])
        self.weightednumber = self._load(["Weightednumber"])
        self.uniform = self._load(["uniform"])
        self.coeff_field = self._load_coeff('field')
        self.coeff_young = self._load_coeff('young')

    def _load(self,values):
        if len(values) > 1:
            output = {}
            for val in values:
                output[val] = self.definition[val]
            #output = tuple(output)
        else:
            output = self.definition[values[0]]
        return output

    def _load_coeff(self,type):
        output = {}
        for datatype in ['all', 'pm', 'dist', 'rv', 'pmdist', 'pmrv', 'distrv']:
            output[datatype] = [self.definition["{}_{}_M".format(type,datatype)],self.definition["{}_{}_S".format(type,datatype)]]

        return output

    def simulate(self):
        # generate random points within the moving group distribution
        a = self.param3d['A'] * np.random.randn(self.n_iter)
        b = self.param3d['B'] * np.random.randn(self.n_iter)
        c = self.param3d['C'] * np.random.randn(self.n_iter)
        # rotate them to account for the skewed distribution
        rotmatrix = ellipse.rotate(self.param3d['UV'],self.param3d['UW'],self.param3d['VW'])
        rotmatrix = rotmatrix.transpose()
        vec = np.dot(rotmatrix,np.asarray([a,b,c]))

        self.sim = [vec[0]+self.param3d['U'], vec[1]+self.param3d['V'], vec[2]+self.param3d['W']]

class StarClass:
    def __init__(self, stardata, n_iter):
        self.has_coord = True
        self.has_pm = True
        self.has_plx = True
        self.has_rv = True

        self.recognizer(stardata)

        if not self.has_coord:
            raise ValueError("LACEwING needs an RA and DEC to run.")

        self.n_iter = n_iter


    def recognizer(self, stardata):
        for item in dataformat:
            if item['name'] in stardata:
                print(type(stardata[item['name']]), stardata[item['name']])
                if type(stardata[item['name']]) == units.Quantity:
                    if stardata[item['name']].value != np.nan:
                        setattr(self, item['name'], stardata[item['name']])
                    else:
                        self._fill_default(stardata, item)
                else:
                    if stardata[item['name']] != np.nan:
                        # if the item is listed, give it the attribute converted to the correct quantity
                        if item['forcalc'] == None: # The "None"s are strings
                            setattr(self, item['name'], stardata[item['name']])
                        else:
                            setattr(self, item['name'], stardata[item['name']] * item['default'])
                    else:
                        self._fill_default(stardata, item)
            else:
                print(item['name'])
                self._fill_default(stardata, item)

    def _fill_default(self, stardata, item):
        # if the item is unlisted and default is a quantity, set the attribute to the default
        if type(item['default']) == units.Quantity:
            setattr(self, item['name'], item['default'].value)
        elif item['forcalc'] == None:
            setattr(self, item['name'], "")
        else:
            # else, set the attribute listed to False
            setattr(self, item['default'], False)

    def simulate(self):
        dec_tmp = self.dec + np.random.randn(self.n_iter) * self.edec
        ra_tmp = self.ra + np.random.randn(self.n_iter) * self.era / np.cos(np.deg2rad(dec_tmp))

        cosd = np.cos(np.deg2rad(dec_tmp))
        sind = np.sin(np.deg2rad(dec_tmp))
        cosa = np.cos(np.deg2rad(ra_tmp))
        sina = np.sin(np.deg2rad(ra_tmp))
            
        a_c = np.array([ [cosa*cosd,-sina,-cosa*sind],
                            [sina*cosd,cosa,-sina*sind],
                            [sind,0,cosd] ]) # rotation matrix for cartesian to spherical

        self.b = np.dot(R_GAL,a_c)

        
    def xyz(self, parallax, eparallax=0.001*units.arcsec):
        # first, get the XYZ position
        tdec = self.dec + np.random.randn(self.n_iter) * self.edec
        tra = self.ra + (np.random.randn(self.n_iter) * self.era) * np.cos(np.deg2rad(tdec))
        tplx = parallax + np.random.randn(self.n_iter) * eparallax

        print(tra.unit,tdec.unit,tplx.unit)
        u,v,w,self.x,self.y,self.z = kinematics.gal_uvwxyz(ra=tra, dec=tdec, plx=tplx, 
                                                           pmra=np.zeros(len(tra))*units.arcsecond/units.year,
                                                           pmdec=np.zeros(len(tra))*units.arcsecond/units.year,
                                                           vrad=np.zeros(len(tra))*units.km/units.second)

class LACEwING:
    def __init__(self, config, iterate, young):
        # load all the moving group identifications with pandas
        grouptable = pd.read_csv(config)
        # convert to records
        groups = grouptable.to_dict(orient='records')
        self.n_iter = iterate
        self.young = young
        self.groupdata = []
        self.mgplist = []
        for group in groups:
            if group['Name'] != "Field":
                mgp = MovingGroup(group, iterate)
                mgp.simulate()
                self.groupdata.append(mgp)
                self.mgplist.append(mgp.name)
        self.columns=["name","ra","dec","probability","gof","sig_pm","p_pmra",
                      "p_epmra","p_pmdec","p_epmdec","pmra","epmra","pmdec","epmdec",
                      "sig_dist","k_dist","k_edist","dist","edist","sig_rv","p_rv",
                      "p_erv","rv","erv","sig_pos","sig_epos","sep","note"]

    def calculate(self,stardata):
        """
        This is the main implementation of the LACEwING algorithm.
        """

        star = StarClass(stardata, self.n_iter)
        star.simulate()
        if star.has_plx:
            star.xyz(star.plx, star.eplx)

        self.result_frame = pd.DataFrame(columns=self.columns, index=self.mgplist)
        for column in self.columns:
            if column in star.__dict__:
                self.result_frame[column] = star.__dict__[column]

        
        for mgp in self.groupdata:
            sigma = []
            weights = []

            # new model: result_frame is a pandas DataFrame of the required size with empty cells. 
            # Each part of the computation will fill and read the table.

            # xyz cannot be a required feature of a star, because it might not have a parallax.

            # spatial must be able to do an xyz computation using the kinematic distance.

            # carry units through the entire computation


            self.converge(mgp,star)
            if star.has_pm:
                sig_pm, weight = self.pmcompare(mgp, star)
                sigma.append(sig_pm)
                weights.append(weight)
            if not star.has_plx and star.has_pm:
                star.xyz(self.resultframe.loc[mgp.name, "k_dist"], self.resultframe.loc[mgp.name, "k_edist"])
                sig_pos, weight = self.spatial(mgp, star)
                sigma.append(sig_pos)
                weights.append(weight)
            if star.has_plx:
                sig_pos, weight = self.spatial(mgp, star)
                sigma.append(sig_pos)
                weights.append(weight)
            if star.has_plx and star.has_pm:
                # step 3: Distances AND proper motions exist. Now we can compute a kinematic distance error.
                # 3a. use the kinematic distance calculated in step 1c, which should have already run if pmexists==1.
                sig_dist, weight = self.distcompare(mgp, star)
                sigma.append(sig_dist)
                weights.append(1)
            # step 4: RVs exist. This is straightforward: We have an estimated RV from the convergence.
            if star.has_rv:
                sig_rv, weight = self.rvcompare(mgp, star)
                sigma.append(sig_rv)
                weights.append(weight)

            datatype = ""
            # create the coefficient name to query
            if star.has_pm:
                datatype += "pm"
            if star.has_plx:
                datatype += "dist"
            if star.has_rv:
                datatype += "rv"
            if datatype == "pmdistrv":
                datatype = "all"

            if datatype != "":
                if self.young:
                    coeff = mgp.coeff_young[datatype]
                else:
                    coeff = mgp.coeff_field[datatype]
                self.result_frame.loc[mgp.name, "gof"] = self._weightadd(sigma,weights)
                self.result_frame.loc[mgp.name, "probability"] = self._gauscdf(self.result_frame.loc[mgp.name, "gof"],coeff[0],coeff[1])*100
            else:
                self.result_frame.loc[mgp.name, "gof"] = np.nan
                self.result_frame.loc[mgp.name, "probability"] = np.nan

        print(self.result_frame.loc[mgp.name, "probability"], star.has_pm, star.has_plx, star.has_rv, sigma)

        return self.result_frame
            

    # Calculate the proper motion vector at the star's RA and DEC and a distance of 10 pc
    def converge(self,mgp,star):

        vec1 = (star.b[0,0] * mgp.sim[0] + star.b[1,0] * mgp.sim[1] + star.b[2,0] * mgp.sim[2])
        vec2 = (star.b[0,1] * mgp.sim[0] + star.b[1,1] * mgp.sim[1] + star.b[2,1] * mgp.sim[2])
        vec3 = (star.b[0,2] * mgp.sim[0] + star.b[1,2] * mgp.sim[1] + star.b[2,2] * mgp.sim[2])

        vrad = vec1
        pmra = vec2 / (K * 10.0)
        pmdec = vec3 / (K * 10.0)
        
        self.result_frame.loc[mgp.name, 'p_rv'] = np.mean(vrad)
        self.result_frame.loc[mgp.name, 'p_erv'] = np.std(vrad,ddof=1)
        self.result_frame.loc[mgp.name, 'p_pmra'] = np.mean(pmra)
        self.result_frame.loc[mgp.name, 'p_epmra'] = np.std(pmra,ddof=1)
        self.result_frame.loc[mgp.name, 'p_pmdec'] = np.mean(pmdec)
        self.result_frame.loc[mgp.name, 'p_epmdec'] = np.std(pmdec,ddof=1)

    # given the ra, dec, and parallax (or expected distance), check the spatial agreement.
    def spatial(self, mgp, star):

        # second: We now need to transform the XYZ position into the coordinate space of the randomly oriented ellipse.
        #  step 1: get the position of the target relative to the moving group
        pos = [star.x-mgp.param3d['X'],star.y-mgp.param3d['Y'],star.z-mgp.param3d['Z']]
        self.result_frame.loc[mgp.name, 'dist'] = np.sqrt(np.mean(pos[0])**2+np.mean(pos[1])**2+np.mean(pos[2])**2)

        #  step 2: rotate that position into the moving group's coordinate system
        rotmatrix = ellipse.rotate(mgp.param3d['XY'],mgp.param3d['XZ'],mgp.param3d['YZ'])
        rotmatrix = rotmatrix.transpose()

        pos = np.dot(rotmatrix,pos)
        # pos now contains the position of the star relative to the center of the moving group, rotated to the frame of the moving group.
        
        # third (and finally): Get a sigma agreement.
        sigmas = np.sqrt((pos[0]/mgp.param3d['D'])**2+(pos[1]/mgp.param3d['E'])**2+(pos[2]/mgp.param3d['F'])**2)
        self.result_frame.loc[mgp.name, 'sig_pos'] = np.mean(sigmas)
        self.result_frame.loc[mgp.name, 'sig_epos'] = np.std(sigmas,ddof=1)

        return self.result_frame.loc[mgp.name, 'sig_pos'], 1

    def pmcompare(self, mgp, star):
        p_pm, p_epm, p_pa, p_epa = astrometry.pmjoin(self.result_frame.loc[mgp.name, "p_pmra"], self.result_frame.loc[mgp.name, "p_epmra"], self.result_frame.loc[mgp.name, "p_pmdec"], self.result_frame.loc[mgp.name, "p_epmdec"])
        # 1a. Calculate the parallel and perpendicular components of the proper motion
        #   First: compute a rotation matrix to transform pmRA,pmDEC into pm= and pm_|_
        cosa = np.cos(p_pa*np.pi/180.0)
        sina = np.sin(p_pa*np.pi/180.0)
        rotmatrix = [ [ cosa, -sina ] , [sina, cosa] ]
        
        #   Second: Apply this rotation matrix to the measured and expected proper 
        #    motion vectors.
        pm_new = np.dot(rotmatrix,[star.pmra,star.pmdec])
        epm_new = np.dot(rotmatrix,[star.epmra,star.epmdec])
        p_pm_new = np.dot(rotmatrix,[self.result_frame.loc[mgp.name, "p_pmra"],self.result_frame.loc[mgp.name, "p_pmdec"]]) # should be [0,pm] if correctly done
        p_epm_new = np.dot(rotmatrix,[self.result_frame.loc[mgp.name, "p_epmra"],self.result_frame.loc[mgp.name, "p_epmdec"]])
        # The 0th element is the perpendicular component, the 1st is parallel
        
        # 1b. Compute a kinematic distance by comparing the magnitude of the 
        #   expected pm vector (which was calculated for a distance of 10 
        #   parsecs) to the measured parallel one
        self.result_frame.loc[mgp.name, "k_dist"] = p_pm_new[1]/pm_new[1] * 10.0
        self.result_frame.loc[mgp.name, "k_edist"] = np.sqrt((p_epm_new[1]/p_pm_new[1])**2 + (epm_new[1]/pm_new[1])**2) * self.result_frame.loc[mgp.name, "k_dist"]

        # 1c. Now we use that kinematic distance to scale the error out. 
        #   Why? If the proper motion was in completely the wrong direction, 
        #   but tiny, it would come out as a match when it shouldn't. 
        #   This scales the error appropriately to stop that from happening.
        p_epm_new = p_epm_new * (10/self.result_frame.loc[mgp.name, "k_edist"])

        # 1d. Compute the proper motion sigma based on the perpendicular 
        #   component (with scaled errors)
        self.result_frame.loc[mgp.name, "sig_pm"] = np.sqrt(pm_new[0]**2)/(np.sqrt(epm_new[0]**2 + p_epm_new[0]**2))

        # 1e. this should weed out anything where the pm vector was flipped (and
        #   distance was negative). I initially thought this should go under the
        #   "distance" section, but I can verify it entirely based on proper
        #   motions even when a distance isn't known, so it should go here.
        if ((abs(p_epm_new[1]) < abs(p_pm_new[1])) & (self.result_frame.loc[mgp.name, "k_dist"] < 0)):
            self.result_frame.loc[mgp.name, "sig_pm"] += 1000.
            # I might also try this relative to a km/sec velocity dispersion, 
            #  pm_new[0]*exp_dist*4.74

        return self.result_frame.loc[mgp.name, "sig_pm"], 1

    def rvcompare(self, mgp, star):
        self.result_frame.loc[mgp.name, "sig_rv"] = np.abs(star.rv - self.result_frame.loc[mgp.name, "p_rv"]) / np.sqrt(star.erv**2 + self.result_frame.loc[mgp.name, "p_erv"]**2)

        return self.result_frame.loc[mgp.name, "sig_rv"], 1

    def distcompare(self, mgp, star):
        self.result_frame.loc[mgp.name, "dist"] = 1/star.plx
        self.result_frame.loc[mgp.name, "edist"] = star.eplx/star.plx**2
        self.result_frame.loc[mgp.name, "sig_dist"] = np.abs(1/star.plx - self.result_frame.loc[mgp.name, "k_dist"]) / np.sqrt((star.eplx/star.plx**2)**2 + self.result_frame.loc[mgp.name, "k_edist"]**2)

        return self.result_frame.loc[mgp.name, "sig_dist"], 1

    def _weightadd(self, sigma, weight):
        sig = 0
        for g in range(len(sigma)):
            sig = sig + weight[g]**2*sigma[g]**2
        sig = np.sqrt(sig)/np.sum(weight)

        return sig

    def _gaus(self, x,a,sigma):
        return a*np.exp(-(x)**2/(2*sigma**2))

    def _gauscdf(self, x,m,sigma):
        return 0.5*(erfc((x-m)/(sigma*np.sqrt(2))))


def csv_loader(filename):
    recognized = []

    datatable = pd.read_csv(filename)
    datatable.columns = map(str.lower, datatable.columns)
    for known in dataformat:
        found = False
        i = 0
        while found == False and i < len(known['keys']):
            key = known['keys'][i]
            if key in datatable:
                recognized.append(known["name"])
                datatable.rename(columns={key: known["name"]}, inplace=True)
                if known['quant'] is not None:
                    datatable[known["name"]] = datatable[known["name"]].apply(lambda x: x*known['quant'])
                found = True
            i += 1

    # if we've gone through all of this and not found an RA or DEC, check for sexagesimal columns
    if "ra" not in recognized:
        if "RAh" in datatable.keys():
            ra = []
            dec = []
            for x in range(len(datatable["RAh"])):
                coord = SkyCoord("{} {} {}".format(datatable['RAh'],datatable['RAm'],datatable['RAs']), "{}{} {} {}".format(datatable['DECs'],datatable['DECd'],datatable['DECm'],datatable['DECs']), units=[units.hourangle, units.degree])
                ra.append(coord.icrs.ra.deg)
                dec.append(coord.icrs.dec.deg)

            datatable["ra"] = ra * units.deg
            datatable["dec"] = dec * units.deg
            recognized.append("rah")
            recognized.append("decd")
    # if we've gone through all of this and not found a parallax, check for distance
    if "plx" not in recognized:
        if "dist" in datatable.keys():
            plx = []
            for x in range(len(datatable["dist"])):
                plx.append(1000/datatable["dist"][x])
            datatable["plx"] = plx
            recognized.append("dist")
    print("Recognized columns:")
    print(", ".join(recognized))
    return datatable


if __name__ == "__main__":
    try:
        infilename = argv[1]
    except IndexError: 
        print('syntax: python lacewing.py inputfile [young] [outfilename] [verbose] [g.o.f.]')
    try:
        if argv[2] == "young":
            young = True
        else:
            young = False
    except IndexError:
        young = False
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

    stardata = csv_loader(infilename)
    print(stardata)

    lacewing = LACEwING('moving_group_all.csv', 10000, young)

    if verbose == "quiet":
        with open(outfilename,'w') as outfile:
            # regular output is a one line per entry summary, also in .csv form
            outfile.write('Name,Note,Group,Probability,Predicted Dist,Predicted Dist uncertainty,Predicted RV,Predicted RV uncertainty,')
            for i in range(len(lacewing.groupdata)):
                outfile.write('{0:},'.format(lacewing.groupdata[i].name))
            outfile.write('\n')
    else:
        if os.path.exists(outfilename):
            os.remove(outfilename)
        


    for cnt in range(len(stardata)):
        star = stardata.iloc[cnt]
        out = lacewing.calculate(star)

        # if verbose output has been selected, print our strings to a file here
        if verbose == "quiet":
            with open(outfilename,'a') as outfile:
                # We've finished processing one star. If regular output was selected, create and print the output string.
                probs = [out['probability'][j] for j in range(len(out))]
                order = np.argsort(probs)[::-1]
                # if even the best match isn't above the threshold of consideration:
                if out['probability'][order[0]] < 20:
                    outfile.write('{0:},{1:},(None),,,,,,'.format(star.name,out['note'][order[0]]))
                else:
                    outfile.write('{0:},{1:},{2:},{3: 5.0f},{4: 8.2f},{5:7.2f},{6:+6.2f},{7: 5.2f},'.format(star.name,out['note'][order[0]],lacewing.mgplist[order[0]],out['probability'][order[0]],out['kin_dist'][order[0]],out['kin_edist'][order[0]],out['kin_rv'][order[0]],out['kin_erv'][order[0]]))
                for k in range(len(out)):
                    outfile.write('{0: 5.0f},'.format(out['probability'][k]))
                outfile.write('\n')
        else:
            out.to_csv(path_or_buf=outfilename, mode="a")
                