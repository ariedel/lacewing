import numpy
from scipy.special import erfc
import ellipse
import kinematics

# Calculate the proper motion vector at the star's RA and DEC and a distance of 10 pc
def converge(assoc,ra,era,dec,edec,n_init):

   ra_tmp = ra + numpy.random.randn(n_init) * era / numpy.cos(dec*180/numpy.pi)
   dec_tmp = dec + numpy.random.randn(n_init) * edec

   # generate random points within the moving group distribution
   a = assoc.A * numpy.random.randn(n_init)
   b = assoc.B * numpy.random.randn(n_init)
   c = assoc.C * numpy.random.randn(n_init)
   # rotate them to account for the skewed distribution
   rotmatrix = ellipse.rotate(assoc.UV,assoc.UW,assoc.VW)
   rotmatrix = rotmatrix.transpose()
   vec = numpy.dot(rotmatrix,numpy.asarray([a,b,c]))

   u = vec[0]+assoc.U
   v = vec[1]+assoc.V
   w = vec[2]+assoc.W

   # conversion between galactic and equatorial coordinates
   k = 4.74047     #Equivalent of 1 A.U/yr in km/s  
   a_g = numpy.array([[-0.0548755604, -0.8734370902,-0.4838350155],
                [+0.4941094279, -0.4448296300, 0.7469822445],
                [-0.8676661490, -0.1980763734, +0.4559837762]]) # rotation matrix for Galactic-->J2000
   radeg = 180/numpy.pi

   cosd = numpy.cos(dec_tmp/radeg)
   sind = numpy.sin(dec_tmp/radeg)
   cosa = numpy.cos(ra_tmp/radeg)
   sina = numpy.sin(ra_tmp/radeg)
    
   a_c = numpy.array([ [cosa*cosd,-sina,-cosa*sind],
                       [sina*cosd,cosa,-sina*sind],
                       [sind,0,cosd] ]) # rotation matrix for cartesian to spherical

   b = numpy.dot(a_g,a_c)

   vec1 = (b[0,0] * u + b[1,0] * v + b[2,0] * w)
   vec2 = (b[0,1] * u + b[1,1] * v + b[2,1] * w)
   vec3 = (b[0,2] * u + b[1,2] * v + b[2,2] * w)

   vrad = vec1
   pmra = vec2 / (k * 10.0)
   pmdec = vec3 / (k * 10.0)
   
   exp_rv = numpy.mean(vrad)
   exp_erv = numpy.std(vrad,ddof=1)
   exp_pmra = numpy.mean(pmra)
   exp_epmra = numpy.std(pmra,ddof=1)
   exp_pmdec = numpy.mean(pmdec)
   exp_epmdec = numpy.std(pmdec,ddof=1)

   return exp_pmra,exp_epmra,exp_pmdec,exp_epmdec,exp_rv,exp_erv


# given the ra, dec, and parallax (or expected distance), check the spatial agreement.
def spatial(mgp,ra,era,dec,edec,pi,epi,n_init):
   
   # first, get the XYZ position
   tra = ra + (numpy.random.randn(n_init)*era)*numpy.cos(dec*numpy.pi/180.)
   tdec = dec + numpy.random.randn(n_init)*edec
   tpi = pi + numpy.random.randn(n_init)*epi
   u,v,w,x,y,z = kinematics.gal_uvwxyz(ra=tra,dec=tdec,plx=tpi,pmra=numpy.zeros_like(tra),pmdec=numpy.zeros_like(tra),vrad=numpy.zeros_like(tra))

   # second: We now need to transform the XYZ position into the coordinate space of the randomly oriented ellipse.
   #  step 1: get the position of the target relative to the moving group
   star = [x-mgp.X,y-mgp.Y,z-mgp.Z]

   #  step 2: rotate that position into the moving group's coordinate system
   rotmatrix = ellipse.rotate(mgp.XY,mgp.XZ,mgp.YZ)
   rotmatrix = rotmatrix.transpose()

   pos = numpy.dot(rotmatrix,star)
   # pos now contains the position of the star relative to the center of the moving group, rotated to the frame of the moving group.
   
   # third (and finally): Get a sigma agreement.
   sigmas = numpy.sqrt((pos[0]/mgp.D)**2+(pos[1]/mgp.E)**2+(pos[2]/mgp.F)**2)
   sigma = numpy.mean(sigmas)
   esigma = numpy.std(sigmas,ddof=1)

   return sigma,esigma

##############################################
# PROPER MOTION to GREAT CIRCLE
##############################################

def greatcircle(ra,dec,pmra,pmdec):
   radeg=180/numpy.pi
   # set up a great circle
   raeq=numpy.arange(0,361,1)/radeg
   deceq=numpy.zeros(361)
                                # to find the pole of the appropriate
                                # great circle, we must take the cross
                                # product of the two unit vectors
                                # (ra,dec) and (ra+pmra/dec,dec+pmdec)
                                # in cartesian coordinates.
   rar=ra/radeg
   decr=dec/radeg
#    print ra,dec,rar,decr
   pmrar=pmra/3600./radeg
   pmdecr=pmdec/3600./radeg
   x1=numpy.cos(decr)*numpy.cos(rar)
   y1=numpy.cos(decr)*numpy.sin(rar)
   z1=numpy.sin(decr)
   x2=numpy.cos(decr+pmdecr)*numpy.cos(rar+pmrar/numpy.cos(decr))
   y2=numpy.cos(decr+pmdecr)*numpy.sin(rar+pmrar/numpy.cos(decr))
   z2=numpy.sin(decr+pmdecr)
   #    print x1,y1,z1,sqrt(x1**2+y1**2+z1**2),x2,y2,z2,sqrt(x2**2+y2**2+z2**2)
   #   print

    # the cross-product:
   x=y1*z2-z1*y2
   y=z1*x2-x1*z2
   z=x1*y2-y1*x2
    # (x,y,z) is the coordinate of the great circle's pole.
   r=numpy.sqrt(x**2+y**2+z**2)

    # get the RA and DEC (in radians, fixing the bounds)
   rap=numpy.arctan2(y,x)
   if rap < 0:
      rap = rap + 2*numpy.pi
   decp=numpy.arcsin(z/r)
#   if decp gt !pi/2. then decp=!pi/2.-decp
   
#   print astrotools.deg2HMS(ra=ra),astrotools.deg2HMS(dec=dec),astrotools.deg2HMS(ra=(ra+pmra/3600./numpy.cos(dec))),astrotools.deg2HMS(dec=dec+pmdec/3600.),astrotools.deg2HMS(ra=rap*radeg),astrotools.deg2HMS(dec=decp*radeg)
    # angular distance between star and equatorial pole (ie, the dec)
   dlon=(90-dec)/radeg
   sdp = numpy.sin(decp)
   cdp = numpy.sqrt(1.0-sdp*sdp)
   
    # stolen from glactc.pro
   sgb = numpy.sin(deceq)
   cgb = numpy.sqrt(1.0-sgb*sgb)
   sdec = sgb*sdp + cgb*cdp*numpy.cos(dlon-raeq)
   decgc = radeg * numpy.arcsin(sdec)
   cdec = numpy.sqrt(1.0-sdec*sdec)
   sinf = cgb * numpy.sin(dlon-raeq) / cdec
   cosf = (sgb-sdp*sdec) / (cdp*cdec)
   ragc = radeg * (rap + numpy.arctan2(sinf,cosf))
   
   ragc[numpy.where(ragc < 0)] = ragc[numpy.where(ragc < 0)]+360

   return ragc,decgc

def weightadd(sigma):
   sig = 0
   for g in range(len(sigma)):
      sig = sig + sigma[g]**2
   sig = numpy.sqrt(sig)

   if len(sigma) > 1:
      sig = sig/(len(sigma))

   return sig

def gaus(x,a,sigma):
    return a*numpy.exp(-(x)**2/(2*sigma**2))

def gauscdf(x,m,sigma):
    return 0.5*(erfc((x-m)/(sigma*numpy.sqrt(2))))
