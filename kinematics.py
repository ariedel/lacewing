import numpy
#import galpy
from scipy import weave

# Functions: 

# gal_uvwxyz: converts equatorial to galactic coordinates.
# gal_rdpmuv: converts galactic coordinates into equatorial.
# gal_tester: checks to make sure gal_uvwxyz and gal_rdpmuv are inverses

# ballistic: computes a straight-line linear traceback of equatorial coordinates to a specified time in the past, for a specified number of monte carlo iterations
# epicyclic: the same as ballistic, but uses an epicyclic approximation of galactic orbits from Makarov et al. (2004).
# potential: the same as ballistic, but uses a galactic gravitational potential from galpy

# random: returns uniform random variables from -1 to 1 (rather than 0 to 1)

# ballistic_uniform: the same as ballistic, but for uniformly distributed uncertainties rather than normally distributed uncertainties.
# epicyclic_uniform: uniform monte carlo version of epicyclic 
# potential_uniform: uniform monte carlo version of potential

# from astrolibpy
# AR 2013.0910: Fixed this script to run with vectors.  Mostly, I put back pieces of 
# http://code.google.com/p/astrolibpy/source/browse/astrolib/gal_uvw.py
# from whence the original translation comes.
def gal_uvwxyz(distance=None, lsr=None, ra=None, dec=None, pmra=None, pmdec=None, vrad=None, plx=None):
   """
    NAME:
        GAL_UVWXYZ
    PURPOSE:
        Calculate the Galactic space velocity (U,V,W) and position (X,Y,Z) of star
    EXPLANATION:
        Calculates the Galactic space velocity U, V, W of star given its
        (1) coordinates, (2) proper motion, (3) distance (or parallax), and
        (4) radial velocity.
    CALLING SEQUENCE:
        GAL_UVW [/LSR, RA=, DEC=, PMRA= ,PMDEC=, VRAD= , DISTANCE=
                 PLX= ]
    OUTPUT PARAMETERS:
         U - Velocity (km/s) positive toward the Galactic *anti*center
         V - Velocity (km/s) positive in the direction of Galactic rotation
         W - Velocity (km/s) positive toward the North Galactic Pole
    REQUIRED INPUT KEYWORDS:
         User must supply a position, proper motion,radial velocity and distance
         (or parallax).    Either scalars or vectors can be supplied.
        (1) Position:
         RA - Right Ascension in *Degrees*
         Dec - Declination in *Degrees*
        (2) Proper Motion
         PMRA = Proper motion in RA in arc units (typically milli-arcseconds/yr)
         PMDEC = Proper motion in Declination (typically mas/yr)
        (3) Radial Velocity
         VRAD = radial velocity in km/s
        (4) Distance or Parallax
         DISTANCE - distance in parsecs
                    or
         PLX - parallax with same distance units as proper motion measurements
               typically milliarcseconds (mas)
   
    OPTIONAL INPUT KEYWORD:
         /LSR - If this keyword is set, then the output velocities will be
                corrected for the solar motion (U,V,W)_Sun = (-8.5, 13.38, 6.49)
                (Coskunoglu et al. 2011 MNRAS) to the local standard of rest.
                Note that the value of the solar motion through the LSR remains
                poorly determined.
     EXAMPLE:
         (1) Compute the U,V,W coordinates for the halo star HD 6755.
             Use values from Hipparcos catalog, and correct to the LSR
         ra = ten(1,9,42.3)*15.    & dec = ten(61,32,49.5)
         pmra = 627.89  &  pmdec = 77.84         ;mas/yr
         dis = 144    &  vrad = -321.4
         gal_uvw,u,v,w,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,vrad=vrad,dis=dis,/lsr
             ===>  u=154  v = -493  w = 97        ;km/s
   
         (2) Use the Hipparcos Input and Output Catalog IDL databases (see
         http://idlastro.gsfc.nasa.gov/ftp/zdbase/) to obtain space velocities
         for all stars within 10 pc with radial velocities > 10 km/s
   
         dbopen,'hipparcos,hic'      ;Need Hipparcos output and input catalogs
         list = dbfind('plx>100,vrad>10')      ;Plx > 100 mas, Vrad > 10 km/s
         dbext,list,'pmra,pmdec,vrad,ra,dec,plx',pmra,pmdec,vrad,ra,dec,plx
         ra = ra*15.                 ;Need right ascension in degrees
         GAL_UVW,u,v,w,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,vrad=vrad,plx = plx
         forprint,u,v,w              ;Display results
    METHOD:
         Follows the general outline of Johnson & Soderblom (1987, AJ, 93,864)
         except that U is positive outward toward the Galactic *anti*center, and
         the J2000 transformation matrix to Galactic coordinates is taken from
         the introduction to the Hipparcos catalog.
    REVISION HISTORY:
         Written, W. Landsman                       December   2000
         fix the bug occuring if the input arrays are longer than 32767
           and update the Sun velocity           Sergey Koposov June 2008
           vectorization of the loop -- performance on large arrays
           is now 10 times higher                Sergey Koposov December 2008
   """

   n_params = 3
   
   if n_params == 0:  
      print 'Syntax - GAL_UVW, U, V, W, [/LSR, RA=, DEC=, PMRA= ,PMDEC=, VRAD='
      print '                  Distance=, PLX='
      print '         U, V, W, X, Y, Z - output Galactic space velocities (km/s) and positions'
      return None
   
   if ra is None or dec is None:  
      raise Exception('ERROR - The RA, Dec (J2000) position keywords must be supplied (degrees)')
   if plx is None and distance is None:
      raise Exception('ERROR - Either a parallax or distance must be specified')
   if distance is not None:
      if numpy.any(distance==0):
         raise Exception('ERROR - All distances must be > 0')
      plx = 1 / distance          #Parallax in arcseconds
   if plx is not None and numpy.any(plx==0):
      raise Exception('ERROR - Parallaxes must be > 0')
   
   cosd = numpy.cos(dec*numpy.pi/180.)
   sind = numpy.sin(dec*numpy.pi/180.)
   cosa = numpy.cos(ra*numpy.pi/180.)
   sina = numpy.sin(ra*numpy.pi/180.)
   
   k = 4.74047     #Equivalent of 1 A.U/yr in km/s  
   a_g = numpy.array([[-0.0548755604, +0.4941094279, -0.8676661490],
                [-0.8734370902, -0.4448296300, -0.1980763734],
                [-0.4838350155, 0.7469822445, +0.4559837762]]) # rotation matrix for J2000 -->Galactic
   pos1 = cosd*cosa
   pos2 = cosd*sina
   pos3 = sind

   #AR 2013.0910: In order to use this with vectors, we need more control over the matrix multiplication

   x = 1/plx * (a_g[0,0] * pos1 + a_g[1,0] * pos2 + a_g[2,0] * pos3)
   y = 1/plx * (a_g[0,1] * pos1 + a_g[1,1] * pos2 + a_g[2,1] * pos3)
   z = 1/plx * (a_g[0,2] * pos1 + a_g[1,2] * pos2 + a_g[2,2] * pos3)

   vec1 = vrad
   vec2 = k * pmra / plx
   vec3 = k * pmdec / plx

   #AR 2013.0910: In order to use this with vectors, we need more control over the matrix multiplication
   u = (a_g[0,0] * cosa * cosd + a_g[1,0] * sina * cosd + a_g[2,0] * sind) * vec1 + (-a_g[0,0] * sina + a_g[1,0] * cosa) * vec2 + (-a_g[0,0] * cosa * sind - a_g[1,0] * sina * sind + a_g[2,0] * cosd) * vec3
   v = (a_g[0,1] * cosa * cosd + a_g[1,1] * sina * cosd + a_g[2,1] * sind) * vec1 + (-a_g[0,1] * sina + a_g[1,1] * cosa) * vec2 + (-a_g[0,1] * cosa * sind - a_g[1,1] * sina * sind + a_g[2,1] * cosd) * vec3
   w = (a_g[0,2] * cosa * cosd + a_g[1,2] * sina * cosd + a_g[2,2] * sind) * vec1 + (-a_g[0,2] * sina + a_g[1,2] * cosa) * vec2 + (-a_g[0,2] * cosa * sind - a_g[1,2] * sina * sind + a_g[2,2] * cosd) * vec3

   lsr_vel = numpy.array([8.5, 13.38, 6.49])
   if (lsr is not None):  
      u = u + lsr_vel[0]
      v = v + lsr_vel[1]
      w = w + lsr_vel[2]
   
   return (u,v,w,x,y,z)

def gal_rdp(u,v,w,x,y,z):

   # conversion between galactic and equatorial coordinates
   k = 4.74047     #Equivalent of 1 A.U/yr in km/s  
   a_g = numpy.array([[-0.0548755604, -0.8734370902,-0.4838350155],
                [+0.4941094279, -0.4448296300, 0.7469822445],
                [-0.8676661490, -0.1980763734, +0.4559837762]]) # rotation matrix for Galactic-->J2000
   radeg = 180/numpy.pi

   #AR 2014.0123: First, rotate galactic xyz back to equatorial xyz
   pos1 = (a_g[0,0] * x + a_g[1,0] * y + a_g[2,0] * z)
   pos2 = (a_g[0,1] * x + a_g[1,1] * y + a_g[2,1] * z)
   pos3 = (a_g[0,2] * x + a_g[1,2] * y + a_g[2,2] * z)

   dist = numpy.sqrt(pos1**2 + pos2**2 + pos3**2)
   ra = numpy.arctan2(pos2,pos1)*radeg
   try:
      flip = numpy.where(ra < 0)
      ra[flip] = 360 + ra[flip]
   except IndexError:
      if ra < 0:
         ra = 360 + ra
   dec = 90.0 - numpy.arccos(pos3/dist)*radeg

   cosd = numpy.cos(dec/radeg)
   sind = numpy.sin(dec/radeg)
   cosa = numpy.cos(ra/radeg)
   sina = numpy.sin(ra/radeg)
    
   a_c = numpy.array([ [cosa*cosd,-sina,-cosa*sind],
                       [sina*cosd,cosa,-sina*sind],
                       [sind,0,cosd] ]) # rotation matrix for cartesian to spherical

   b = numpy.dot(a_g,a_c)
   #vec = numpy.dot(vec,b)

   vec1 = (b[0,0] * u + b[1,0] * v + b[2,0] * w)
   vec2 = (b[0,1] * u + b[1,1] * v + b[2,1] * w)
   vec3 = (b[0,2] * u + b[1,2] * v + b[2,2] * w)

   vrad = vec1
   pmra = vec2 / (k * dist)
   pmdec = vec3 / (k * dist)

   return ra,dec,dist,pmra,pmdec,vrad

def gal_tester():
   for i in xrange(1000):
      ra1 = numpy.random.rand() * 360
      dec1 = (numpy.random.rand()-0.5) * 180
      pi1 = numpy.random.rand() * 100
      pmra1 = (numpy.random.rand()-0.5) * 20
      pmdec1 = (numpy.random.rand()-0.5) * 20
      vrad1 = (numpy.random.rand()-0.5) * 200

      u,v,w,x,y,z = gal_uvwxyz(distance=pi1,ra=ra1,dec=dec1,pmra=pmra1,pmdec=pmdec1,vrad=vrad1)
      ra2,dec2,pi2,pmra2,pmdec2,vrad2 = gal_rdp(u,v,w,x,y,z)

      print ra1,dec1,1/pi1,pmra1,pmdec1,vrad1
      print ra2,dec2,pi2,pmra2,pmdec2,vrad2
      print ra1-ra2,dec1-dec2,1/pi1-pi2,pmra1-pmra2,pmdec1-pmdec2,vrad1-vrad2
      print 

def ballistic(ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv,timespan,timestep,n_int):
   n_time = timespan/timestep
   times = numpy.arange(0, timespan,timestep)

   px = numpy.zeros((n_int,n_time))
   py = numpy.zeros((n_int,n_time))
   pz = numpy.zeros((n_int,n_time))
   pc = []
   
   for j in xrange(n_int):
      bra = ra + numpy.random.randn()*era
      bdec = dec + numpy.random.randn()*edec
      bdist = dist + numpy.random.randn()*edist
      bpmra = pmra + numpy.random.randn()*epmra
      bpmdec = pmdec + numpy.random.randn()*epmdec
      brv = rv + numpy.random.randn()*erv
      tu,tv,tw,tx,ty,tz = gal_uvwxyz(ra=bra,dec=bdec,distance=bdist,pmra=bpmra,pmdec=bpmdec,vrad=brv)
      
      # backtrack this one simulation through all of time
      px[j] = tx + tu * times * 1.0226
      py[j] = ty + tv * times * 1.0226
      pz[j] = tz + tw * times * 1.0226

   return px,py,pz

def epicyclic(ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv,timespan,timestep,n_int):
   n_time = timespan/timestep
   times = numpy.arange(0,timespan,timestep)
   px = numpy.zeros((n_int,n_time))
   py = numpy.zeros((n_int,n_time))
   pz = numpy.zeros((n_int,n_time))
   pc = []

   # Oort constants from Bobylev et al. 2010 2010MNRAS.408.1788B
   A = 0.0178
   B = -0.0132
   vos = 2*numpy.pi/85 # Vertical Oscillations
   kap = numpy.sqrt(-4*B*(A-B)) # planar epicyclic velocity
   
   for j in xrange(n_int):
      bra = ra + (numpy.random.randn()*era)*numpy.cos(dec*numpy.pi/180.)
      bdec = dec + numpy.random.randn()*edec
      bdist = dist + numpy.random.randn()*edist
      bpmra = pmra + numpy.random.randn()*epmra
      bpmdec = pmdec + numpy.random.randn()*epmdec
      brv = rv + numpy.random.randn()*erv
      tu,tv,tw,tx,ty,tz = gal_uvwxyz(ra=bra,dec=bdec,distance=bdist,pmra=bpmra,pmdec=bpmdec,vrad=brv)
      
      # backtrack this one simulation through all of time
      # Calculations taken from Section 4 of Makarov, Olling, and Teuben (2004) 2004MNRAS.352.1199M
      # AR 2013.1110 I'm multiplying tx and ty by ones because Python is having trouble broadcasting.
      px[j] = tx + tu*kap**(-1)*numpy.sin(kap*times) + (tv - 2*A*tx)*(1-numpy.cos(kap*times))*(2*B)**(-1)
      py[j] = ty - tu*(1-numpy.cos(kap*times))*(2*B)**(-1) + tv*(A*kap*times - (A - B)*numpy.sin(kap*times))*(kap*B)**(-1) - tx*2*A*(A-B)*(kap*times - numpy.sin(kap*times))*(kap*B)**(-1)
      pz[j] = tz*numpy.cos(vos*times) + tw*(vos)**(-1)*numpy.sin(vos*times)


   return px,py,pz

#def potential(ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv,timespan,timestep,n_int):
#   galacticpotential = galpy.potential.MWPotential2014()
#   times = numpy.arange(0,timespan,timestep)
#
#   for j in xrange(n_int):
#      bra = ra + numpy.random.randn()*era
#      bdec = dec + numpy.random.randn()*edec
#      bdist = dist + numpy.random.randn()*edist
#      bpmra = pmra + numpy.random.randn()*epmra
#      bpmdec = pmdec + numpy.random.randn()*epmdec
#      brv = rv + numpy.random.randn()*erv
#      
#      orbit = galpy.orbit.Orbit(vxvv=[bra,bdec,bdist/1000,bpmra,bpmdec,brv],radec=True)
#      # integrate the orbit
#      orbit.integrate(times,galacticpotential)
#      # get the orbit
#      orbitresult = orbit.getorbit()
#
#   return px,py,pz

def random():
    return (numpy.random.rand()*2)-1

def ballistic_uniform(ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv,timespan,timestep,n_int):
   n_time = timespan/timestep
   times = numpy.arange(0, timespan,timestep)

   px = numpy.zeros((n_int,n_time))
   py = numpy.zeros((n_int,n_time))
   pz = numpy.zeros((n_int,n_time))
   pc = []
   
   for j in xrange(n_int):
      bra = ra + random()*era
      bdec = dec + random()*edec
      bdist = dist + random()*edist
      bpmra = pmra + random()*epmra
      bpmdec = pmdec + random()*epmdec
      brv = rv + random()*erv
      tu,tv,tw,tx,ty,tz = gal_uvwxyz(ra=bra,dec=bdec,distance=bdist,pmra=bpmra,pmdec=bpmdec,vrad=brv)
      
      # backtrack this one simulation through all of time
      px[j] = tx + tu * times * 1.0226
      py[j] = ty + tv * times * 1.0226
      pz[j] = tz + tw * times * 1.0226

   return px,py,pz

def epicyclic_uniform(ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv,timespan,timestep,n_int):
   n_time = timespan/timestep
   times = numpy.arange(0,timespan,timestep)
   px = numpy.zeros((n_int,n_time))
   py = numpy.zeros((n_int,n_time))
   pz = numpy.zeros((n_int,n_time))
   pc = []

   # Oort constants from Bobylev et al. 2010 2010MNRAS.408.1788B
   A = 0.0178
   B = -0.0132
   vos = 2*numpy.pi/85 # Vertical Oscillations
   kap = numpy.sqrt(-4*B*(A-B)) # planar epicyclic velocity
   
   for j in xrange(n_int):
      bra = ra + random()*era
      bdec = dec + random()*edec
      bdist = dist + random()*edist
      bpmra = pmra + random()*epmra
      bpmdec = pmdec + random()*epmdec
      brv = rv + random()*erv
      tu,tv,tw,tx,ty,tz = gal_uvwxyz(ra=bra,dec=bdec,distance=bdist,pmra=bpmra,pmdec=bpmdec,vrad=brv)
      
      # backtrack this one simulation through all of time
      # Calculations taken from Section 4 of Makarov, Olling, and Teuben (2004) 2004MNRAS.352.1199M
      # AR 2013.1110 I'm multiplying tx and ty by ones because Python is having trouble broadcasting.
      px[j] = tx + tu*kap**(-1)*numpy.sin(kap*times) + (tv - 2*A*tx)*(1-numpy.cos(kap*times))*(2*B)**(-1)
      py[j] = ty - tu*(1-numpy.cos(kap*times))*(2*B)**(-1) + tv*(A*kap*times - (A - B)*numpy.sin(kap*times))*(kap*B)**(-1) - tx*2*A*(A-B)*(kap*times - numpy.sin(kap*times))*(kap*B)**(-1)
      pz[j] = tz*numpy.cos(vos*times) + tw*(vos)**(-1)*numpy.sin(vos*times)


   return px,py,pz

#def potential_uniform(ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv,timespan,timestep,n_int):
#   galacticpotential = galpy.potential.MWPotential2014()
#   times = numpy.arange(0,timespan,timestep)
#
#   for j in xrange(n_int):
#      bra = ra + random()*era
#      bdec = dec + random()*edec
#      bdist = dist + random()*edist
#      bpmra = pmra + random()*epmra
#      bpmdec = pmdec + random()*epmdec
#      brv = rv + random()*erv
#      
#      orbit = galpy.orbit.Orbit(vxvv=[bra,bdec,bdist/1000,bpmra,bpmdec,brv],radec=True)
#      # integrate the orbit
#      orbit.integrate(times,galacticpotential)
#      # get the orbit
#      orbitresult = orbit.getorbit()
#
#   return px,py,pz

## from astrolibpy
## AR 2013.0910: Fixed this script to run with vectors.  Mostly, I put back pieces of 
## http://code.google.com/p/astrolibpy/source/browse/astrolib/gal_uvw.py
## from whence the original translation comes.
#def gal_uvwxyz_weave(distance=None, lsr=None, ra=None, dec=None, pmra=None, pmdec=None, vrad=None, plx=None):
#    n_params = 3
#   
#    if n_params == 0:  
#        print 'Syntax - GAL_UVW, U, V, W, [/LSR, RA=, DEC=, PMRA= ,PMDEC=, VRAD='
#        print '                  Distance=, PLX='
#        print '         U, V, W, X, Y, Z - output Galactic space velocities (km/s) and positions'
#        return None
#   
#    if ra is None or dec is None:  
#        raise Exception('ERROR - The RA, Dec (J2000) position keywords must be supplied (degrees)')
#    if plx is None and distance is None:
#        raise Exception('ERROR - Either a parallax or distance must be specified')
#    if distance is not None:
#        if numpy.any(distance==0):
#            raise Exception('ERROR - All distances must be > 0')
#        plx = 1 / distance          #Parallax in arcseconds
#    if plx is not None and numpy.any(plx==0):
#        raise Exception('ERROR - Parallaxes must be > 0')
#
#    length = len(ra)
#    u = numpy.zeros(length)
#    v = numpy.zeros(length)
#    w = numpy.zeros(length)
#    x = numpy.zeros(length)
#    y = numpy.zeros(length)
#    z = numpy.zeros(length)
#    uvwxyz = numpy.zeros((6,length))
#
#    support = "#include <math.h>"
#    code = """
#   
#    double k = 4.74047;
#    double a_g00 = -0.0548755604;
#    double a_g01 = +0.4941094279;
#    double a_g02 = -0.8676661490;
#    double a_g10 = -0.8734370902;
#    double a_g11 = -0.4448296300;
#    double a_g12 = -0.1980763734;
#    double a_g20 = -0.4838350155;
#    double a_g21 = +0.7469822445;
#    double a_g22 = +0.4559837762;
#
#    //PyArrayObject* uvwxyz((6,length));
#    //double uvwxyz [6][length];
#
#    for(int i = 0; i < length; i++){
#       double cosd = cos(dec[i]*3.1415926535/180.0);
#       double sind = sin(dec[i]*3.1415926535/180.0);
#       double cosa = cos(ra[i]*3.1415926535/180.0);
#       double sina = sin(ra[i]*3.1415926535/180.0);
#
#       double pos1 = cosd*cosa;
#       double pos2 = cosd*sina;
#       double pos3 = sind;
#
#       uvwxyz[3][i] = 1/plx[i] * (a_g00 * pos1 + a_g10 * pos2 + a_g20 * pos3);
#       uvwxyz(4,i) = 1/plx[i] * (a_g01 * pos1 + a_g11 * pos2 + a_g21 * pos3);
#       uvwxyz(5,i) = 1/plx[i] * (a_g02 * pos1 + a_g12 * pos2 + a_g22 * pos3);
#
#       double vec1 = vrad[i];
#       double vec2 = k * pmra[i] / plx[i];
#       double vec3 = k * pmdec[i] / plx[i];
#
#       uvwxyz[0][i] = (a_g00 * cosa * cosd + a_g10 * sina * cosd + a_g20 * sind) * vec1 + (-a_g00 * sina + a_g10 * cosa) * vec2 + (-a_g00 * cosa * sind - a_g10 * sina * sind + a_g20 * cosd) * vec3;
#       uvwxyz[1][i] = (a_g01 * cosa * cosd + a_g11 * sina * cosd + a_g21 * sind) * vec1 + (-a_g01 * sina + a_g11 * cosa) * vec2 + (-a_g01 * cosa * sind - a_g11 * sina * sind + a_g21 * cosd) * vec3;
#       uvwxyz[2][i] = (a_g02 * cosa * cosd + a_g12 * sina * cosd + a_g22 * sind) * vec1 + (-a_g02 * sina + a_g12 * cosa) * vec2 + (-a_g02 * cosa * sind - a_g12 * sina * sind + a_g22 * cosd) * vec3;
#
#    return_val = uvwxyz;
#    }
#
#    
#    """
#    result = weave.inline(code, ['ra','dec','plx','pmra','pmdec','vrad','length','uvwxyz'], support_code = support, libraries = ['m'],type_factories = blitz_type_factories)
#    print result
#
#    lsr_vel = numpy.array([8.5, 13.38, 6.49])
#    if (lsr is not None):  
#        u = u + lsr_vel[0]
#        v = v + lsr_vel[1]
#        w = w + lsr_vel[2]
#   
#    return (u,v,w,x,y,z)

