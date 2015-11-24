import numpy

# Weighted Standard Devation taken from http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
def wstddev(x,w):
    average = numpy.average(x,weights=w)
    variance = numpy.dot(w,(x-average)**2)/w.sum()
    return average,numpy.sqrt(variance)

def isnumber(x):
    #print type(x).__name__
    result = False
    if type(x).__name__ in ["str","string_"]:
        try:
            float(x)/2.0
            result = True
        except ValueError:
            result = False
    if type(x).__name__ in ["int", "int64", "long", "float", "float64", "double"]:
        return True
        
    if type(x).__name__ in ["list","ndarray","Column","MaskedColumn"]:
        result=[]
        for i in xrange(len(x)):
            #print x[i]
            try:
                if isinstance(x[i],numpy.ma.core.MaskedConstant):
                    raise ValueError
                float(x[i])/2. # should induce non-numbers to throw a 'ValueError'
                result.append(True)
            except ValueError,TypeError:
                result.append(False)
    return result


def meanpi(pi,epi,fpi):
    # as seen on www.denseproject.com/25pc/
    npi = len(pi)
    
    top = 0
    bottom = 0
    number = 0
    for i in range(npi):
        if fpi[i] == 'P':
            top = top + (pi[i] / (epi[i]**2))
            bottom = bottom + (1.0 / (epi[i]**2))
            number = number + 1

    if number > 0:
        wpi = top/bottom
        ewpi = 1/numpy.sqrt(bottom)
    else:
        wpi = 0
        ewpi = 0

    return wpi,ewpi

# Astrometry routines
def sixty(x):
    hour = int(x)
    minute = int((x - hour)*60.0)
    second = (((x - hour)*60.0)-minute)*60.0
    return hour,abs(minute),abs(second)

def ten(x):
    #x = [x]
    size = len(x)
    #print size,x
    try:
        len(x[0])
        x = numpy.asarray(x)
        val = numpy.zeros(len(x[0]))
    except TypeError:
        v = 0
                      
    if size == 1:
        val = x[0]
    elif size == 2:
        val = abs(x[0]) + x[1]/60.
    elif size == 3:
        val = abs(x[0]) + x[1]/60. + x[2]/3600.

    try:
        len(x[0]) # will raise TYPEERROR if this is not an array
        flip = numpy.where(x[0] < 0)
        val[flip] = val[flip] * -1
    except TypeError:
        if x[0] < 0:
            val = val * -1

    #print val
    return val

###################################
# PM join
###################################

def pmjoin(pmra,epmra,pmdec,epmdec):
    pm=numpy.sqrt(pmra**2+pmdec**2)
    epm = numpy.sqrt(epmra**2 + epmdec**2)
    
    pa = numpy.arctan2(pmra,pmdec)*180./numpy.pi
    if pa < 0:
        pa = 360 + pa
        
    epa = abs(pa - numpy.arctan2(pmra-epmra,pmdec-epmdec)*180./numpy.pi)
    if epa > 180:
        epa = 360-epa

         
    return pm,epm,pa,epa

###################################
# PM split
###################################
def pmsplit(pm,epm,pa,epa):
   pa = (pa % 360) * numpy.pi/180.
   pi2 = numpy.pi / 2.

   if ((pa >= 0.0) and (pa < pi2)):
      pmra = pm * numpy.sin(pa)
      pmdec = pm * numpy.cos(pa)

   if ((pa >= pi2) and (pa < 2*pi2)):
      pa = pa - pi2
      pmra = pm * numpy.cos(pa)
      pmdec = - (pm * numpy.sin(pa))
      
   if (pa >= 2*pi2) and (pa < 3*pi2):
      pa = pa - (2*pi2)
      pmra =  - (pm * numpy.sin(pa)) 
      pmdec = - (pm * numpy.cos(pa))

   if (pa >= 3*pi2) and (pa < 4*pi2):
      pa = pa - (3*pi2)
      pmra =  - (pm * numpy.cos(pa))
      pmdec = pm * numpy.sin(pa)

   return pmra,numpy.sqrt(epm),pmdec,numpy.sqrt(epm)

def addpm(ra,era,dec,edec,pmra,epmra,pmdec,epmdec,time):
    ra = ra + pmra*numpy.cos(dec*numpy.pi/180)/3600. * time/365.25
    dec = dec + pmdec/3600. * time/365.25

    era = numpy.sqrt((era**2) + (epmra*numpy.cos(dec*numpy.pi/180)/3600. * time/365.25)**2)
    edec = numpy.sqrt((edec**2) + (epmdec/3600. * time/365.25)**2)


    return ra,era,dec,edec

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

