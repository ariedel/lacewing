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
