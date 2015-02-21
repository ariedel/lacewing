import numpy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot
#from matplotlib import cm
from matplotlib.patches import Ellipse
import kinematics
import ellipse
import sys
from astropy.io import ascii

def traceback(argv=None):
    if argv is None:
        argv = sys.argv
        
    readtable = ascii.get_reader(Reader=ascii.Basic)
    readtable.header.splitter.delimiter = ','
    readtable.data.splitter.delimiter = ','
    readtable.header.start_line = 1
    readtable.data.start_line = 3

    young = readtable.read(argv[1])
    mgpname = argv[2]
    n_init = int(argv[3])

    # How many stars are we fitting? Find all that are bona fide members (quality is "Good")
    n_stars=0

    n_stars = numpy.where(numpy.bitwise_and(kinematics.isnumber(young['rv']),kinematics.isnumber(young['pi'])) & (young['GROUP'] == mgpname) & (young['GROUP quality'] == "Good"))[0]

    print n_stars

    ra = []
    era = []
    dec = []
    edec = []
    pi = []
    epi = []
    pmra = []
    epmra = []
    pmdec = []
    epmdec = []
    rv = []
    erv = []

    #print young

    for i in n_stars:
        ra.append(float(young[i]['RA']))
        era.append(float(young[i]['eRA'])/3600000)
        dec.append(float(young[i]['DEC']))
        edec.append(float(young[i]['eDEC'])/3600000)
        pi.append(float(young[i]['pi'])/1000.)
        epi.append(float(young[i]['epi'])/1000.)
        pmra.append(float(young[i]['pmRA'])/1000.)
        epmra.append(float(young[i]['epmRA'])/1000.)
        pmdec.append(float(young[i]['pmDEC'])/1000.)
        epmdec.append(float(young[i]['epmDEC'])/1000.)
        rv.append(float(young[i]['rv']))
        erv.append(float(young[i]['erv']))
        print len(ra), young[i]['Name'],ra[-1],era[-1],dec[-1],edec[-1],pi[-1],epi[-1],pmra[-1],epmra[-1],pmdec[-1],epmdec[-1],rv[-1],erv[-1]
    ra = numpy.asarray(ra)
    era = numpy.asarray(era)
    dec = numpy.asarray(dec)
    edec = numpy.asarray(edec)
    pi = numpy.asarray(pi)
    epi = numpy.asarray(epi)
    pmra = numpy.asarray(pmra)
    epmra = numpy.asarray(epmra)
    pmdec = numpy.asarray(pmdec)
    epmdec = numpy.asarray(epmdec)
    rv = numpy.asarray(rv)
    erv = numpy.asarray(erv)
    print len(ra),ra

    uvwlist = []
    xyzlist = []
    outfile = open("Moving_Group_{0:}.dat".format(mgpname),"wb")
    outfile.write("Time X eX Y eY Z eZ A eA B eB C eC XY eXY XZ eXZ YZ eYZ\n")
    for j in xrange(n_init):
        tra = ra + (numpy.random.randn()*era)*numpy.cos(dec*numpy.pi/180.)
        tdec = dec + numpy.random.randn()*edec
        tpi = pi + numpy.random.randn()*epi
        tpmra = pmra + numpy.random.randn()*epmra 
        tpmdec = pmdec + numpy.random.randn()*epmdec 
        trv = rv + numpy.random.randn()*erv
        #print len(tra),tra
        tu,tv,tw,tx,ty,tz = kinematics.gal_uvwxyz(ra=tra,dec=tdec,plx=tpi,pmra=tpmra,pmdec=tpmdec,vrad=trv)
        obj = ellipse.fitellipse(tx,ty,tz)
        xyzlist.append(obj)
        obj = ellipse.fitellipse(tu,tv,tw)
        uvwlist.append(obj)

    u = numpy.mean([uvwlist[m]['x'] for m in range(n_init)])
    v = numpy.mean([uvwlist[m]['y'] for m in range(n_init)])
    w = numpy.mean([uvwlist[m]['z'] for m in range(n_init)])
    uv = numpy.mean([uvwlist[m]['xy'] for m in range(n_init)])
    uw = numpy.mean([uvwlist[m]['xz'] for m in range(n_init)])
    vw= numpy.mean([uvwlist[m]['yz'] for m in range(n_init)])
    a = numpy.mean([uvwlist[m]['a'] for m in range(n_init)])
    b = numpy.mean([uvwlist[m]['b'] for m in range(n_init)])
    c = numpy.mean([uvwlist[m]['c'] for m in range(n_init)])
    eu = numpy.std([uvwlist[m]['x'] for m in range(n_init)],ddof=1)
    ev = numpy.std([uvwlist[m]['y'] for m in range(n_init)],ddof=1)
    ew = numpy.std([uvwlist[m]['z'] for m in range(n_init)],ddof=1)
    euv = numpy.std([uvwlist[m]['xy'] for m in range(n_init)],ddof=1)
    euw = numpy.std([uvwlist[m]['xz'] for m in range(n_init)],ddof=1)
    evw= numpy.std([uvwlist[m]['yz'] for m in range(n_init)],ddof=1)
    ea = numpy.std([uvwlist[m]['a'] for m in range(n_init)],ddof=1)
    eb = numpy.std([uvwlist[m]['b'] for m in range(n_init)],ddof=1)
    ec = numpy.std([uvwlist[m]['c'] for m in range(n_init)],ddof=1)
    x = numpy.mean([xyzlist[m]['x'] for m in range(n_init)])
    y = numpy.mean([xyzlist[m]['y'] for m in range(n_init)])
    z = numpy.mean([xyzlist[m]['z'] for m in range(n_init)])
    xy = numpy.mean([xyzlist[m]['xy'] for m in range(n_init)])
    xz = numpy.mean([xyzlist[m]['xz'] for m in range(n_init)])
    yz= numpy.mean([xyzlist[m]['yz'] for m in range(n_init)])
    d = numpy.mean([xyzlist[m]['a'] for m in range(n_init)])
    e = numpy.mean([xyzlist[m]['b'] for m in range(n_init)])
    f = numpy.mean([xyzlist[m]['c'] for m in range(n_init)])
    ex = numpy.std([xyzlist[m]['x'] for m in range(n_init)],ddof=1)
    ey = numpy.std([xyzlist[m]['y'] for m in range(n_init)],ddof=1)
    ez = numpy.std([xyzlist[m]['z'] for m in range(n_init)],ddof=1)
    exy = numpy.std([xyzlist[m]['xy'] for m in range(n_init)],ddof=1)
    exz = numpy.std([xyzlist[m]['xz'] for m in range(n_init)],ddof=1)
    eyz= numpy.std([xyzlist[m]['yz'] for m in range(n_init)],ddof=1)
    ed = numpy.std([xyzlist[m]['a'] for m in range(n_init)],ddof=1)
    ee = numpy.std([xyzlist[m]['b'] for m in range(n_init)],ddof=1)
    ef = numpy.std([xyzlist[m]['c'] for m in range(n_init)],ddof=1)

        #Output all the particulars of the moving group at this timestep T.
    outfile.write("{0:8.1f}     {1:12.3f} {2:12.3f}  {3:12.3f} {4:12.3f}  {5:12.3f} {6:12.3f}   {7:12.4f} {8:12.4f}  {9:12.4f} {10:12.4f}  {11:12.4f} {12:12.4f}   {13:12.4f} {14:12.4f}  {15:12.4f} {16:12.4f}  {17:12.4f} {18:12.4f}  {19:12.3f} {20:12.3f}  {21:12.3f} {22:12.3f}  {23:12.3f} {24:12.3f}   {25:12.4f} {26:12.4f}  {27:12.4f} {28:12.4f}  {29:12.4f} {30:12.4f}   {31:12.4f} {32:12.4f}  {33:12.4f} {34:12.4f}  {35:12.4f} {36:12.4f}\n".format(0,u,eu,v,ev,w,ew,a,ea,b,eb,c,ec,uv,euv,uw,euw,vw,evw,x,ex,y,ey,z,ez,d,ed,e,ee,f,ef,xy,exy,xz,exz,yz,eyz))

    outfile.close()


if __name__ ==  "__main__":
    traceback()
