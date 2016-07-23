import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import pyplot
#from matplotlib import cm
from matplotlib.patches import Ellipse
import kinematics
import ellipse
import lacewing
import sys
from astropy.io import ascii

def traceback(argv=None):
    if argv is None:
        argv = sys.argv
        
    name,coords,era,edec,pmra,epmra,pmdec,epmdec,rv,erv,pi,epi,note = lacewing.csv_loader(argv[1])
    mgpname = argv[2]
    n_init = int(argv[3])
    #print name,era
    era = np.asarray(era)
    edec = np.asarray(era)
    pmra = np.asarray(pmra)
    epmra = np.asarray(epmra)
    pmdec = np.asarray(pmdec)
    epmdec = np.asarray(epmdec)
    pi = np.asarray(pi)
    epi = np.asarray(epi)
    rv = np.asarray(rv)
    erv = np.asarray(erv)


    # How many stars are we fitting?
    goodlist = [x for x in range(len(pi)) if (pi[x] is not None) and (rv[x] is not None) and (pmra[x] is not None) and (pmdec[x] is not None)]

    ra = np.array(())
    dec = np.array(())
    for b in goodlist:
        ra = np.append(ra,coords[b].ra.degree)
        dec = np.append(dec,coords[b].dec.degree)

    era = era[goodlist]
    edec = edec[goodlist]
    pmra = pmra[goodlist]
    epmra = epmra[goodlist]
    pmdec = pmdec[goodlist]
    epmdec = epmdec[goodlist]
    pi = pi[goodlist]
    epi = epi[goodlist]
    rv = rv[goodlist]
    erv = erv[goodlist]

    #print goodlist
    print len(goodlist)
    n_stars = len(goodlist)

    uvlist = []
    uwlist = []
    vwlist = []
    xylist = []
    xzlist = []
    yzlist = []
    outfile = open("Moving_Group_{0:}_2d.dat".format(mgpname).replace(' ','_'),"wb")
    outfile.write("U,V,W,A,B,C,UV,UW,VW,X,Y,Z,D,E,G,XY,XZ,YZ\n")
    for j in xrange(n_init):
        tra = ra + (np.random.randn(n_stars)*era)*np.cos(dec*np.pi/180.)
        tdec = dec + np.random.randn(n_stars)*edec
        tpi = pi + np.random.randn(n_stars)*epi
        tpmra = pmra + np.random.randn(n_stars)*epmra 
        tpmdec = pmdec + np.random.randn(n_stars)*epmdec 
        trv = rv + np.random.randn(n_stars)*erv
        #print j
        tu,tv,tw,tx,ty,tz = kinematics.gal_uvwxyz(ra=tra,dec=tdec,plx=tpi,pmra=tpmra,pmdec=tpmdec,vrad=trv)
        tu = np.array(tu,dtype=np.float64)
        tv = np.array(tv,dtype=np.float64)
        tw = np.array(tw,dtype=np.float64)
        tx = np.array(tx,dtype=np.float64)
        ty = np.array(ty,dtype=np.float64)
        tz = np.array(tz,dtype=np.float64)
        objxy = ellipse.fitellipse2d(tx,ty)
        objxz = ellipse.fitellipse2d(tx,tz)
        objyz = ellipse.fitellipse2d(ty,tz)
        xylist.append(objxy)
        xzlist.append(objxz)
        yzlist.append(objyz)
        objuv = ellipse.fitellipse2d(tu,tv)
        objuw = ellipse.fitellipse2d(tu,tw)
        objvw = ellipse.fitellipse2d(tv,tw)
        uvlist.append(objuv)
        uwlist.append(objuw)
        vwlist.append(objvw)

    u = np.mean([uvlist[m]['x'] for m in range(n_init)])
    v = np.mean([vwlist[m]['x'] for m in range(n_init)])
    w = np.mean([uwlist[m]['y'] for m in range(n_init)])
    uv = np.mean([uvlist[m]['xy'] for m in range(n_init)])
    uw = np.mean([vwlist[m]['xy'] for m in range(n_init)])
    vw= np.mean([uwlist[m]['xy'] for m in range(n_init)])
    a = np.mean([uvlist[m]['a'] for m in range(n_init)])
    b = np.mean([vwlist[m]['a'] for m in range(n_init)])
    c = np.mean([uwlist[m]['b'] for m in range(n_init)])
    eu = np.std([uvlist[m]['x'] for m in range(n_init)],ddof=1)
    ev = np.std([vwlist[m]['x'] for m in range(n_init)],ddof=1)
    ew = np.std([uwlist[m]['y'] for m in range(n_init)],ddof=1)
    euv = np.std([uvlist[m]['xy'] for m in range(n_init)],ddof=1)
    euw = np.std([vwlist[m]['xy'] for m in range(n_init)],ddof=1)
    evw= np.std([uwlist[m]['xy'] for m in range(n_init)],ddof=1)
    ea = np.std([uvlist[m]['a'] for m in range(n_init)],ddof=1)
    eb = np.std([vwlist[m]['a'] for m in range(n_init)],ddof=1)
    ec = np.std([uwlist[m]['b'] for m in range(n_init)],ddof=1)
    x = np.mean([xylist[m]['x'] for m in range(n_init)])
    y = np.mean([yzlist[m]['x'] for m in range(n_init)])
    z = np.mean([xzlist[m]['y'] for m in range(n_init)])
    xy = np.mean([xylist[m]['xy'] for m in range(n_init)])
    xz = np.mean([xzlist[m]['xy'] for m in range(n_init)])
    yz= np.mean([yzlist[m]['xy'] for m in range(n_init)])
    d = np.mean([xylist[m]['a'] for m in range(n_init)])
    e = np.mean([yzlist[m]['a'] for m in range(n_init)])
    f = np.mean([xzlist[m]['b'] for m in range(n_init)])
    ex = np.std([xylist[m]['x'] for m in range(n_init)],ddof=1)
    ey = np.std([yzlist[m]['x'] for m in range(n_init)],ddof=1)
    ez = np.std([xzlist[m]['y'] for m in range(n_init)],ddof=1)
    exy = np.std([xylist[m]['xy'] for m in range(n_init)],ddof=1)
    exz = np.std([xzlist[m]['xy'] for m in range(n_init)],ddof=1)
    eyz = np.std([yzlist[m]['xy'] for m in range(n_init)],ddof=1)
    ed = np.std([xylist[m]['a'] for m in range(n_init)],ddof=1)
    ee = np.std([yzlist[m]['a'] for m in range(n_init)],ddof=1)
    ef = np.std([xzlist[m]['b'] for m in range(n_init)],ddof=1)

        #Output all the particulars of the moving group at this timestep T.
    outfile.write("{0:12.3f},{1:12.3f},{2:12.3f},{3:12.3f},{4:12.3f},{5:12.3f},{6:12.4f},{7:12.4f},{8:12.4f},{9:12.4f},{10:12.4f},{11:12.4f},{12:12.4f},{13:12.4f},{14:12.4f},{15:12.4f},{16:12.4f},{17:12.4f}\n".format(u,v,w,a,b,c,uv,uw,vw,x,y,z,d,e,f,xy,xz,yz))
    outfile.write("{0:12.3f},{1:12.3f},{2:12.3f},{3:12.3f},{4:12.3f},{5:12.3f},{6:12.4f},{7:12.4f},{8:12.4f},{9:12.4f},{10:12.4f},{11:12.4f},{12:12.4f},{13:12.4f},{14:12.4f},{15:12.4f},{16:12.4f},{17:12.4f}\n".format(eu,ev,ew,ea,eb,ec,euv,euw,evw,ex,ey,ez,ed,ee,ef,exy,exz,eyz))

    outfile.close()


if __name__ ==  "__main__":
    traceback()
