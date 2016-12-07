from astropy.io import ascii
from matplotlib import pyplot
from matplotlib.patches import Polygon
import numpy as np
import lacewing
import kinematics
import sys

def traceback(argv=None):
    if argv is None:
        argv = sys.argv
    
    name,coord,era,edec,pmra,epmra,pmdec,epmdec,rv,erv,plx,eplx,note = lacewing.csv_loader(argv[1])
    withmgp = argv[2].replace(' ','_')
    method = argv[3]
    mgpage = argv[4]
    mgpage2 = argv[5]

    timespan = float(argv[6])
    timestep = -0.1
    n_int = np.int(argv[7])

    # load the association for comparison
    mgp = ascii.read('Moving_Group_'+withmgp+'_'+method+'.dat')

    #print timespan/timestep

    time = mgp['Time'][0:int(timespan/timestep)]
    mgpx = mgp['X'][0:int(timespan/timestep)]
    mgpy = mgp['Y'][0:int(timespan/timestep)]
    mgpz = mgp['Z'][0:int(timespan/timestep)]
    mgpa = mgp['A'][0:int(timespan/timestep)]
    mgpea = mgp['eA'][0:int(timespan/timestep)]
    mgpb = mgp['B'][0:int(timespan/timestep)]
    mgpeb = mgp['eB'][0:int(timespan/timestep)]
    mgpc = mgp['C'][0:int(timespan/timestep)]
    mgpec = mgp['eC'][0:int(timespan/timestep)]
        
    # AR 2013.1122 The equivalent radius should be a RADIUS (4/3 pi a b c=4/3 pi r^3)
    mgprad = (mgpa * mgpb * mgpc)**(1./3.)
    mgpradmin = ((mgpa - mgpea) * (mgpb - mgpeb) * (mgpc - mgpec))**(1./3.)
    mgpradmax = ((mgpa + mgpea) * (mgpb + mgpeb) * (mgpc + mgpec))**(1./3.)
    # AR 2014.0321 This is to prevent ill-behaved associations like Argus from screwing up everything.
    mgpradmin[np.where(np.logical_not(np.isfinite(mgpradmin)))] = 0.0

    good_stars = [x for x in xrange(len(coord)) if ((pmra[x] is not None) & (pmdec[x] is not None) & (plx[x] is not None) & (rv[x] is not None))]
    n_stars = len(good_stars)

    for i in good_stars:
        #print name[i],coord[i].ra.degree,era[i],coord[i].dec.degree,edec[i],plx[i],eplx[i],pmra[i],epmra[i],pmdec[i],epmdec[i],rv[i],erv[i]
        print '({0:2d}) {1:16} {2:08.4f} {3:+07.4f} {4:6.2f} {5:+.4f} {6:+.4f} {7:+6.2f}'.format(i,name[i],coord[i].ra.degree,coord[i].dec.degree,1/plx[i],pmra[i],pmdec[i],rv[i])
        sigmin = np.zeros((4,len(time)))+99999.99
        sigmax = np.zeros((4,len(time)))

        if method == 'ballistic':
            px,py,pz = kinematics.ballistic_uniform(coord[i].ra.degree,0,coord[i].dec.degree,0,1/plx[i],0,pmra[i],0,pmdec[i],0,rv[i],0,timespan,timestep,1)
        elif method == 'epicyclic':
            px,py,pz = kinematics.epicyclic_uniform(coord[i].ra.degree,0,coord[i].dec.degree,0,1/plx[i],0,pmra[i],0,pmdec[i],0,rv[i],0,timespan,timestep,1)
        elif method == 'potential':
            px,py,pz = kinematics.potential_uniform(coord[i].ra.degree,0,coord[i].dec.degree,0,1/plx[i],0,pmra[i],0,pmdec[i],0,rv[i],0,timespan,timestep,1)

        distance = np.sqrt((px-mgpx)**2 + (py-mgpy)**2 + (pz-mgpz)**2)
        sigmin[0] = distance
        sigmax[0] = distance

        # now run monte carlos in three distributions
        for k in [1,2,3]:
            if method == 'ballistic':
                px,py,pz = kinematics.ballistic_uniform(coord[i].ra.degree,era[i]*k,coord[i].dec.degree,edec[i]*k,1/plx[i],k*eplx[i]/(plx[i]**2),pmra[i],epmra[i]*k,pmdec[i],epmdec[i]*k,rv[i],erv[i]*k,timespan,timestep,n_int)
            elif method == 'epicyclic':
                px,py,pz = kinematics.epicyclic_uniform(coord[i].ra.degree,era[i]*k,coord[i].dec.degree,edec[i]*k,1/plx[i],k*eplx[i]/(plx[i]**2),pmra[i],epmra[i]*k,pmdec[i],epmdec[i]*k,rv[i],erv[i]*k,timespan,timestep,n_int)
            elif method == 'potential':
                px,py,pz = kinematics.potential_uniform(coord[i].ra.degree,era[i]*k,coord[i].dec.degree,edec[i]*k,1/plx[i],k*eplx[i]/(plx[i]**2),pmra[i],epmra[i]*k,pmdec[i],epmdec[i]*k,rv[i],erv[i]*k,timespan,timestep,n_int)

            # We must rotate these so we are slicing across time, not different stars
            px = np.rot90(px,3)
            py = np.rot90(py,3)
            pz = np.rot90(pz,3)
            # loop through time
            for j in range(len(time)):
                distance = np.sqrt((px[j]-mgpx[j])**2 + (py[j]-mgpy[j])**2 + (pz[j]-mgpz[j])**2)
                sigmin[k,j] = np.amin(distance)
                sigmax[k,j] = np.amax(distance)

        #print sigmin[1]
        x = np.concatenate((time,time[::-1],[time[0]]))
        y1 = np.concatenate((sigmin[1],sigmax[1][::-1],[sigmin[1][0]]))
        y2 = np.concatenate((sigmin[2],sigmax[2][::-1],[sigmin[2][0]]))
        y3 = np.concatenate((sigmin[3],sigmax[3][::-1],[sigmin[3][0]]))
        
        mg = np.concatenate((np.zeros_like(mgpradmin),mgpradmax[::-1],[0]))
        
        fig = pyplot.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        
        poly1 = Polygon(zip(x,y1),facecolor="#000000",edgecolor='none',alpha=0.25)
        poly2 = Polygon(zip(x,y2),facecolor="#1F1F1F",edgecolor='none',alpha=0.25)
        poly3 = Polygon(zip(x,y3),facecolor="#3F3F3F",edgecolor='none',alpha=0.25)
        polym = Polygon(zip(x,mg),facecolor="#FF0000",edgecolor='none',alpha=0.3)
        polya = Polygon(zip([mgpage,mgpage,mgpage2,mgpage2,mgpage],[0,500,500,0,0]),facecolor="#0000FF",edgecolor='none',alpha=0.3)

        ax.add_patch(poly3)
        ax.add_patch(poly2)
        ax.add_patch(poly1)
        ax.add_patch(polym)
        ax.add_patch(polya)
        ax.plot(time,sigmin[0],color="#000000")
        ax.plot(time,mgprad,color="#FF0000")
    
        ax.set_xlim(0,timespan)
        ax.set_ylim(0,200)
        ax.set_xlabel('Time (Myr)')
        ax.set_ylabel('Distance between Star and Moving Group (pc)')
        #ax.set_title('Traceback for {0:} and {3:}, {4:}'.format(name,ra,dec,withmgp,method))
        #ax.xaxis.set_ticks(np.arange(0,timespan,20))
        ax.yaxis.set_ticks(np.arange(0,200,50))
        ax.grid(b=True,which='minor', color="#EFEFEF")
        ax.text(timespan/20.,170,'Real Values',color="#000000")
        ax.text(timespan/20.,160,'1 $\sigma$',color="#9F9F9F")
        ax.text(timespan/20.,150,'2 $\sigma$',color="#BFBFBF")
        ax.text(timespan/20.,140,'3 $\sigma$',color="#DFDFDF")
        ax.text(timespan/20.,130,'Moving Group Volume-Radius (pc)',color="#FF0000")
        ax.text(timespan/20.,120,'Group Age Spread',color="#0000FF")
        #ax.vlines([mgpage],0,500)
    
        pyplot.savefig("{3:}/traceback_{0:}_{1:07.3f}_{2:+07.3f}_{3:}_{4:}.png".format(name[i].replace(' ','_'),coord[i].ra.degree,coord[i].dec.degree,withmgp,method),dpi=100)
        pyplot.clf()
        pyplot.close()
    
if __name__ ==  "__main__":
   if len(sys.argv) == 0:
      print "tracewing.py <inputfile> <group> <method> <minage> <maxage> <maxplotage> <iterations>"
   else:
      traceback()
