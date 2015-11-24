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
    
    readtable = ascii.get_reader(Reader=ascii.Basic)
    readtable.header.splitter.delimiter = ','
    readtable.data.splitter.delimiter = ','
    readtable.header.start_line = 0
    readtable.data.start_line = 1

    inputdata = readtable.read(argv[1])
    withmgp = argv[2]
    method = argv[3]
    mgpage = argv[5]

    timespan = float(argv[4])
    timestep = -0.1
    n_int = int(argv[6])

    for i in np.arange(1,len(inputdata),1):
        # parse all the input file
        star = inputdata[i]
        #print star
        name = star['Name']
        try:
            # all inputs are assumed to be in milliarcseconds, except for RA, DEC, and RV
            ra =   float(star['RA'])
            era =  float(star['eRA'])/3600000
            dec =  float(star['DEC'])
            edec = float(star['eDEC'])/3600000
            plx = float(star['pi'])/1000.
            eplx = float(star['epi'])/1000.
            pmra = float(star['pmRA'])/1000.
            epmra = float(star['epmRA'])/1000.
            pmdec = float(star['pmDEC'])/1000.
            epmdec = float(star['epmDEC'])/1000.
            rv =   float(star['rv'])
            erv =  float(star['erv'])
        except ValueError:
            print name
            continue

        print name,ra,era,dec,edec,plx,eplx,pmra,epmra,pmdec,epmdec,rv,erv

        # load the association for comparison
        mgp = ascii.read('Moving_Group_'+withmgp+'_'+method+'.dat')

        #print association
        print timespan/timestep

        time = mgp['Time'][0:timespan/timestep]
        mgpx = mgp['X'][0:timespan/timestep]
        mgpy = mgp['Y'][0:timespan/timestep]
        mgpz = mgp['Z'][0:timespan/timestep]
        mgpa = mgp['A'][0:timespan/timestep]
        mgpea = mgp['eA'][0:timespan/timestep]
        mgpb = mgp['B'][0:timespan/timestep]
        mgpeb = mgp['eB'][0:timespan/timestep]
        mgpc = mgp['C'][0:timespan/timestep]
        mgpec = mgp['eC'][0:timespan/timestep]
        
        # AR 2013.1122 The equivalent radius should be a RADIUS (4/3 pi a b c=4/3 pi r^3)
        mgprad = (mgpa * mgpb * mgpc)**(1./3.)
        mgpradmin = ((mgpa - mgpea) * (mgpb - mgpeb) * (mgpc - mgpec))**(1./3.)
        mgpradmax = ((mgpa + mgpea) * (mgpb + mgpeb) * (mgpc + mgpec))**(1./3.)
        # AR 2014.0321 This is to prevent ill-behaved associations like Argus from screwing up everything.
        mgpradmin[np.where(np.logical_not(np.isfinite(mgpradmin)))] = 0.0

        sigmin = np.zeros((4,len(time)))+99999.99
        sigmax = np.zeros((4,len(time)))

        if method == 'ballistic':
            px,py,pz = kinematics.ballistic_uniform(ra,0,dec,0,1/plx,0,pmra,0,pmdec,0,rv,0,timespan,timestep,1)
        elif method == 'epicyclic':
            px,py,pz = kinematics.epicyclic_uniform(ra,0,dec,0,1/plx,0,pmra,0,pmdec,0,rv,0,timespan,timestep,1)
        elif method == 'potential':
            px,py,pz = kinematics.potential_uniform(ra,0,dec,0,1/plx,0,pmra,0,pmdec,0,rv,0,timespan,timestep,1)

        distance = np.sqrt((px-mgpx)**2 + (py-mgpy)**2 + (pz-mgpz)**2)
        sigmin[0] = distance
        sigmax[0] = distance

        # now run monte carlos in three distributions
        for k in [1,2,3]:
            if method == 'ballistic':
                px,py,pz = kinematics.ballistic_uniform(ra,era*k,dec,edec*k,1/plx,k*eplx/(plx**2),pmra,epmra*k,pmdec,epmdec*k,rv,erv*k,timespan,timestep,n_int)
            elif method == 'epicyclic':
                px,py,pz = kinematics.epicyclic_uniform(ra,era*k,dec,edec*k,1/plx,k*eplx/(plx**2),pmra,epmra*k,pmdec,epmdec*k,rv,erv*k,timespan,timestep,n_int)
            elif method == 'potential':
                px,py,pz = kinematics.potential_uniform(ra,era*k,dec,edec*k,1/plx,k*eplx/(plx**2),pmra,epmra*k,pmdec,epmdec*k,rv,erv*k,timespan,timestep,n_int)

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

        mg = np.concatenate((mgpradmin,mgpradmax[::-1],[mgpradmin[0]]))
        
        fig = pyplot.figure(figsize=(7,5))
        ax = fig.add_subplot(111)
        
        poly1 = Polygon(zip(x,y1),facecolor="#000000",edgecolor='none',alpha=0.25)
        poly2 = Polygon(zip(x,y2),facecolor="#1F1F1F",edgecolor='none',alpha=0.25)
        poly3 = Polygon(zip(x,y3),facecolor="#3F3F3F",edgecolor='none',alpha=0.25)
        polym = Polygon(zip(x,mg),facecolor="#FF0000",edgecolor='none',alpha=0.3)

        ax.add_patch(polym)
        ax.add_patch(poly3)
        ax.add_patch(poly2)
        ax.add_patch(poly1)
        ax.plot(time,sigmin[0],color="#000000")
        ax.plot(time,mgprad,color="#FF0000")
        
        ax.set_xlim(0,timespan)
        ax.set_ylim(0,160)
        ax.set_xlabel('Myr ago')
        ax.set_ylabel('Distance between Star and Moving Group')
        ax.set_title('Traceback for {0:} and {3:}, {4:}'.format(name,ra,dec,withmgp,method))
#        ax.xaxis.set_ticks(np.arange(0,timespan,20))
        ax.yaxis.set_ticks(np.arange(0,100,20))
        ax.grid(b=True,which='minor', color="#EFEFEF")
        ax.text(-10,140,'Real Values',color="#000000")
        ax.text(-10,135,'1 $\sigma$',color="#9F9F9F")
        ax.text(-10,130,'2 $\sigma$',color="#BFBFBF")
        ax.text(-10,125,'3 $\sigma$',color="#DFDFDF")
        ax.text(-10,115,'Moving Group Volume-Radius',color="#FF0000")
        ax.vlines([mgpage],0,500)
    
        pyplot.savefig("traceback_{0:}_{1:07.3f}_{2:+07.3f}_{3:}_{4:}.png".format(name,ra,dec,withmgp,method),dpi=100)
        pyplot.clf()
        pyplot.close()
    
if __name__ ==  "__main__":
   traceback()
