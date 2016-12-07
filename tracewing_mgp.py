import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot
#from matplotlib import cm
from matplotlib.patches import Ellipse
from matplotlib.patches import Polygon
from matplotlib import patches
from matplotlib import _png
import kinematics
import ellipse
import astrometry
import lacewing
import sys
from astropy.io import ascii

# Normally, matplotlib just builds up a list of plot elements until you 
#  call savefig(), when it renders them all at once. That uses a LOT 
#  of memory. This function replaces savefig() with something that plots
#  every element as it comes in.
def save(fig,filename):
    #We have to work around 'fig.canvas.print_png', etc calling 'draw'
    renderer = fig.canvas.renderer
    with open(filename,'w') as outfile:
        _png.write_png(renderer._renderer,outfile,fig.dpi)

def traceback(argv=None):
    if argv is None:
        argv = sys.argv
        
    mgpname = argv[2]
    method = argv[3]
    mgpage = np.float(argv[4])
    mgpage2 = np.float(argv[5])

    timespan = np.float(argv[6])
    timestep = -0.1
    n_int = int(1000)
    full_timespan = -800
   
    name,coord,era,edec,pmra,epmra,pmdec,epmdec,rv,erv,plx,eplx,note = lacewing.csv_loader(argv[1])
   
    # How many stars are we fitting? 
    good_stars = [x for x in xrange(len(coord)) if ((pmra[x] is not None) & (pmdec[x] is not None) & (plx[x] is not None) & (rv[x] is not None))]
    n_stars = len(good_stars)
    print 'Number of stars in solution: {0:}'.format(n_stars)

    # it saves time and memory to make arrays in advance, even in python
    # Set dtype=np.float32 to save memory
    mgp_x = np.zeros((n_stars,n_int,np.int(np.ceil(full_timespan/timestep))),dtype=np.float32)
    mgp_y = np.zeros((n_stars,n_int,np.int(np.ceil(full_timespan/timestep))),dtype=np.float32)
    mgp_z = np.zeros((n_stars,n_int,np.int(np.ceil(full_timespan/timestep))),dtype=np.float32)
    # These are for the 3D movie version
    mgp_size = []
    #mgp_color = []
    mgp_n = []

    n = 0
    print "  #  Name             RA        DEC     Dist.   pmRA    pmDEC    RV"
    for i in good_stars:
        print '({0:2d}) {1:16} {2:08.4f} {3:+07.4f} {4:6.2f} {5:+.4f} {6:+.4f} {7:+6.2f}'.format(i,name[i],coord[i].ra.degree,coord[i].dec.degree,1/plx[i],pmra[i],pmdec[i],rv[i])
        ###############################################################
        ### We now have the particulars about one star.  We are now ###
        ###  going to run I monte carlo iterations through the      ###
        ###  specified traceback method. The traceback is going to  ###
        ###  compute UVW points and run them back in time. We could ###
        ###  then fit an ellipse to this star and save only the     ###
        ###  ellipse parameters, but then information would be lost ###
        ###  (or need to be re-created) when we want to determine   ###
        ###  the shape of the moving group itself. Instead, we'll   ###
        ###  save and use every single monte carlo iteration.       ###
        ###############################################################
                    
        if method == 'ballistic':
            px,py,pz = kinematics.ballistic(coord[i].ra.degree,era[i],coord[i].dec.degree,edec[i],1/plx[i],eplx[i]/(plx[i]**2),pmra[i],epmra[i],pmdec[i],epmdec[i],rv[i],erv[i],full_timespan,timestep,n_int)
        elif method == 'epicyclic':
            px,py,pz = kinematics.epicyclic(coord[i].ra.degree,era[i],coord[i].dec.degree,edec[i],1/plx[i],eplx[i]/(plx[i]**2),pmra[i],epmra[i],pmdec[i],epmdec[i],rv[i],erv[i],full_timespan,timestep,n_int)
        elif method == 'potential':
            px,py,pz = kinematics.potential(coord[i].ra.degree,era[i],coord[i].dec.degree,edec[i],1/plx[i],eplx[i]/(plx[i]**2),pmra[i],epmra[i],pmdec[i],epmdec[i],rv[i],erv[i],full_timespan,timestep,n_int)
        # store these iterations
        mgp_x[n] = px
        mgp_y[n] = py
        mgp_z[n] = pz
        mgp_n.append(name[i])
            
        #mgp_color.extend([color]*n_int)
        n = n+1
    # remove spaces from name - helps with programming later on
    mgpname = mgpname.replace(' ', '_')

    ################################################
    ###  At this point in the program, we have   ###
    ### an NxIxT grid of positions constituting  ###
    ### the positions of I iterations of N stars ###
    ### at T times. Now we must fit ellipses to  ###
    ###   those values at every time T.          ###
    ################################################

    # AR 2014.0319: Based on an idea from Adrian Price-Whelan, rather than calculating the full tracebacks of n*1000 stars over X Myr, I'm going to calculate 1000 tracebacks of n stars over X Myr.  
    mgp_x = np.asarray(mgp_x,dtype=np.float32)
    mgp_y = np.asarray(mgp_y,dtype=np.float32)
    mgp_z = np.asarray(mgp_z,dtype=np.float32)
    times = np.arange(0,full_timespan,timestep)
    ## output positions of individual stars as a function of time.
    #for s in range(len(times)):
    #    outfile = open("mgp_{0:}_{1:}_{2:}.csv".format(mgpname,method,times[s]),"wb")
    #    outfile.write("Name,X,Y,Z,A,B,C,XY,XZ,YZ\n")
    #    for t in range(len(mgp_x)):
    #        obj = ellipse.fitellipse(mgp_x[t,:,s],mgp_y[t,:,s],mgp_z[t,:,s])
    #
    #        outfile.write("{0:}, {1: 12.8f}, {2: 12.8f}, {3: 12.8f}, {4: 12.8f}, {5: 12.8f}, {6: 12.8f}, {7: 12.8f}, {8: 12.8f}, {9: 12.8f}\n".format(mgp_n[t],obj['x'],obj['y'],obj['z'],obj['a'],obj['b'],obj['c'],obj['xy'],obj['xz'],obj['yz']))
    #    outfile.close()

    mgpmaster = []
    outfile = open("Moving_Group_{0:}_{1:}.dat".format(mgpname,method),"wb")
    outfile.write("Time X eX Y eY Z eZ A eA B eB C eC XY eXY XZ eXZ YZ eYZ\n")
    for k in xrange(len(times)):
        objlist = []
        for j in xrange(n_int):
            # this is one ellipse per monte carlo iteration
            obj = ellipse.fitellipse(mgp_x[:,j,k],mgp_y[:,j,k],mgp_z[:,j,k])
            objlist.append(obj)
        x = np.mean([objlist[m]['x'] for m in range(n_int)])
        y = np.mean([objlist[m]['y'] for m in range(n_int)])
        z = np.mean([objlist[m]['z'] for m in range(n_int)])
        xy = np.mean([objlist[m]['xy'] for m in range(n_int)])
        xz = np.mean([objlist[m]['xz'] for m in range(n_int)])
        yz= np.mean([objlist[m]['yz'] for m in range(n_int)])
        a = np.mean([objlist[m]['a'] for m in range(n_int)])
        b = np.mean([objlist[m]['b'] for m in range(n_int)])
        c = np.mean([objlist[m]['c'] for m in range(n_int)])
        ex = np.std([objlist[m]['x'] for m in range(n_int)],ddof=1)
        ey = np.std([objlist[m]['y'] for m in range(n_int)],ddof=1)
        ez = np.std([objlist[m]['z'] for m in range(n_int)],ddof=1)
        exy = np.std([objlist[m]['xy'] for m in range(n_int)],ddof=1)
        exz = np.std([objlist[m]['xz'] for m in range(n_int)],ddof=1)
        eyz= np.std([objlist[m]['yz'] for m in range(n_int)],ddof=1)
        ea = np.std([objlist[m]['a'] for m in range(n_int)],ddof=1)
        eb = np.std([objlist[m]['b'] for m in range(n_int)],ddof=1)
        ec = np.std([objlist[m]['c'] for m in range(n_int)],ddof=1)

        # We're re-saving a dictionary of one ellipse per TIMESTEP so that we can make a 3D plot of it later.
        mgpmaster.append({'x':x,'ex':ex,'y':y,'ey':ey,'z':z,'ez':ez,'xy':xy,'exy':exy,'xz':xz,'exz':exz,'yz':yz,'eyz':eyz,'a':a,'ea':ea,'b':b,'eb':eb,'c':c,'ec':ec})
        #Output all the particulars of the moving group at this timestep T.
        outfile.write("{0:8.1f}     {1:12.3f} {2:12.3f}  {3:12.3f} {4:12.3f}  {5:12.3f} {6:12.3f}   {7:12.4f} {8:12.4f}  {9:12.4f} {10:12.4f}  {11:12.4f} {12:12.4f}   {13:12.4f} {14:12.4f}  {15:12.4f} {16:12.4f}  {17:12.4f} {18:12.4f}\n".format(times[k],mgpmaster[k]['x'],mgpmaster[k]['ex'],mgpmaster[k]['y'],mgpmaster[k]['ey'],mgpmaster[k]['z'],mgpmaster[k]['ez'],mgpmaster[k]['a'],mgpmaster[k]['ea'],mgpmaster[k]['b'],mgpmaster[k]['eb'],mgpmaster[k]['c'],mgpmaster[k]['ec'],mgpmaster[k]['xy'],mgpmaster[k]['exy'],mgpmaster[k]['xz'],mgpmaster[k]['exz'],mgpmaster[k]['yz'],mgpmaster[k]['eyz']))

    outfile.close()

    #####################################################
    ### Flatten everything! Now we can plot all our   ###
    ### iterations at once (for the waterfall diagram ###
    ### and the 3D explosion plot)                    ###
    #####################################################

    # flatten by one dimension; we now have N*I stars at every time T.
    mgp_x = np.reshape(mgp_x,(n_stars*n_int,np.ceil(full_timespan/timestep)))
    mgp_y = np.reshape(mgp_y,(n_stars*n_int,np.ceil(full_timespan/timestep)))
    mgp_z = np.reshape(mgp_z,(n_stars*n_int,np.ceil(full_timespan/timestep)))
    
    ## rotate so that each strip contains n_stars*n_int elements for a given time.
    ## (it's easier to plot)
    #mgp_x = np.rot90(mgp_x,3)
    #mgp_y = np.rot90(mgp_y,3)
    #mgp_z = np.rot90(mgp_z,3)
    #times = np.rot90([times],3)
    #mgp_color = np.reshape(mgp_color,-1)

    #######################################################
    ### Draw a traceback plot of all N*I stars relative ###
    ### to their computed center, to show the behavior  ###
    ### of the moving group itself. This is the         ###
    ###             "waterfall diagram"                 ###
    #######################################################

    fig2 = pyplot.figure(figsize=(9.6,5.4),dpi=600)
    ax2 = fig2.add_subplot(111)
    ax2.set_ylim((0,200))
    ax2.set_xlim((0,timespan))
    ax2.set_xlabel('Time (Myr)')
    ax2.set_ylabel('(pc)')

    fig2.canvas.draw()

    line = ax2.plot([0,1],[0,1],color=(1,1,1,0.1),linewidth=1)
    line = line[0]
    poly = ax2.add_patch(patches.Polygon([[0,1],[1,2],[2,0]],closed=True,facecolor=(1,1,1,0.15), edgecolor='none'))

    #ax2.set_title('{0:} {1:} traceback'.format(mgpname,method))
    a = []
    b = []
    c = []
    x = []
    y = []
    z = []
    for q in range(len(times)):
        a.append(mgpmaster[q]["a"])
        b.append(mgpmaster[q]["b"])
        c.append(mgpmaster[q]["c"])
        x.append(mgpmaster[q]["x"])
        y.append(mgpmaster[q]["y"])
        z.append(mgpmaster[q]["z"])
    a = np.asarray(a)
    b = np.asarray(b)
    c = np.asarray(c)
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    print a
    # Control the crowding/density of the plots. Only plot 30,000 curves regardless of how many there really are.
    if n_stars*n_int > 30000:
        p = np.asarray(np.ceil(np.random.rand(30000)*(n_stars*n_int-1)),np.int)
    else:
        p = np.arange(0,(n_stars*n_int-1),1)
    mgpdist = np.sqrt((mgp_x[p,:] - x)**2 + (mgp_y[p,:] - y)**2 + (mgp_z[p,:] - z)**2)
    for cntr in range(len(mgpdist)):
        line.set_data(times,mgpdist[cntr])
        line.set_linewidth(0.2)
        line.set_color((0,0,0,0.05))
        ax2.draw_artist(line)
    line.set_data(times,(a*b*c)**(1/3.))
    line.set_color((1,0,0,0.5))
    line.set_linewidth(1)
    ax2.draw_artist(line)

    x = np.concatenate((times,times[::-1],[times[0]]))
        
    mgprad = (a*b*c)**(1/3.)
    mg = np.concatenate((mgprad,np.zeros_like(mgprad),[mgprad[0]]))
    poly.set_xy(zip(*(x,mg)))
    poly.set_facecolor((1,0,0,0.2))
    ax2.draw_artist(poly)

    poly.set_xy(zip(*([mgpage,mgpage,mgpage2,mgpage2,mgpage],[0,500,500,0,0])))
    poly.set_facecolor((0,0,1,0.2))
    ax2.draw_artist(poly)

    save(fig2,'Trace_{0:}_{1:}.png'.format(mgpname,method))
    fig2.clf()
    pyplot.close()

if __name__ ==  "__main__":
   if len(sys.argv) == 0:
      print "tracewing_mgp.py <inputfile> <group> <method> <minage> <maxage> <maxplotage>"
   else:
      traceback()
