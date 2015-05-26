import numpy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot
#from matplotlib import cm
from matplotlib.patches import Ellipse
from matplotlib.patches import Polygon
from matplotlib import patches
from matplotlib import _png
import kinematics
import ellipse
import astrometry
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
        _png.write_png(renderer._renderer.buffer_rgba(),renderer.width,renderer.height,outfile,fig.dpi)

def traceback(argv=None):
    if argv is None:
        argv = sys.argv
        
    mgpname = argv[2]
    method = argv[3]
    mgpage = numpy.float(argv[4])

    timespan = numpy.float(argv[5])
    timestep = -0.1
    n_int = 5000
   
    readtable = ascii.get_reader(Reader=ascii.Basic)
    readtable.header.splitter.delimiter = ','
    readtable.data.splitter.delimiter = ','
    readtable.header.start_line = 1
    readtable.data.start_line = 3

    young = readtable.read(argv[1])
   
    # How many stars are we fitting? Find all that are bona fide members (quality is "Good")
    n_stars=0

    n_stars = len(numpy.where(numpy.bitwise_and(astrometry.isnumber(young['rv']),astrometry.isnumber(young['pi'])) & (young['GROUP'] == mgpname) & (young['GROUP quality'] == "Good"))[0])

    print n_stars

    # it saves time and memory to make arrays in advance, even in python
    # Set dtype=numpy.float32 to save memory
    mgp_x = numpy.zeros((n_stars,n_int,6000),dtype=numpy.float32)
    mgp_y = numpy.zeros((n_stars,n_int,6000),dtype=numpy.float32)
    mgp_z = numpy.zeros((n_stars,n_int,6000),dtype=numpy.float32)
    # These are for the 3D movie version
    mgp_size = []
    mgp_color = []
    mgp_n = []

    print mgp_x.shape

    #print young

    n = 0
    for i in xrange(len(young)):
        if isinstance(young[i]['RA'], numpy.ma.core.MaskedConstant):
            continue
        else:
            #print young[i]
            if young[i]['GROUP'] == mgpname:
                # key assumption: All objects that are Good have all the relevant data filled in.
                if young[i]['GROUP quality'] == "Good":
                    #print n
                    try:
                        if isinstance(young[i]['SpType'], numpy.ma.core.MaskedConstant):
                            if isinstance(young[i]['V-K'], numpy.ma.core.MaskedConstant):
                                size=20
                                color='#000000'
                            else:
                                if young[i]['V-K'] <= -0.3:  #O
                                    size = 60
                                    color= "#3F00FF"
                                elif young[i]['V-K'] <= 0 and young[i]['V-K'] > -0.3: #B
                                    size = 50
                                    color= "#3F3FFF"
                                elif young[i]['V-K'] <= 0.5 and young[i]['V-K'] > 0: #A
                                    size = 40
                                    color= "#7F7FFF"
                                elif young[i]['V-K'] <= 1.0 and young[i]['V-K'] > 0.5: #F
                                    size = 30
                                    color= "#CFCFFF"
                                elif young[i]['V-K'] <= 2.0 and young[i]['V-K'] > 1.0: #G
                                    size = 20
                                    color= "#CFCFCF"
                                elif young[i]['V-K'] <= 3.8 and young[i]['V-K'] > 2.0: #K
                                    size = 15
                                    color= "#FFDF1F"
                                elif young[i]['V-K'] <= 9.5 and young[i]['V-K'] > 3.8: #M
                                    size = 10
                                    color= "#FF3F00"
                                elif young[i]['V-K'] > 9.5: #L
                                    size = 5
                                    color = "#FF0000"
                                    
                        else:
                            if len(young[i]['SpType']) == 0:
                                size = 20
                                color= "#000000"
                            elif young[i]['SpType'][0] == 'O':
                                size = 60
                                color= "#3F00FF"
                            elif young[i]['SpType'][0] == 'B':
                                size = 50
                                color= "#3F3FFF"
                            elif young[i]['SpType'][0] == 'A':
                                size = 40
                                color= "#7F7FFF"
                            elif young[i]['SpType'][0] == 'F':
                                size = 30
                                color= "#CFCFFF"
                            elif young[i]['SpType'][0] == 'G':
                                size = 20
                                color= "#CFCFCF"
                            elif young[i]['SpType'][0] == 'K':
                                size = 15
                                color= "#FFDF1F"
                            elif young[i]['SpType'][0] == 'M':
                                size = 10
                                color= "#FF3F00"
                            elif young[i]['SpType'][0] == 'L':
                                size = 5
                                color = "#FF0000"
                        ra =   float(young[i]['RA'])
                        era =  float(young[i]['eRA'])/3600000
                        dec =  float(young[i]['DEC'])
                        edec = float(young[i]['eDEC'])/3600000
                        dist = 1000/float(young[i]['pi'])
                        edist = float(young[i]['epi'])/float(young[i]['pi'])*dist
                        pmra = float(young[i]['pmRA'])/1000.
                        epmra = float(young[i]['epmRA'])/1000.
                        pmdec = float(young[i]['pmDEC'])/1000.
                        epmdec = float(young[i]['epmDEC'])/1000.
                        rv =   float(young[i]['rv'])
                        erv =  float(young[i]['erv'])
                        print n, young[i]['Name'],ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv
                    except ValueError:
                        continue
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
                        px,py,pz = kinematics.ballistic(ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv,-600,timestep,n_int)
                    elif method == 'epicyclic':
                        px,py,pz = kinematics.epicyclic(ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv,-600,timestep,n_int)
                    elif method == 'potential':
                        px,py,pz = kinematics.potential(ra,era,dec,edec,dist,edist,pmra,epmra,pmdec,epmdec,rv,erv,-600,timestep,n_int)
                    # store these iterations
                    mgp_x[n] = px
                    mgp_y[n] = py
                    mgp_z[n] = pz
                    mgp_n.append(young[i]['Name'])
                    
                    mgp_color.extend([color]*n_int)
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
    mgp_x = numpy.asarray(mgp_x,dtype=numpy.float32)
    mgp_y = numpy.asarray(mgp_y,dtype=numpy.float32)
    mgp_z = numpy.asarray(mgp_z,dtype=numpy.float32)
    times = numpy.arange(0,-600,timestep)
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
        x = numpy.mean([objlist[m]['x'] for m in range(n_int)])
        y = numpy.mean([objlist[m]['y'] for m in range(n_int)])
        z = numpy.mean([objlist[m]['z'] for m in range(n_int)])
        xy = numpy.mean([objlist[m]['xy'] for m in range(n_int)])
        xz = numpy.mean([objlist[m]['xz'] for m in range(n_int)])
        yz= numpy.mean([objlist[m]['yz'] for m in range(n_int)])
        a = numpy.mean([objlist[m]['a'] for m in range(n_int)])
        b = numpy.mean([objlist[m]['b'] for m in range(n_int)])
        c = numpy.mean([objlist[m]['c'] for m in range(n_int)])
        ex = numpy.std([objlist[m]['x'] for m in range(n_int)],ddof=1)
        ey = numpy.std([objlist[m]['y'] for m in range(n_int)],ddof=1)
        ez = numpy.std([objlist[m]['z'] for m in range(n_int)],ddof=1)
        exy = numpy.std([objlist[m]['xy'] for m in range(n_int)],ddof=1)
        exz = numpy.std([objlist[m]['xz'] for m in range(n_int)],ddof=1)
        eyz= numpy.std([objlist[m]['yz'] for m in range(n_int)],ddof=1)
        ea = numpy.std([objlist[m]['a'] for m in range(n_int)],ddof=1)
        eb = numpy.std([objlist[m]['b'] for m in range(n_int)],ddof=1)
        ec = numpy.std([objlist[m]['c'] for m in range(n_int)],ddof=1)

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
    mgp_x = numpy.reshape(mgp_x,(n_stars*n_int,6000))
    mgp_y = numpy.reshape(mgp_y,(n_stars*n_int,6000))
    mgp_z = numpy.reshape(mgp_z,(n_stars*n_int,6000))
    
    # rotate so that each strip contains n_stars*n_int elements for a given time.
    # (it's easier to plot)
    mgp_x = numpy.rot90(mgp_x,3)
    mgp_y = numpy.rot90(mgp_y,3)
    mgp_z = numpy.rot90(mgp_z,3)
    times = numpy.rot90([times],3)
    mgp_color = numpy.reshape(mgp_color,-1)

    #######################################################
    ### Draw a traceback plot of all N*I stars relative ###
    ### to their computed center, to show the behavior  ###
    ### of the moving group itself. This is the         ###
    ###             "waterfall diagram"                 ###
    #######################################################

    fig2 = pyplot.figure(figsize=(9.6,5.4),dpi=200)
    ax2 = fig2.add_subplot(111)
    ax2.set_ylim((0,100))
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
    a = numpy.asarray(a)
    b = numpy.asarray(b)
    c = numpy.asarray(c)
    x = numpy.asarray(x)
    y = numpy.asarray(y)
    z = numpy.asarray(z)
    for p in numpy.arange(0,n_stars*n_int,3):
        mgpdist = numpy.sqrt((mgp_x[:,p] - x)**2 + (mgp_y[:,p] - y)**2 + (mgp_z[:,p] - z)**2)
        line.set_data(times,mgpdist)
        line.set_linewidth(0.2)
        line.set_color((0,0,0,0.05))
        ax2.draw_artist(line)
    line.set_data(times,(a*b*c)**(1/3.))
    line.set_color((1,0,0,1))
    line.set_linewidth(1)
    ax2.draw_artist(line)

    x = numpy.concatenate((times,times[::-1],[times[0]]))
        
    mgprad = (a*b*c)**(1/3.)
    mg = numpy.concatenate((mgprad,numpy.zeros_like(mgprad),[mgprad[0]]))
    poly.set_xy(zip(*(x,mg)))
    poly.set_facecolor((1,0,0,0.2))
    ax2.draw_artist(poly)

    line.set_data([mgpage,mgpage],[0,500])
    line.set_color((1,0,0,1))
    line.set_linewidth(1)
    ax2.draw_artist(line)
    save(fig2,'Trace_{0:}_{1:}.png'.format(mgpname,method))
    fig2.clf()
    pyplot.close()

#    #####################################################
#    ### The 3D plot. We're going to plot all of the   ###
#    ###  iterations of the moving group so we can see ###
#    ###  the full effect. Also, it looks like an      ###
#    ###  explosion, and that's always cool.           ###
#    #####################################################
#    
#    print mgp_x[0]   
#    for k in xrange(len(times)):
#        print times[k]
#        fig = pyplot.figure(figsize=(9.6,5.5))
#        ax = fig.add_subplot(111,projection='3d', aspect='equal')
#
#        mgp_tx,mgp_ty,mgp_tz = ellipse.triaxial(mgpmaster[k])
#
#        #Add the Sun
#        #ax.scatter(0-mgpmaster[k]['x'],0-mgpmaster[k]['y'],0-mgpmaster[k]['z'],s=50,color="#FFCF00")
#
#        ax.scatter(mgp_x[k]-mgpmaster[k]['x'],mgp_y[k]-mgpmaster[k]['y'],mgp_z[k]-mgpmaster[k]['z'],c=mgp_color,s=1,linewidth=0.02,alpha=0.5)
#        ax.plot_surface(mgp_tx-mgpmaster[k]['x'],mgp_ty-mgpmaster[k]['y'],mgp_tz-mgpmaster[k]['z'],rstride=4,cstride=4,color="#000000",alpha=0.1,linewidth=0.02)
#
#        ax.set_xlim(-200,200)
#        ax.set_ylim(-200,200)
#        ax.set_zlim(-200,200)
#        ax.plot([-200,200],[0,0],[0,0],"k")
#        ax.plot([0,0],[-200,200],[0,0],"k")
#        ax.plot([0,0],[0,0],[-200,200],"k")
#        ax.set_xlabel('$\Delta$X (pc)')
#        ax.set_ylabel('$\Delta$Y (pc)')
#        ax.set_zlabel('$\Delta$Z (pc)')
#        ax.set_title("{0: 4.1f}".format(times[k][0])+" Myr")
#        ax.text(-500,-350,400,argv[2],color="#000000")
#
#        pyplot.savefig('{0:}-{1:}/traceback_{2:04d}_gagne.png'.format(method,mgpname,k), dpi=200)
#        pyplot.clf()
#        pyplot.close()
#      
#        print mgpmaster[k]

if __name__ ==  "__main__":
    traceback()
