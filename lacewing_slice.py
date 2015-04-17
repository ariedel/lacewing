import numpy
import scipy
from matplotlib import pyplot
from matplotlib import cm
from matplotlib.patches import Ellipse
from matplotlib import font_manager
from scipy.stats import norm
import astropy.io.ascii as ascii
from sys import argv
import astrometry
import kinematics
import ellipse

###################################
###################################
### MAIN ROUTINE
###################################
###################################

#prop = font_manager.FontProperties(fname='/home/riedel/Arial.ttf')

try:
    infilename = argv[1]
except IndexError: 
    print 'syntax: python lacewing_slice.py inputfile XYZ outfile'
try:
   XYZ = bool(argv[2])
except IndexError:
   XYZ = False
try:
    outname = argv[3]
except IndexError:
    outname = 'lacewing.uvwxyz.output'
#try:
#   threed = argv[4]
#except IndexError:
#   threed = False

class Mgp:
    def __init__(self,mgp):
        self.name = mgp["Name"]
        self.U = numpy.float(mgp["U"])
        self.A = numpy.float(mgp["A"])
        self.V = numpy.float(mgp["V"])
        self.B = numpy.float(mgp["B"])
        self.W = numpy.float(mgp["W"])
        self.C = numpy.float(mgp["C"])
        self.UV = numpy.float(mgp["UV"])
        self.UW = numpy.float(mgp["UW"])
        self.VW = numpy.float(mgp["VW"])
        self.X = numpy.float(mgp["X"])
        self.D = numpy.float(mgp["D"])
        self.Y = numpy.float(mgp["Y"])
        self.E = numpy.float(mgp["E"])
        self.Z = numpy.float(mgp["Z"])
        self.F = numpy.float(mgp["F"])
        self.XY = numpy.float(mgp["XY"])
        self.XZ = numpy.float(mgp["XZ"])
        self.YZ = numpy.float(mgp["YZ"])

        self.A2 = numpy.float(mgp["A2"])
        self.B2 = numpy.float(mgp["B2"])
        self.C2 = numpy.float(mgp["C2"])
        self.UV2 = numpy.float(mgp["UV2"])
        self.UW2 = numpy.float(mgp["UW2"])
        self.VW2 = numpy.float(mgp["VW2"])
        self.D2 = numpy.float(mgp["D2"])
        self.E2 = numpy.float(mgp["E2"])
        self.F2 = numpy.float(mgp["F2"])
        self.XY2 = numpy.float(mgp["XY2"])
        self.XZ2 = numpy.float(mgp["XZ2"])
        self.YZ2 = numpy.float(mgp["YZ2"])


        self.color = [mgp["Red"],mgp["Green"],mgp["Blue"]]

        self.coeff_all = [ [mgp['field_all_M'],mgp['field_all_S']],[mgp['field_pm_M'],mgp['field_pm_S']],[mgp['field_dist_M'],mgp['field_dist_S']],[mgp['field_rv_M'],mgp['field_rv_S']],[mgp['field_pmdist_M'],mgp['field_pmdist_S']],[mgp['field_pmrv_M'],mgp['field_pmrv_M']],[mgp['field_distrv_M'],mgp['field_distrv_S']]]
        self.coeff_young = [ [mgp['young_all_M'],mgp['young_all_S']],[mgp['young_pm_M'],mgp['young_pm_S']],[mgp['young_dist_M'],mgp['young_dist_S']],[mgp['young_rv_M'],mgp['young_rv_S']],[mgp['young_pmdist_M'],mgp['young_pmdist_S']],[mgp['young_pmrv_M'],mgp['young_pmrv_S']],[mgp['young_distrv_M'],mgp['young_distrv_S']]]

# read in ellipse parameters for associations

filename = argv[1]
file = open('Moving_Group_all.csv','rb')
readtable = ascii.get_reader(Reader=ascii.Basic)
readtable.header.splitter.delimiter = ','
readtable.data.splitter.delimiter = ','
readtable.header.start_line = 0
readtable.data.start_line = 1
groups = readtable.read(file)
file.close()

outfile = open(outname,'wb')
moving_groups = []

moving_groups.append(Mgp(groups[0]))
moving_groups.append(Mgp(groups[1]))
moving_groups.append(Mgp(groups[2]))
moving_groups.append(Mgp(groups[3]))
moving_groups.append(Mgp(groups[4]))
moving_groups.append(Mgp(groups[5]))
moving_groups.append(Mgp(groups[6]))
moving_groups.append(Mgp(groups[7]))
moving_groups.append(Mgp(groups[8]))
moving_groups.append(Mgp(groups[9]))
moving_groups.append(Mgp(groups[10]))
moving_groups.append(Mgp(groups[11]))
moving_groups.append(Mgp(groups[12]))
moving_groups.append(Mgp(groups[13]))

filename = argv[1]
file = open(filename,'rb')
readtable = ascii.get_reader(Reader=ascii.Basic)
readtable.header.splitter.delimiter = ','
readtable.data.splitter.delimiter = ','
readtable.header.start_line = 1
readtable.data.start_line = 3
star = readtable.read(file)
file.close()

for i in numpy.arange(0,len(star)):
    distlen = 1
    rvlen = 1
    monte = 10000000
    
    name = star[i]['Name']
    try:
        ra = float(star[i]['RA'])
        dec = float(star[i]['DEC'])
    except ValueError:
        try:
            ra = astrometry.ten((float(star[i]['Rah']),float(star[i]['Ram']),float(star[i]['Ras'])))*15.
            dec = astrometry.ten((numpy.abs(float(star[i]['DECd'])),float(star[i]['DECm']),float(star[i]['DECs'])))
            if star[i]['DECf'] == '-':
                dec = dec * -1.0
        except ValueError:
            ra = astrometry.ten((float(star[i]['Rah']),float(star[i]['Ram']),float(star[i]['Ras'])))*15.
            dec = astrometry.ten((numpy.abs(float(star[i]['DECd'])),float(star[i]['DECm'])))
            if star[i]['DECf'] == '-':
                dec = dec * -1.0

    try:
        era = numpy.float(star[i]['eRA'])/1000./3600.
        edec = numpy.float(star[i]['eDEC'])/1000./3600.
    except (IndexError,ValueError):
        era = 1.0/3600.
        edec = 1.0/3600.

    pmexists = 0
    try:
        pmra = numpy.float(star[i]['pmRA'])/1000.
        epmra = numpy.float(star[i]['epmRA'])/1000.
        pmdec = numpy.float(star[i]['pmDEC'])/1000.
        epmdec = numpy.float(star[i]['epmDEC'])/1000.
        pm,epm,pa,epa = astrometry.pmjoin(pmra,epmra,pmdec,epmdec)
        pmexists = 1
    except (IndexError,ValueError):
        try:
            pmra = numpy.float(star[i]['pmRA'])/1000.
            epmra = 0.01
            pmdec = numpy.float(star[i]['pmDEC'])/1000.
            epmdec = 0.01
            pm,epm,pa,epa = astrometry.pmjoin(pmra,epmra,pmdec,epmdec)
            pmexists = 1
        except (IndexError,ValueError):
            pmexists = 0

    distexists = 0
    rvexists = 0
    try:
        pi = numpy.float(star[i]['pi'])/1000.
        pi = [pi]
        epi = numpy.float(star[i]['epi'])/1000.
        epi = [epi]
        distlen = 1
        distexists = 1
    except (IndexError,ValueError):
        pi = 1/numpy.arange(1.0,150.0,2.0)
        epi = numpy.ones(len(pi))*0.001
        distlen = len(pi)
        monte = monte/100.
        distexists = 0
    try:
        rv = numpy.float(star[i]['rv'])
        rv = [rv]
        erv = numpy.float(star[i]['erv'])
        erv = [erv]
        rvlen = 1
        rvexists = 1
    except (IndexError,ValueError):
        rv = numpy.arange(-80,80,2)
        erv = [1]*len(rv)
        rvlen = len(rv)
        rvexists = 0
        monte = monte/100
        
    print "{8:}  {0:09.5f} {1:7.5f}  {2:+09.5f} {3:7.5f}  {4:+9.5f} {5:8.5f}  {6:+9.5f} {7:8.5f}  {9:5.4f} {10:06.2f}".format(ra,era,dec,edec,pmra,epmra,pmdec,epmdec,name,pm,pa)
    result = numpy.zeros((distlen,rvlen,18))
    #print len(result)

    #print dist,edist
    #print monte,distlen,rvlen

    print distlen,rvlen,result.shape
    for j in xrange(distlen):
        for k in xrange(rvlen):
            tra = ra + numpy.random.randn(monte)*era/numpy.cos(dec*numpy.pi/180)
            tdec = dec + numpy.random.randn(monte)*edec
            tpi = pi[j] + numpy.random.randn(monte)*epi[j]
            tpmra = pmra + numpy.random.randn(monte)*epmra
            tpmdec = pmdec + numpy.random.randn(monte)*epmdec
            trv = rv[k] + numpy.random.randn(monte)*erv[k]

            bu,bv,bw,bx,by,bz = kinematics.gal_uvwxyz(ra=tra,dec=tdec,plx=tpi,pmra=tpmra,pmdec=tpmdec,vrad=trv)

            objXY = ellipse.fitellipse2d(bx,by)
            objXZ = ellipse.fitellipse2d(bx,bz)
            objYZ = ellipse.fitellipse2d(by,bz)
            
            objUV = ellipse.fitellipse2d(bu,bv)
            objUW = ellipse.fitellipse2d(bu,bw)
            objVW = ellipse.fitellipse2d(bv,bw)
            
            # if this is useful, there will only be one of these
            U = numpy.mean(bu)
            V = numpy.mean(bv)
            W = numpy.mean(bw)
            X = numpy.mean(bx)
            Y = numpy.mean(by)
            Z = numpy.mean(bz)

            eU = numpy.std(bu,ddof=1)
            eV = numpy.std(bv,ddof=1)
            eW = numpy.std(bw,ddof=1)
            eX = numpy.std(bx,ddof=1)
            eY = numpy.std(by,ddof=1)
            eZ = numpy.std(bz,ddof=1)

            result[j][k][0] = objUV['x'] 
            result[j][k][1] = objUV['a'] 
            result[j][k][2] = objVW['x'] 
            result[j][k][3] = objVW['a'] 
            result[j][k][4] = objUW['y'] 
            result[j][k][5] = objUW['b'] 
            result[j][k][6] = objUV['xy']
            result[j][k][7] = objUW['xy']
            result[j][k][8] = objVW['xy']
            result[j][k][9] = objXY['x'] 
            result[j][k][10] = objXY['a']  
            result[j][k][11] = objYZ['x']  
            result[j][k][12] = objYZ['a']  
            result[j][k][13] = objXZ['y']  
            result[j][k][14] = objXZ['b']  
            result[j][k][15] = objXY['xy'] 
            result[j][k][16] = objXZ['xy'] 
            result[j][k][17] = objYZ['xy'] 


    #    print result[0][:][1]
    result = zip(*result)
    if distlen == 1 and rvlen == 1:
        print 'U={0:+6.3f} {1:6.3f} V={2:+6.3f} {3:6.3f} W={4:+6.3f} {5:6.3f}'.format(U,eU,V,eV,W,eW)
        outfile.write( 'U={0:+6.3f} {1:6.3f} V={2:+6.3f} {3:6.3f} W={4:+6.3f} {5:6.3f}\n'.format(U,eU,V,eV,W,eW))
    if distlen == 1:
        print 'X={0:+6.3f} {1:6.3f} Y={2:+6.3f} {3:6.3f} Z={4:+6.3f} {5:6.3f}'.format(X,eX,Y,eY,Z,eZ)
        outfile.write('X={0:+6.3f} {1:6.3f} Y={2:+6.3f} {3:6.3f} Z={4:+6.3f} {5:6.3f}\n'.format(X,eX,Y,eY,Z,eZ))
    


      
##############################
# Now find association matches
##############################


# figure for plotting
    fig = pyplot.figure()
    if XYZ:
        fig.set_size_inches(12,6)
        #fig.suptitle(name+' UVW matches', fontproperties=prop, size='x-large')
        ax1 = fig.add_subplot(231, aspect='equal')
        ax2 = fig.add_subplot(232, aspect='equal')
        ax3 = fig.add_subplot(233, aspect='equal')
        ax4 = fig.add_subplot(234, aspect='equal')
        ax5 = fig.add_subplot(235, aspect='equal')
        ax6 = fig.add_subplot(236, aspect='equal')
    else:
        fig.set_size_inches(12,3)
        #fig.suptitle(name+' UVW matches', fontproperties=prop, size='x-large')
        # UV plane
        ax1 = fig.add_subplot(131, aspect='equal')
        ax2 = fig.add_subplot(132, aspect='equal')
        ax3 = fig.add_subplot(133, aspect='equal')
    ax1.set_xlabel('U (km s$^{-1}$)')#,fontproperties=prop)
    ax1.set_ylabel('V (km s$^{-1}$)')#,fontproperties=prop)
    ax1.set_xlim(-25,20)
    ax1.set_ylim(-30,5)
    ax1.set_xticklabels(numpy.asarray(ax1.get_xticks(),dtype=numpy.int))#,fontproperties=prop)
    ax1.set_yticklabels(numpy.asarray(ax1.get_yticks(),dtype=numpy.int))#,fontproperties=prop)
    ax2.set_xlabel('U (km s$^{-1}$)')#,fontproperties=prop)
    ax2.set_ylabel('W (km s$^{-1}$)')#,fontproperties=prop)
    ax2.set_xlim(-25,20)
    ax2.set_ylim(-25,10)
    ax2.set_xticklabels(numpy.asarray(ax2.get_xticks(),dtype=numpy.int))#,fontproperties=prop)
    ax2.set_yticklabels(numpy.asarray(ax2.get_yticks(),dtype=numpy.int))#,fontproperties=prop)
    ax3.set_xlabel('V (km s$^{-1}$)')#,fontproperties=prop)
    ax3.set_ylabel('W (km s$^{-1}$)')#,fontproperties=prop)
    ax3.set_xlim(-35,10)
    ax3.set_ylim(-20,15)
    ax3.set_xticklabels(numpy.asarray(ax3.get_xticks(),dtype=numpy.int))#,fontproperties=prop)
    ax3.set_yticklabels(numpy.asarray(ax3.get_yticks(),dtype=numpy.int))#,fontproperties=prop)
    if XYZ:
        ax4.set_xlabel('X (pc)')#,fontproperties=prop)
        ax4.set_ylabel('Y (pc)')#,fontproperties=prop)
        ax4.set_xlim(-150,100)
        ax4.set_ylim(-150,50)
        ax4.set_xticklabels(numpy.asarray(ax4.get_xticks(),dtype=numpy.int))#,fontproperties=prop)
        ax4.set_yticklabels(numpy.asarray(ax4.get_yticks(),dtype=numpy.int))#,fontproperties=prop)
        ax5.set_xlabel('X (pc)')#,fontproperties=prop)
        ax5.set_ylabel('Z (pc)')#,fontproperties=prop)
        ax5.set_xlim(-150,100)
        ax5.set_ylim(-100,100)
        ax5.set_xticklabels(numpy.asarray(ax5.get_xticks(),dtype=numpy.int))#,fontproperties=prop)
        ax5.set_yticklabels(numpy.asarray(ax5.get_yticks(),dtype=numpy.int))#,fontproperties=prop)
        ax6.set_xlabel('Y (pc)')#,fontproperties=prop)
        ax6.set_ylabel('Z (pc)')#,fontproperties=prop)
        ax6.set_xlim(-150,100)
        ax6.set_ylim(-100,100)
        ax6.set_xticklabels(numpy.asarray(ax6.get_xticks(),dtype=numpy.int))#,fontproperties=prop)
        ax6.set_yticklabels(numpy.asarray(ax6.get_yticks(),dtype=numpy.int))#,fontproperties=prop)
    
    a = 0
    for i in range(len(moving_groups)):
        assocellipse1 = Ellipse(xy=(moving_groups[i].U,moving_groups[i].V),width=moving_groups[i].A2*2,height=moving_groups[i].B2*2, angle=moving_groups[i].UV2*180/numpy.pi,linewidth=0.5)
        assocellipse2 = Ellipse(xy=(moving_groups[i].U,moving_groups[i].W),width=moving_groups[i].A*2,height=moving_groups[i].C2*2, angle=moving_groups[i].UW2*180/numpy.pi,linewidth=0.5)
        assocellipse3 = Ellipse(xy=(moving_groups[i].V,moving_groups[i].W),width=moving_groups[i].B2*2,height=moving_groups[i].C2*2, angle=moving_groups[i].VW2*180/numpy.pi,linewidth=0.5)
        ax1.add_artist(assocellipse1)
        assocellipse1.set_alpha(.6)
        assocellipse1.set_clip_box(ax1.bbox)
        assocellipse1.set_facecolor(moving_groups[i].color)
        ax2.add_artist(assocellipse2)
        assocellipse2.set_alpha(.6)
        assocellipse2.set_clip_box(ax2.bbox)
        assocellipse2.set_facecolor(moving_groups[i].color)
        ax3.add_artist(assocellipse3)
        assocellipse3.set_alpha(.6)
        assocellipse3.set_clip_box(ax3.bbox)
        assocellipse3.set_facecolor(moving_groups[i].color)	
        if XYZ:
            assocellipse4 = Ellipse(xy=(moving_groups[i].X,moving_groups[i].Y),width=moving_groups[i].D2*2,height=moving_groups[i].E2*2, angle=moving_groups[i].XY2*180/numpy.pi,linewidth=0.5)
            assocellipse5 = Ellipse(xy=(moving_groups[i].X,moving_groups[i].Z),width=moving_groups[i].D2*2,height=moving_groups[i].F2*2, angle=moving_groups[i].XZ2*180/numpy.pi,linewidth=0.5)
            assocellipse6 = Ellipse(xy=(moving_groups[i].Y,moving_groups[i].Z),width=moving_groups[i].E2*2,height=moving_groups[i].F2*2, angle=moving_groups[i].YZ2*180/numpy.pi,linewidth=0.5)
            ax4.add_artist(assocellipse4)
            assocellipse4.set_alpha(.6)
            assocellipse4.set_clip_box(ax4.bbox)
            assocellipse4.set_facecolor(moving_groups[i].color)
            ax5.add_artist(assocellipse5)
            assocellipse5.set_alpha(.6)
            assocellipse5.set_clip_box(ax5.bbox)
            assocellipse5.set_facecolor(moving_groups[i].color)
            ax6.add_artist(assocellipse6)
            assocellipse6.set_alpha(.6)
            assocellipse6.set_clip_box(ax6.bbox)
            assocellipse6.set_facecolor(moving_groups[i].color)	
            
        ax2.text(5,8+a*-2.5,moving_groups[i].name,color=moving_groups[i].color)#,fontproperties=prop)
        a=a+1

    for n in xrange(rvlen):
        for p in xrange(distlen):
            print n,p
            pts1 = Ellipse(xy = (result[n][p][0],result[n][p][2]), width=result[n][p][1]*2,height=result[n][p][3]*2, angle=result[n][p][6]*180/numpy.pi,linewidth=3)
            pts2 = Ellipse(xy = (result[n][p][0],result[n][p][4]), width=result[n][p][1]*2,height=result[n][p][5]*2, angle=result[n][p][7]*180/numpy.pi,linewidth=3)
            pts3 = Ellipse(xy = (result[n][p][2],result[n][p][4]), width=result[n][p][3]*2,height=result[n][p][5]*2, angle=result[n][p][8]*180/numpy.pi,linewidth=3)
            ax1.add_artist(pts1)
            pts1.set_clip_box(ax1.bbox)
            pts1.set_alpha(.60)
            pts1.set_facecolor('#000000')
            pts1.set_edgecolor('#FF0000')
            ax2.add_artist(pts2)
            pts2.set_clip_box(ax2.bbox)
            pts2.set_alpha(.60)
            pts2.set_facecolor('#000000')
            pts2.set_edgecolor('#FF0000')
            ax3.add_artist(pts3)
            pts3.set_clip_box(ax3.bbox)
            pts3.set_alpha(.60)
            pts3.set_facecolor('#000000')
            pts3.set_edgecolor('#FF0000')
            if XYZ:
                pts4 = Ellipse(xy = (result[n][p][9],result[n][p][11]), width=result[n][p][10]*2,height=result[n][p][12]*2, angle=result[n][p][15]*180/numpy.pi,linewidth=3)
                pts5 = Ellipse(xy = (result[n][p][9],result[n][p][13]), width=result[n][p][10]*2,height=result[n][p][14]*2, angle=result[n][p][16]*180/numpy.pi,linewidth=3)
                pts6 = Ellipse(xy = (result[n][p][11],result[n][p][13]), width=result[n][p][12]*2,height=result[n][p][14]*2, angle=result[n][p][17]*180/numpy.pi,linewidth=3)
                ax4.add_artist(pts4)
                pts4.set_clip_box(ax4.bbox)
                pts4.set_alpha(.60)
                pts4.set_facecolor('#000000')
                pts4.set_edgecolor('#FF0000')
                ax5.add_artist(pts5)
                pts5.set_clip_box(ax5.bbox)
                pts5.set_alpha(.60)
                pts5.set_facecolor('#000000')
                pts5.set_edgecolor('#FF0000')
                ax6.add_artist(pts6)
                pts6.set_clip_box(ax6.bbox)
                pts6.set_alpha(.60)
                pts6.set_facecolor('#000000')
                pts6.set_edgecolor('#FF0000')
            
    if ((distlen == 1) & (XYZ)):
        ax4.scatter(result[n][p][9],result[n][p][11],marker='*',s=50,zorder=100,c='#FF0000')
        ax5.scatter(result[n][p][9],result[n][p][13],marker='*',s=50,zorder=100,c='#FF0000')
        ax6.scatter(result[n][p][11],result[n][p][13],marker='*',s=50,zorder=100,c='#FF0000')

    # Isometric 3D plot
    #ax = fig.gca(projection='3d')
    #for i in result[:,:,0]:
    #X,Y,Z = triaxial(result[i,0],result[i,2],result[i,4],result[i,1],result[i,3],result[i,5])
    outname = 'mgp_figures/{0:}.png'.format(name)
    pyplot.tight_layout()
    fig.savefig(outname,dpi=300,transparent=False)
    pyplot.clf()
    pyplot.close()
