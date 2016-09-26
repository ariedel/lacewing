import numpy as np
from matplotlib import pyplot
from matplotlib.patches import Ellipse
import astropy.io.ascii as ascii
from sys import argv
import kinematics
import ellipse
import lacewing

###################################
###################################
### MAIN ROUTINE
###################################
###################################

try:
    infilename = argv[1]
except IndexError: 
    print 'syntax: python lacewing_uvwxyz.py inputfile XYZ outfile'
try:
    outname = argv[2]
except IndexError:
    outname = 'lacewing_uvwxyz.csv'
try:
   XYZ = argv[3]
   if XYZ == "XYZ":
       XYZ = True
except IndexError:
   XYZ = False

moving_groups = lacewing.moving_group_loader()
moving_groups = [moving_groups[x] for x in range(len(moving_groups)) if moving_groups[x].name != "Field"]

name,coord,era,edec,pmra,epmra,pmdec,epmdec,rv,erv,plx,eplx,note = lacewing.csv_loader(infilename)
print name

outfile = open(outname,'wb')

for i in xrange(len(coord)):
    if pmra[i] == None:
        continue
    monte = 1000000
    
    if rv[i] == None:
        prv = np.arange(-70,70,2)
        perv = np.ones_like(prv)*1
        monte = monte / 100.
    else:
        prv = [rv[i]]
        perv = [erv[i]]
    if plx[i] == None:
        pplx = np.arange(0,1.0000,0.01)
        peplx = np.ones_like(pplx)*0.01
        monte = monte / 100.
    else:
        pplx = [plx[i]]
        peplx = [eplx[i]]
    rvlen = len(prv)
    distlen = len(pplx)
    result = np.zeros((distlen,rvlen,18))

    for j in xrange(distlen):
        for k in xrange(rvlen):
            tra = coord[i].ra.degree + np.random.randn(monte)*era[i]/np.cos(coord[i].dec.degree*np.pi/180)
            tdec = coord[i].dec.degree + np.random.randn(monte)*edec[i]
            tplx = pplx[j] + np.random.randn(monte)*peplx[j]
            tpmra = pmra[i] + np.random.randn(monte)*epmra[i]
            tpmdec = pmdec[i] + np.random.randn(monte)*epmdec[i]
            trv = prv[k] + np.random.randn(monte)*perv[k]

            bu,bv,bw,bx,by,bz = kinematics.gal_uvwxyz(ra=tra,dec=tdec,plx=tplx,pmra=tpmra,pmdec=tpmdec,vrad=trv)

            objXY = ellipse.fitellipse2d(bx,by)
            objXZ = ellipse.fitellipse2d(bx,bz)
            objYZ = ellipse.fitellipse2d(by,bz)
            
            objUV = ellipse.fitellipse2d(bu,bv)
            objUW = ellipse.fitellipse2d(bu,bw)
            objVW = ellipse.fitellipse2d(bv,bw)
            
            # if this is useful, there will only be one of these
            U = np.mean(bu)
            V = np.mean(bv)
            W = np.mean(bw)
            X = np.mean(bx)
            Y = np.mean(by)
            Z = np.mean(bz)

            eU = np.std(bu,ddof=1)
            eV = np.std(bv,ddof=1)
            eW = np.std(bw,ddof=1)
            eX = np.std(bx,ddof=1)
            eY = np.std(by,ddof=1)
            eZ = np.std(bz,ddof=1)

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


    result = zip(*result)
    if distlen == 1 and rvlen != 1:
        print '{6:},X=,{0:+6.3f},{1:6.3f},Y=,{2:+6.3f},{3:6.3f},Z=,{4:+6.3f},{5:6.3f}'.format(X,eX,Y,eY,Z,eZ,name[i])
        outfile.write('{6:},X=,{0:+6.3f},{1:6.3f},Y=,{2:+6.3f},{3:6.3f},Z=,{4:+6.3f},{5:6.3f},\n'.format(X,eX,Y,eY,Z,eZ,name[i]))
    if distlen == 1 and rvlen == 1:
        print '{6:},X=,{0:+6.3f},{1:6.3f},Y=,{2:+6.3f},{3:6.3f},Z=,{4:+6.3f},{5:6.3f},'.format(X,eX,Y,eY,Z,eZ,name[i]),
        outfile.write('{6:},X=,{0:+6.3f},{1:6.3f},Y=,{2:+6.3f},{3:6.3f},Z=,{4:+6.3f},{5:6.3f},'.format(X,eX,Y,eY,Z,eZ,name[i]))
        print 'U=,{0:+6.2f},{1:6.2f},V=,{2:+6.2f},{3:6.2f},W=,{4:+6.2f},{5:6.2f}'.format(U,eU,V,eV,W,eW)
        outfile.write('U=,{0:+6.2f},{1:6.2f},V=,{2:+6.2f},{3:6.2f},W=,{4:+6.2f},{5:6.2f},\n'.format(U,eU,V,eV,W,eW))
    
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
    ax1.set_xlabel('U (km s$^{-1}$)')
    ax1.set_ylabel('V (km s$^{-1}$)')
    ax1.set_xlim(-25,20)
    ax1.set_ylim(-30,5)
    ax1.set_xticklabels(np.asarray(ax1.get_xticks(),dtype=np.int))
    ax1.set_yticklabels(np.asarray(ax1.get_yticks(),dtype=np.int))
    ax2.set_xlabel('U (km s$^{-1}$)')
    ax2.set_ylabel('W (km s$^{-1}$)')
    ax2.set_xlim(-25,20)
    ax2.set_ylim(-25,10)
    ax2.set_xticklabels(np.asarray(ax2.get_xticks(),dtype=np.int))
    ax2.set_yticklabels(np.asarray(ax2.get_yticks(),dtype=np.int))
    ax3.set_xlabel('V (km s$^{-1}$)')
    ax3.set_ylabel('W (km s$^{-1}$)')
    ax3.set_xlim(-35,10)
    ax3.set_ylim(-20,15)
    ax3.set_xticklabels(np.asarray(ax3.get_xticks(),dtype=np.int))
    ax3.set_yticklabels(np.asarray(ax3.get_yticks(),dtype=np.int))
    if XYZ:
        ax4.set_xlabel('X (pc)')
        ax4.set_ylabel('Y (pc)')
        ax4.set_xlim(-150,100)
        ax4.set_ylim(-150,50)
        ax4.set_xticklabels(np.asarray(ax4.get_xticks(),dtype=np.int))
        ax4.set_yticklabels(np.asarray(ax4.get_yticks(),dtype=np.int))
        ax5.set_xlabel('X (pc)')
        ax5.set_ylabel('Z (pc)')
        ax5.set_xlim(-150,100)
        ax5.set_ylim(-100,100)
        ax5.set_xticklabels(np.asarray(ax5.get_xticks(),dtype=np.int))
        ax5.set_yticklabels(np.asarray(ax5.get_yticks(),dtype=np.int))
        ax6.set_xlabel('Y (pc)')
        ax6.set_ylabel('Z (pc)')
        ax6.set_xlim(-150,100)
        ax6.set_ylim(-100,100)
        ax6.set_xticklabels(np.asarray(ax6.get_xticks(),dtype=np.int))
        ax6.set_yticklabels(np.asarray(ax6.get_yticks(),dtype=np.int))
    
    a = 0
    for m in range(len(moving_groups)):
        assocellipse1 = Ellipse(xy=(moving_groups[m].U,moving_groups[m].V),width=moving_groups[m].A2*2,height=moving_groups[m].B2*2, angle=moving_groups[m].UV2*180/np.pi,linewidth=0.5)
        assocellipse2 = Ellipse(xy=(moving_groups[m].U,moving_groups[m].W),width=moving_groups[m].A2*2,height=moving_groups[m].C2*2, angle=moving_groups[m].UW2*180/np.pi,linewidth=0.5)
        assocellipse3 = Ellipse(xy=(moving_groups[m].V,moving_groups[m].W),width=moving_groups[m].B2*2,height=moving_groups[m].C2*2, angle=moving_groups[m].VW2*180/np.pi,linewidth=0.5)
        ax1.add_artist(assocellipse1)
        assocellipse1.set_alpha(.6)
        assocellipse1.set_clip_box(ax1.bbox)
        assocellipse1.set_facecolor(moving_groups[m].color)
        ax2.add_artist(assocellipse2)
        assocellipse2.set_alpha(.6)
        assocellipse2.set_clip_box(ax2.bbox)
        assocellipse2.set_facecolor(moving_groups[m].color)
        ax3.add_artist(assocellipse3)
        assocellipse3.set_alpha(.6)
        assocellipse3.set_clip_box(ax3.bbox)
        assocellipse3.set_facecolor(moving_groups[m].color)	
        if XYZ:
            assocellipse4 = Ellipse(xy=(moving_groups[m].X,moving_groups[m].Y),width=moving_groups[m].D2*2,height=moving_groups[m].E2*2, angle=moving_groups[m].XY2*180/np.pi,linewidth=0.5)
            assocellipse5 = Ellipse(xy=(moving_groups[m].X,moving_groups[m].Z),width=moving_groups[m].D2*2,height=moving_groups[m].F2*2, angle=moving_groups[m].XZ2*180/np.pi,linewidth=0.5)
            assocellipse6 = Ellipse(xy=(moving_groups[m].Y,moving_groups[m].Z),width=moving_groups[m].E2*2,height=moving_groups[m].F2*2, angle=moving_groups[m].YZ2*180/np.pi,linewidth=0.5)
            ax4.add_artist(assocellipse4)
            assocellipse4.set_alpha(.6)
            assocellipse4.set_clip_box(ax4.bbox)
            assocellipse4.set_facecolor(moving_groups[m].color)
            ax5.add_artist(assocellipse5)
            assocellipse5.set_alpha(.6)
            assocellipse5.set_clip_box(ax5.bbox)
            assocellipse5.set_facecolor(moving_groups[m].color)
            ax6.add_artist(assocellipse6)
            assocellipse6.set_alpha(.6)
            assocellipse6.set_clip_box(ax6.bbox)
            assocellipse6.set_facecolor(moving_groups[m].color)	
            
        ax2.text(5,8+a*-1.7,moving_groups[m].name,color=moving_groups[m].color,fontsize="small")
        a=a+1

    for n in xrange(rvlen):
        for p in xrange(distlen):
            pts1 = Ellipse(xy = (result[n][p][0],result[n][p][2]), width=result[n][p][1]*2,height=result[n][p][3]*2, angle=result[n][p][6]*180/np.pi,linewidth=3)
            pts2 = Ellipse(xy = (result[n][p][0],result[n][p][4]), width=result[n][p][1]*2,height=result[n][p][5]*2, angle=result[n][p][7]*180/np.pi,linewidth=3)
            pts3 = Ellipse(xy = (result[n][p][2],result[n][p][4]), width=result[n][p][3]*2,height=result[n][p][5]*2, angle=result[n][p][8]*180/np.pi,linewidth=3)
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
                pts4 = Ellipse(xy = (result[n][p][9],result[n][p][11]), width=result[n][p][10]*2,height=result[n][p][12]*2, angle=result[n][p][15]*180/np.pi,linewidth=3)
                pts5 = Ellipse(xy = (result[n][p][9],result[n][p][13]), width=result[n][p][10]*2,height=result[n][p][14]*2, angle=result[n][p][16]*180/np.pi,linewidth=3)
                pts6 = Ellipse(xy = (result[n][p][11],result[n][p][13]), width=result[n][p][12]*2,height=result[n][p][14]*2, angle=result[n][p][17]*180/np.pi,linewidth=3)
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

    outname = '{0:}.png'.format(name[i])
    pyplot.tight_layout()
    fig.savefig(outname,dpi=300,transparent=False)
    pyplot.clf()
    pyplot.close()
