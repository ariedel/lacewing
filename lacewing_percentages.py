import numpy
from scipy.special import erfc
from matplotlib import pyplot as plt
from sys import argv
from scipy.optimize import curve_fit
import kinematics
import lacewing
import os

#########################################################
#########################################################
### MAIN ROUTINE
### ARR 2015-10-16
### 1.3: Now compatible with LACEwING v1.3, and uses the 
###      lacewing.moving_group_loader() function
#########################################################
#########################################################

def gaus(x,a,sigma):
    return a*numpy.exp(-(x)**2/(2*sigma**2))

def gauscdf(x,m,sigma):
    return 0.5*(erfc((x-m)/(sigma*numpy.sqrt(2))))

stop = 24
step = 0.1

sigrange = numpy.arange(0.0,stop,step)
acc_all = numpy.zeros_like(sigrange)
acc_pm = numpy.zeros_like(sigrange)
acc_dist = numpy.zeros_like(sigrange)
acc_rv = numpy.zeros_like(sigrange)
acc_pmdist = numpy.zeros_like(sigrange)
acc_pmrv = numpy.zeros_like(sigrange)
acc_distrv = numpy.zeros_like(sigrange)
total_all = numpy.zeros_like(sigrange)
total_pm = numpy.zeros_like(sigrange)
total_dist = numpy.zeros_like(sigrange)
total_rv = numpy.zeros_like(sigrange)
total_pmdist = numpy.zeros_like(sigrange)
total_pmrv = numpy.zeros_like(sigrange)
total_distrv = numpy.zeros_like(sigrange)
good_all = numpy.zeros((15,len(sigrange)))
good_pm = numpy.zeros((15,len(sigrange)))
good_dist = numpy.zeros((15,len(sigrange)))
good_rv = numpy.zeros((15,len(sigrange)))
good_pmdist = numpy.zeros((15,len(sigrange)))
good_pmrv = numpy.zeros((15,len(sigrange)))
good_distrv = numpy.zeros((15,len(sigrange)))

moving_groups = lacewing.moving_group_loader()

groups = []
for i in xrange(len(moving_groups)):
    groups.append(moving_groups[i].name.replace(' ','_'))
groups = numpy.asarray(groups)

youth=""
if argv[2] == 'young':
    youth=".youngonly"

if not os.path.exists('montecarlo'):
    os.mkdir('montecarlo')

i = 0
# The input file contains every star as matched to one of the moving groups.
# This loop reads in every line in the file, one by one, and builds histograms
# Each line contains the goodness-of-fit values for the 7 different data cases, relative to one moving group.
# For each of the 7 data cases, there are arrays: One for total number of points
#  with that particular binned goodness-of-fit score, and 14 to separate out objects that were "really" 
#  members of each of the moving groups, plus one for objects that were "really" field stars.
with open(argv[1], 'rb') as f:
    for line in f:
        if i % 100000 == 0:
            print i
        entry = line.split(",")

        if argv[2] == 'young':
            if entry[9] == "Field":
                continue
        
        #print entry[9],numpy.where(entry[9] == groups)
        grp = numpy.where(entry[9] == groups)[0][0]
        #print grp

        location = float(entry[1])/step
        if location < stop/step:
            total_all[location] += 1
            good_all[grp][location] += 1
                    
        location = float(entry[2])/step
        if location < stop/step:
            total_pm[location] += 1
            good_pm[grp][location] += 1

        location = float(entry[3])/step
        if location < stop/step:
            total_dist[location] += 1
            good_dist[grp][location] += 1

        location = float(entry[4])/step
        if location < stop/step:
            total_rv[location] += 1
            good_rv[grp][location] += 1

        location = float(entry[5])/step
        if location < stop/step:
            total_pmdist[location] += 1
            good_pmdist[grp][location] += 1

        location = float(entry[6])/step
        if location < stop/step:
            total_pmrv[location] += 1
            good_pmrv[grp][location] += 1

        location = float(entry[7])/step
        if location < stop/step:
            total_distrv[location] += 1
            good_distrv[grp][location] += 1

        i+= 1

# Now identify the moving group we're interested in...
print entry[0]
real = numpy.where(entry[0].replace(' ','_') == groups)[0][0]
print real

groupname = entry[0].replace(' ','_')

outfile = open('montecarlo/{0:}{1:}.values'.format(groupname,youth),'wb')

outfile.write('sigval,tot_all,good_all,tot_pmdist,good_pmdist,tot_pmrv,good_pmrv,tot_distrv,good_distrv,tot_pm,good_pm,tot_dist,good_dist,tot_rv,good_rv\n')

for m in range(len(sigrange)):

    print sigrange[m], total_all[m], good_all[real][m], total_pmdist[m], good_pmdist[real][m], total_pmrv[m], good_pmrv[real][m], total_distrv[m],good_distrv[real][m],total_pm[m],good_pm[real][m],total_dist[m],good_dist[real][m],total_rv[m],good_rv[real][m]
    outfile.write('{0:},{1:},{2:},{3:},{4:},{5:},{6:},{7:},{8:},{9:},{10:},{11:},{12:},{13:},{14:}\n'.format(sigrange[m], total_all[m], good_all[real][m], total_pmdist[m], good_pmdist[real][m], total_pmrv[m], good_pmrv[real][m], total_distrv[m],good_distrv[real][m],total_pm[m],good_pm[real][m],total_dist[m],good_dist[real][m],total_rv[m],good_rv[real][m]))

outfile.close()

# First plot: Percentage of total stars with that score that are "really" members, as a function of the goodness-of-fit value

fig1 = plt.figure(figsize=(6,4))
ax = fig1.add_subplot(111)

others = range(len(groups))
others.remove(real)
print others

# Create arrays of the fraction of real members, as a function of goodness-of-fit
fraction_all = numpy.asarray(good_all[real]/total_all)
fraction_pmdist = numpy.asarray(good_pmdist[real]/total_pmdist)
fraction_pmrv = numpy.asarray(good_pmrv[real]/total_pmrv)
fraction_distrv = numpy.asarray(good_distrv[real]/total_distrv)
fraction_pm = numpy.asarray(good_pm[real]/total_pm)
fraction_dist = numpy.asarray(good_dist[real]/total_dist)
fraction_rv = numpy.asarray(good_rv[real]/total_rv)

pltall = ax.plot(sigrange,fraction_all*100,linewidth=4,color='#000000',label="ALL",drawstyle='steps-mid')
pltpmrv = ax.plot(sigrange,fraction_pmrv*100,linewidth=2,marker='s',color='#FF00FF',label="PM+RV",drawstyle='steps-mid')
pltdistrv = ax.plot(sigrange,fraction_distrv*100,linewidth=2,marker='*',color='#00AAFF',label="DIST+RV",drawstyle='steps-mid')
pltpmdist = ax.plot(sigrange,fraction_pmdist*100,linewidth=2,marker='h',color='#FFAA00',label="PM+DIST",drawstyle='steps-mid')
pltrv = ax.plot(sigrange,fraction_rv*100,linewidth=1,marker='o',color='#0000FF',label="RV",drawstyle='steps-mid')
pltpm = ax.plot(sigrange,fraction_pm*100,linewidth=1,marker='v',color='#FF0000',label="PM",drawstyle='steps-mid')
pltdist = ax.plot(sigrange,fraction_dist*100,linewidth=1,marker='+',color='#00AA00',label="DIST",drawstyle='steps-mid')

plt.legend(loc='upper right',numpoints=1)

cumulative_all = good_all[real]
cumulative_pm = good_pm[real]
cumulative_dist = good_dist[real]
cumulative_rv = good_rv[real]
cumulative_pmdist = good_pmdist[real]
cumulative_pmrv = good_pmrv[real]
cumulative_distrv = good_distrv[real]

##Uncomment this to plot the members of other groups as well, rather than just the fraction that are members of the "real" group.
#for i in others:
#    cumulative_all  +=  good_all[i]
#    cumulative_pm   +=  good_pm[i]
#    cumulative_dist +=  good_dist[i]
#    cumulative_rv   +=  good_rv[i]
#    cumulative_pmdist += good_pmdist[i]
#    cumulative_pmrv +=  good_pmrv[i]
#    cumulative_distrv += good_distrv[i]
#
#    pltall = axall.plot(sigrange,numpy.asarray(cumulative_all,dtype=numpy.float64)/numpy.asarray(total_all,dtype=numpy.float64))
#    pltpmrv = axpmrv.plot(sigrange,numpy.asarray(cumulative_pmrv,dtype=numpy.float64)/numpy.asarray(total_pmrv,dtype=numpy.float64))
#    pltdistrv = axdistrv.plot(sigrange,numpy.asarray(cumulative_distrv,dtype=numpy.float64)/numpy.asarray(total_distrv,dtype=numpy.float64))
#    pltpmdist = axpmdist.plot(sigrange,numpy.asarray(cumulative_pmdist,dtype=numpy.float64)/numpy.asarray(total_pmdist,dtype=numpy.float64))
#    pltrv = axrv.plot(sigrange,numpy.asarray(cumulative_rv,dtype=numpy.float64)/numpy.asarray(total_rv,dtype=numpy.float64))
#    pltpm = axpm.plot(sigrange,numpy.asarray(cumulative_pm,dtype=numpy.float64)/numpy.asarray(total_pm,dtype=numpy.float64))
#    pltdist = axdist.plot(sigrange,numpy.asarray(cumulative_dist,dtype=numpy.float64)/numpy.asarray(total_dist,dtype=numpy.float64))
#    axall.text(0,cumulative_all[0],groups[i])
#    axpm.text(0,cumulative_pm[0],groups[i])
#    axdist.text(0,cumulative_dist[0],groups[i])
#    axrv.text(0,cumulative_rv[0],groups[i])
#    axpmdist.text(0,cumulative_pmdist[0],groups[i])
#    axpmrv.text(0,cumulative_pmrv[0],groups[i])
#    axdistrv.text(0,cumulative_distrv[0],groups[i])

#for axis in [axall,axpm,axdist,axrv,axpmdist,axpmrv,axdistrv]:
#    axis.set_ylim((0,1.1))
#    axis.set_xlim((0,4.0))
ax.set_xlim(0,3.0)

#adjustprops = dict(left=0.2,bottom=0.2,right=0.8,top=0.8,wspace=0.2,hspace=0.2)
ax.set_xlabel('Goodness of Fit')
ax.set_ylabel('Percentage of "real" members')
plt.tight_layout()

# Create the curve fits to the membership fractions above, and write them to a file
#  Open the file
outfile = open('montecarlo/{0:}{1:}.percentages'.format(groupname,youth),'wb')

fractionlist = numpy.asarray([fraction_pm,fraction_dist,fraction_rv,fraction_pmdist,fraction_pmrv,fraction_distrv,fraction_all])
#fractionlist = numpy.asarray([fraction_pm,fraction_dist,fraction_rv,fraction_pmdist,fraction_pmrv,fraction_distrv,fraction_all])

#axislist = [axpm,axdist,axrv,axpmdist,axpmrv,axdistrv,axall]

mean_param = []
sig_param = []
fitcolors = ['r:','g:','b:','y:','m:','c:','k:']

for a in range(len(fractionlist)): 
#for combination in [fraction_all,fraction_pa,fraction_dist,fraction_rv,fraction_padist,fraction_parv,fraction_distrv]:
    sift1 = numpy.where(numpy.asarray(numpy.isnan(fractionlist[a])))
    #print sift1

    if len(sift1) > 0:
        fractionlist[a][sift1] = 1

    sift9 = numpy.where(sigrange < 0)[0]

    if len(sift9) > 0:
        sigcut = sigrange[sift9] 
        percentage = fractionlist[a][sift9]
    else:
        sigcut=sigrange
        percentage=fractionlist[a]

    #print percentage
    
    sigma = 1
    popt,pcov = curve_fit(gauscdf,sigcut,percentage,p0=[0.5,sigma],maxfev=10000000)

    outfile.write('{0:},{1:},'.format(popt[0],popt[1]))
    mean_param.append(popt[0])
    sig_param.append(popt[1])

    ax.plot(sigrange,gauscdf(sigrange,*popt)*100,fitcolors[a],label='fit')


plt.savefig('montecarlo/{0:}{1:}.eps'.format(groupname,youth))
plt.clf()
plt.close()

# plot 2: cumulative number of "real" moving group members recovered, as a function of maximum goodness-of-fit value.

cumulative_all = numpy.cumsum(good_all[real])
cumulative_pmdist = numpy.cumsum(good_pmdist[real])
cumulative_pmrv = numpy.cumsum(good_pmrv[real])
cumulative_distrv = numpy.cumsum(good_distrv[real])
cumulative_pm = numpy.cumsum(good_pm[real])
cumulative_dist = numpy.cumsum(good_dist[real])
cumulative_rv = numpy.cumsum(good_rv[real])

cumfalse_all = numpy.cumsum(total_all-good_all[real])
cumfalse_pmdist = numpy.cumsum(total_pmdist-good_pmdist[real])
cumfalse_pmrv = numpy.cumsum(total_pmrv-good_pmrv[real])
cumfalse_distrv = numpy.cumsum(total_distrv-good_distrv[real])
cumfalse_pm = numpy.cumsum(total_pm-good_pm[real])
cumfalse_dist = numpy.cumsum(total_dist-good_dist[real])
cumfalse_rv = numpy.cumsum(total_rv-good_rv[real])

cumtotal_all = numpy.cumsum(total_all)
cumtotal_pmdist = numpy.cumsum(total_pmdist)
cumtotal_pmrv = numpy.cumsum(total_pmrv)
cumtotal_distrv = numpy.cumsum(total_distrv)
cumtotal_pm = numpy.cumsum(total_pm)
cumtotal_dist = numpy.cumsum(total_dist)
cumtotal_rv = numpy.cumsum(total_rv)

for x in range(len(cumfalse_all)):
    print sigrange[x],cumulative_all[x],cumfalse_all[x],cumtotal_all[x]

fig2 = plt.figure(figsize=(6,4))
ax2 = fig2.add_subplot(111)
pltall = ax2.plot(sigrange,cumulative_all/cumulative_all[-1]*100,linewidth=4,color = '#000000',label="ALL")
pltpmrv = ax2.plot(sigrange,cumulative_pmrv/cumulative_pmrv[-1]*100,linewidth=2,marker='s', color='#FF00FF',label="PM+RV")
pltdistrv = ax2.plot(sigrange,cumulative_distrv/cumulative_distrv[-1]*100,linewidth=2,marker='*',color='#00AAFF',label="DIST+RV")
pltpmdist = ax2.plot(sigrange,cumulative_pmdist/cumulative_pmdist[-1]*100,linewidth=2,marker='h',color='#FFAA00',label="PM+DIST")
pltrv = ax2.plot(sigrange,cumulative_rv/cumulative_rv[-1]*100,linewidth=1,marker='o',color='#0000FF',label="RV")
pltpm = ax2.plot(sigrange,cumulative_pm/cumulative_pm[-1]*100,linewidth=1,marker='v',color='#FF0000',label="PM")
pltdist = ax2.plot(sigrange,cumulative_dist/cumulative_dist[-1]*100,linewidth=1,marker='+',color='#00AA00',label="DIST")

ax2.set_ylim((0,100.0))
ax2.set_xlim((0,3.0))
ax2.set_xlabel('Maximum Goodness of Fit Parameter')
ax2.set_ylabel('Fraction of actual members recovered')

plt.legend(loc='lower right',numpoints=1)

#plt.title('{0:}'.format(entry[0]))

plt.savefig('montecarlo/{0:}_cumulative.eps'.format(groupname,youth))
plt.clf()
outfile.close()

# Plot 3: Number of stars recovered as a function of minimum percentage accepted

percent_pm = 100*gauscdf(sigrange,mean_param[0],sig_param[0])
percent_dist = 100*gauscdf(sigrange,mean_param[1],sig_param[1])
percent_rv = 100*gauscdf(sigrange,mean_param[2],sig_param[2])
percent_pmdist = 100*gauscdf(sigrange,mean_param[3],sig_param[3])
percent_pmrv = 100*gauscdf(sigrange,mean_param[4],sig_param[4])
percent_distrv = 100*gauscdf(sigrange,mean_param[5],sig_param[5])
percent_all = 100*gauscdf(sigrange,mean_param[6],sig_param[6])


fig3 = plt.figure(figsize=(6,4))
ax3 = fig3.add_subplot(111)
pltall = ax3.plot(percent_all,cumulative_all/cumulative_all[-1]*100.,linewidth=4,color = '#000000',label="ALL")
pltpmrv = ax3.plot(percent_pmrv,cumulative_pmrv/cumulative_pmrv[-1]*100,linewidth=2,marker='s', color='#FF00FF',label="PM+RV")
pltdistrv = ax3.plot(percent_distrv,cumulative_distrv/cumulative_distrv[-1]*100,linewidth=2,marker='*',color='#00AAFF',label="DIST+RV")
pltpmdist = ax3.plot(percent_pmdist,cumulative_pmdist/cumulative_pmdist[-1]*100,linewidth=2,marker='h',color='#FFAA00',label="PM+DIST")
pltrv = ax3.plot(percent_rv,cumulative_rv/cumulative_rv[-1]*100,linewidth=1,marker='o',color='#0000FF',label="RV")
pltpm = ax3.plot(percent_pm,cumulative_pm/cumulative_pm[-1]*100,linewidth=1,marker='v',color='#FF0000',label="PM")
pltdist = ax3.plot(percent_dist,cumulative_dist/cumulative_dist[-1]*100,linewidth=1,marker='+',color='#00AA00',label="DIST")

ax3.set_ylim((0,100))
ax3.set_xlim((0,100))
ax3.set_xlabel('Minimum Percentage Accepted')
ax3.set_ylabel('Percent of Members Recovered')

plt.legend(loc='upper right',numpoints=1)

#plt.title('{0:}'.format(entry[0]))

plt.savefig('montecarlo/{0:}{1:}_percentages.eps'.format(groupname,youth))
plt.close()

fig4 = plt.figure(figsize=(6,4))
ax4 = fig4.add_subplot(111)
pltall = ax4.plot(percent_all,cumfalse_all/cumtotal_all*100.,linewidth=4,color = '#000000',label="ALL")
pltpmrv = ax4.plot(percent_pmrv,cumfalse_pmrv/cumtotal_pmrv*100,linewidth=2,marker='s', color='#FF00FF',label="PM+RV")
pltdistrv = ax4.plot(percent_distrv,cumfalse_distrv/cumtotal_distrv*100,linewidth=2,marker='*',color='#00AAFF',label="DIST+RV")
pltpmdist = ax4.plot(percent_pmdist,cumfalse_pmdist/cumtotal_pmdist*100,linewidth=2,marker='h',color='#FFAA00',label="PM+DIST")
pltrv = ax4.plot(percent_rv,cumfalse_rv/cumtotal_rv*100,linewidth=1,marker='o',color='#0000FF',label="RV")
pltpm = ax4.plot(percent_pm,cumfalse_pm/cumtotal_pm*100,linewidth=1,marker='v',color='#FF0000',label="PM")
pltdist = ax4.plot(percent_dist,cumfalse_dist/cumtotal_dist*100,linewidth=1,marker='+',color='#00AA00',label="DIST")

ax4.set_ylim((0,100))
ax4.set_xlim((0,100))
ax4.set_xlabel('Minimum Percentage Accepted')
ax4.set_ylabel('False positive contamination')

plt.legend(loc='upper right',numpoints=1)

#plt.title('{0:}'.format(entry[0]))
plt.savefig('montecarlo/{0:}{1:}_falsepositive.eps'.format(groupname,youth))
plt.close()

fig5 = plt.figure(figsize=(7,5))
ax5 = fig5.add_subplot(111)
pltall = ax5.plot(cumulative_all/cumulative_all[-1]*100,cumfalse_all/cumtotal_all*100.,linewidth=4,color = '#000000',label="ALL")
pltpmrv = ax5.plot(cumulative_pmrv/cumulative_pmrv[-1]*100,cumfalse_pmrv/cumtotal_pmrv*100,linewidth=2,marker='s', color='#FF00FF',label="PM+RV")
pltdistrv = ax5.plot(cumulative_distrv/cumulative_distrv[-1]*100,cumfalse_distrv/cumtotal_distrv*100,linewidth=2,marker='*',color='#00AAFF',label="DIST+RV")
pltpmdist = ax5.plot(cumulative_pmdist/cumulative_pmdist[-1]*100,cumfalse_pmdist/cumtotal_pmdist*100,linewidth=2,marker='h',color='#FFAA00',label="PM+DIST")
pltrv = ax5.plot(cumulative_rv/cumulative_rv[-1]*100,cumfalse_rv/cumtotal_rv*100,linewidth=1,marker='o',color='#0000FF',label="RV")
pltpm = ax5.plot(cumulative_pm/cumulative_pm[-1]*100,cumfalse_pm/cumtotal_pm*100,linewidth=1,marker='v',color='#FF0000',label="PM")
pltdist = ax5.plot(cumulative_dist/cumulative_dist[-1]*100,cumfalse_dist/cumtotal_dist*100,linewidth=1,marker='+',color='#00AA00',label="DIST")

ax5.set_ylim((0,100))
ax5.set_xlim((0,100))
ax5.set_xlabel('Percent of Members Recovered')
ax5.set_ylabel('False positive contamination')

plt.legend(loc='lower right')

#plt.title('{0:}'.format(entry[0]))
plt.savefig('montecarlo/{0:}{1:}_falserecovery.eps'.format(groupname,youth))
plt.close()
