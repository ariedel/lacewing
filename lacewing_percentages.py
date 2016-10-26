import numpy as np
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
### ARR 2016-10-26
### 1.3: Now compatible with LACEwING v1.3, and uses the 
###      lacewing.moving_group_loader() function
### 1.4: Reads the new smaller files.
### 1.5: Fieldflip allows us to let in exactly as many 
###      field stars as young stars, to account for a 
###      young field option.
### 1.6: Reorganized and commented entire code
#########################################################
#########################################################

#########################################################
## The Program Flow
#########################################################
## 1. Load Moving Group names
## 2a. Set up histograms
## 2b. Read in cluster file and count up stars
## 3. Print out results
## 4. Compute the Gaussian Cumulative Distribution Functions that fit
##    the data, and print them to a file for copy+pasting into 
##    Moving_Groups_all.csv
## 5. First plot (<group>.eps): Percentage of total stars with that 
##    score that are "really" members, as a function of the 
##    goodness-of-fit value
## 6. Second Plot (<group>_cumulative.eps): The cumulative 
##    percentage of "real" moving group members recovered, as a 
##    function of maximum goodness-of-fit value accepted.
## 7. Third plot (<group>_percentages.eps) Number of stars recovered 
##    as a function of minimum membership percentage accepted. (basically
##    the same as Plot 2 <group>_cumulative.eps but with percentages)
## 8. Fourth plot (<group>_falsepositive.eps): The false positive rate
##    for a sample with a given membership percentage lower limit
## 9. Fifth Plot (<group>_falserecovery.eps): The false positive 
##    rate as a function of the recovery rate of stars. Combination of 
##    Plot 3 and 4, showing the amount of false positives necessary to 
##    recover a percentage of real members.
#########################################################

def gaus(x,a,sigma):
    return a*np.exp(-(x)**2/(2*sigma**2))

def gauscdf(x,m,sigma):
    return 0.5*(erfc((x-m)/(sigma*np.sqrt(2))))


#########################################################
## 1. Load Moving Group names
#########################################################

moving_groups = lacewing.moving_group_loader()

groups = []
for i in xrange(len(moving_groups)):
    groups.append(moving_groups[i].name)
groups = np.asarray(groups)

youth=""
if len(argv) > 2:
    if argv[2] == 'young':
        youth=".youngonly"

if not os.path.exists('montecarlo'):
    os.mkdir('montecarlo')

#########################################################
## 2a. Set up histograms
#########################################################

stop = 24
step = 0.1

sigrange = np.arange(0.0,stop,step)
acc_all = np.zeros_like(sigrange)
acc_pm = np.zeros_like(sigrange)
acc_dist = np.zeros_like(sigrange)
acc_rv = np.zeros_like(sigrange)
acc_pmdist = np.zeros_like(sigrange)
acc_pmrv = np.zeros_like(sigrange)
acc_distrv = np.zeros_like(sigrange)
total_all = np.zeros_like(sigrange)
total_pm = np.zeros_like(sigrange)
total_dist = np.zeros_like(sigrange)
total_rv = np.zeros_like(sigrange)
total_pmdist = np.zeros_like(sigrange)
total_pmrv = np.zeros_like(sigrange)
total_distrv = np.zeros_like(sigrange)
good_all = np.zeros((len(groups),len(sigrange)))
good_pm = np.zeros((len(groups),len(sigrange)))
good_dist = np.zeros((len(groups),len(sigrange)))
good_rv = np.zeros((len(groups),len(sigrange)))
good_pmdist = np.zeros((len(groups),len(sigrange)))
good_pmrv = np.zeros((len(groups),len(sigrange)))
good_distrv = np.zeros((len(groups),len(sigrange)))


#########################################################
## 2b. Read in cluster file and count up stars
#########################################################

i = 0
ystars = 0
# The input file contains every star as matched to one of the moving groups.
# This loop reads in every line in the file, one by one, and builds histograms
# Each line contains the goodness-of-fit values for the 7 different data cases, relative to one moving group.
# For each of the 7 data cases, there are arrays: One for total number of points
#  with that particular binned goodness-of-fit score, and N to separate out objects that were "really" 
#  members of each of the moving groups, plus one for objects that were "really" field stars.
with open(argv[1], 'rb') as f:
    for line in f:
        if i % 100000 == 0:
            print i
        entry = line.split(",")

        if argv[2] == 'young':
            if entry[8] == "Field":
                if ystars >= 1: # ystars will increment each time a moving group member is found, and decrement for every field star. Thus, we'll arrive at a 1:1 field:young ratio.
                    ystars -= 1
                else:
                    continue
            else:
                ystars = ystars +1

        grp = np.int(np.where(entry[8] == groups)[0][0])

        location = np.int(float(entry[1])/step)
        if location < stop/step:
            total_all[location] += 1
            good_all[grp][location] += 1
                    
        location = np.int(float(entry[2])/step)
        if location < stop/step:
            total_pm[location] += 1
            good_pm[grp][location] += 1

        location = np.int(float(entry[3])/step)
        if location < stop/step:
            total_dist[location] += 1
            good_dist[grp][location] += 1

        location = np.int(float(entry[4])/step)
        if location < stop/step:
            total_rv[location] += 1
            good_rv[grp][location] += 1

        location = np.int(float(entry[5])/step)
        if location < stop/step:
            total_pmdist[location] += 1
            good_pmdist[grp][location] += 1

        location = np.int(float(entry[6])/step)
        if location < stop/step:
            total_pmrv[location] += 1
            good_pmrv[grp][location] += 1

        location = np.int(float(entry[7])/step)
        if location < stop/step:
            total_distrv[location] += 1
            good_distrv[grp][location] += 1

        i+= 1

#########################################################
## 3. Print out results
#########################################################


# Now identify the moving group we're interested in...
print entry[0]
real = np.where(entry[0] == groups)[0][0]
print real

groupname = entry[0].replace(' ','_')

outfile = open('montecarlo/{0:}{1:}.values'.format(groupname,youth),'wb')

outfile.write('sigval,tot_all,good_all,tot_pmdist,good_pmdist,tot_pmrv,good_pmrv,tot_distrv,good_distrv,tot_pm,good_pm,tot_dist,good_dist,tot_rv,good_rv\n')

for m in range(len(sigrange)):

    print sigrange[m], total_all[m], good_all[real][m], total_pmdist[m], good_pmdist[real][m], total_pmrv[m], good_pmrv[real][m], total_distrv[m],good_distrv[real][m],total_pm[m],good_pm[real][m],total_dist[m],good_dist[real][m],total_rv[m],good_rv[real][m]
    outfile.write('{0:},{1:},{2:},{3:},{4:},{5:},{6:},{7:},{8:},{9:},{10:},{11:},{12:},{13:},{14:}\n'.format(sigrange[m], total_all[m], good_all[real][m], total_pmdist[m], good_pmdist[real][m], total_pmrv[m], good_pmrv[real][m], total_distrv[m],good_distrv[real][m],total_pm[m],good_pm[real][m],total_dist[m],good_dist[real][m],total_rv[m],good_rv[real][m]))

outfile.close()

#########################################################
## 4. Compute the Gaussian Cumulative Distribution Functions that fit
##    the data, and print them to a file for copy+pasting into 
##    Moving_Groups_all.csv
#########################################################

# Create arrays of the fraction of real members, as a function of goodness-of-fit
fraction_all = np.asarray(good_all[real]/total_all)
fraction_pmdist = np.asarray(good_pmdist[real]/total_pmdist)
fraction_pmrv = np.asarray(good_pmrv[real]/total_pmrv)
fraction_distrv = np.asarray(good_distrv[real]/total_distrv)
fraction_pm = np.asarray(good_pm[real]/total_pm)
fraction_dist = np.asarray(good_dist[real]/total_dist)
fraction_rv = np.asarray(good_rv[real]/total_rv)


#  Open the file
outfile = open('montecarlo/{0:}{1:}.percentages'.format(groupname,youth),'wb')

fractionlist = np.asarray([fraction_pm,fraction_dist,fraction_rv,fraction_pmdist,fraction_pmrv,fraction_distrv,fraction_all])

mean_param = []
sig_param = []

# Shift this so the curves fit to the middle of the bins rather than a corner
sigrange = sigrange+0.05

for a in range(len(fractionlist)):
    # replace NaN with 1 (educated guess here, for situations with too few objects in the first bins)
    sift1 = np.where(np.asarray(np.isnan(fractionlist[a])))
    #print sift1

    if len(sift1) > 0:
        fractionlist[a][sift1] = 1

    # Don't use the first bin at all. This tends to have behavior different from all the others, and small number statistics.
    sift9 = np.where(sigrange > 0.1)[0]

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



#########################################################
## 5. First plot (<group>.eps): Percentage of total stars with that 
##    score that are "really" members, as a function of the 
##    goodness-of-fit value
#########################################################

fig1 = plt.figure(figsize=(7,5))
ax = fig1.add_subplot(111)

others = range(len(groups))
others.remove(real)
print others

pltall = ax.plot(sigrange,fraction_all*100,linewidth=4,color='#000000',label="ALL",drawstyle='steps-mid')
pltpmrv = ax.plot(sigrange,fraction_pmrv*100,linewidth=2,marker='s',color='#FF00FF',label="PM+RV",drawstyle='steps-mid')
pltdistrv = ax.plot(sigrange,fraction_distrv*100,linewidth=2,marker='*',color='#00AAFF',label="DIST+RV",drawstyle='steps-mid')
pltpmdist = ax.plot(sigrange,fraction_pmdist*100,linewidth=2,marker='h',color='#FFAA00',label="PM+DIST",drawstyle='steps-mid')
pltrv = ax.plot(sigrange,fraction_rv*100,linewidth=1,marker='o',color='#0000FF',label="RV",drawstyle='steps-mid')
pltpm = ax.plot(sigrange,fraction_pm*100,linewidth=1,marker='v',color='#FF0000',label="PM",drawstyle='steps-mid')
pltdist = ax.plot(sigrange,fraction_dist*100,linewidth=1,marker='+',color='#00AA00',label="DIST",drawstyle='steps-mid')

# Plot the Gaussian CDF fits to the histograms
fitcolors = ['r:','g:','b:','y:','m:','c:','k:']
for a in range(len(fractionlist)):
    ax.plot(sigrange,gauscdf(sigrange,mean_param[a],sig_param[a])*100,fitcolors[a])

plt.legend(loc='upper right',numpoints=1)

ax.set_xlim(0,3.0)
ax.set_ylim(0,100)

ax.set_xlabel('Goodness of Fit')
ax.set_ylabel('Percentage of "Real" Members')
plt.tight_layout()

plt.savefig('montecarlo/{0:}{1:}.eps'.format(groupname,youth))
plt.clf()
plt.close()

#########################################################
## 6. Second Plot (<group>_cumulative.eps): The cumulative 
##    percentage of "real" moving group members recovered, as a 
##    function of maximum goodness-of-fit value accepted.
#########################################################

# generate cumulative values of members, nonmembers, and total. 
cumulative_all = np.cumsum(good_all[real])
cumulative_pmdist = np.cumsum(good_pmdist[real])
cumulative_pmrv = np.cumsum(good_pmrv[real])
cumulative_distrv = np.cumsum(good_distrv[real])
cumulative_pm = np.cumsum(good_pm[real])
cumulative_dist = np.cumsum(good_dist[real])
cumulative_rv = np.cumsum(good_rv[real])

#for x in range(len(cumfalse_all)):
#    print sigrange[x],cumulative_all[x],cumfalse_all[x],cumtotal_all[x]

fig2 = plt.figure(figsize=(7,5))
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
ax2.set_ylabel('Percentage of "Real" Members Recovered')
plt.legend(loc='lower right',numpoints=1)

plt.tight_layout()
plt.savefig('montecarlo/{0:}_cumulative.eps'.format(groupname,youth))
plt.clf()
outfile.close()

#########################################################
## 7. Third plot (<group>_percentages.eps) Number of stars recovered 
##    as a function of minimum membership percentage accepted. (basically
##    the same as Plot 2 <group>_cumulative.eps but with percentages)
#########################################################

# Generate the percentage values that map to the goodness-of-fit values
percent_pm = 100*gauscdf(sigrange,mean_param[0],sig_param[0])
percent_dist = 100*gauscdf(sigrange,mean_param[1],sig_param[1])
percent_rv = 100*gauscdf(sigrange,mean_param[2],sig_param[2])
percent_pmdist = 100*gauscdf(sigrange,mean_param[3],sig_param[3])
percent_pmrv = 100*gauscdf(sigrange,mean_param[4],sig_param[4])
percent_distrv = 100*gauscdf(sigrange,mean_param[5],sig_param[5])
percent_all = 100*gauscdf(sigrange,mean_param[6],sig_param[6])


fig3 = plt.figure(figsize=(7,5))
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
ax3.set_xlabel('Minimum Membership Percentage Accepted')
ax3.set_ylabel('Percentage of "Real" Members Recovered')

plt.legend(loc='upper right',numpoints=1)

#plt.title('{0:}'.format(entry[0]))
plt.tight_layout()
plt.savefig('montecarlo/{0:}{1:}_percentages.eps'.format(groupname,youth))
plt.close()

#########################################################
## 8. Fourth plot (<group>_falsepositive.eps): The false positive rate
##    for a sample with a given membership percentage lower limit
#########################################################

cumfalse_all = np.cumsum(total_all-good_all[real])
cumfalse_pmdist = np.cumsum(total_pmdist-good_pmdist[real])
cumfalse_pmrv = np.cumsum(total_pmrv-good_pmrv[real])
cumfalse_distrv = np.cumsum(total_distrv-good_distrv[real])
cumfalse_pm = np.cumsum(total_pm-good_pm[real])
cumfalse_dist = np.cumsum(total_dist-good_dist[real])
cumfalse_rv = np.cumsum(total_rv-good_rv[real])

cumtotal_all = np.cumsum(total_all)
cumtotal_pmdist = np.cumsum(total_pmdist)
cumtotal_pmrv = np.cumsum(total_pmrv)
cumtotal_distrv = np.cumsum(total_distrv)
cumtotal_pm = np.cumsum(total_pm)
cumtotal_dist = np.cumsum(total_dist)
cumtotal_rv = np.cumsum(total_rv)


fig4 = plt.figure(figsize=(7,5))
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
ax4.set_xlabel('Minimum Membership Percentage Accepted')
ax4.set_ylabel('False positive contamination (%)')

plt.legend(loc='upper right',numpoints=1)
plt.tight_layout()
#plt.title('{0:}'.format(entry[0]))
plt.savefig('montecarlo/{0:}{1:}_falsepositive.eps'.format(groupname,youth))
plt.close()

#########################################################
## 9. Fifth Plot (<group>_falserecovery.eps): The false positive 
##    rate as a function of the recovery rate of stars. Combination of 
##    Plot 3 and 4, showing the amount of false positives necessary to 
##    recover a percentage of real members.
#########################################################

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
ax5.set_xlabel('Percent of "Real" Members Recovered')
ax5.set_ylabel('False positive contamination (%)')

plt.legend(loc='lower right')
plt.tight_layout()
#plt.title('{0:}'.format(entry[0]))
plt.savefig('montecarlo/{0:}{1:}_falserecovery.eps'.format(groupname,youth))
plt.close()
