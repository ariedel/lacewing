import numpy
import scipy
from astropy.io import ascii
from matplotlib import pyplot as plt
from sys import argv
import astrometry
import lacewing

def starsort(comparisons,sortquantity,sortkey,i,span):
    if astrometry.isnumber(comparisons[i][sortkey]):
        sig = sortquantity[0]
        loc = comparisons[i:i+span]['Group'][sortquantity] == comparisons[i]['Note'].split(";")[0]
        if numpy.sum(loc) == 0:
            loc = -1
        else:
            loc = numpy.argmax(loc)
        val = comparisons[i:i+span][sortquantity[0]][sortkey]
    else:
        sig = -1
        loc = -1
        val = -1
    return sig,loc,val

filename = argv[1]
file = open(filename,'rb')
readtable = ascii.get_reader(Reader=ascii.Basic)
readtable.header.splitter.delimiter = ','
readtable.data.splitter.delimiter = ','
readtable.header.start_line = 0
readtable.data.start_line = 1
comparisons = readtable.read(file)
file.close()

foldername = argv[2]

#print comparisons['sig_dist']

moving_groups = lacewing.moving_group_loader()

moving_groups = [moving_groups[x] for x in range(len(moving_groups)) if moving_groups[x].name != "Field"]

span = len(moving_groups)

groups = []
for i in range(len(moving_groups)):
    groups.append(moving_groups[i].name)

mgp_real = []
mgp_note = []
val_final = []
sig_final = []
sig_pm = []
sig_dist = []
sig_rv = []
sig_pos = []
loc_final = []
loc_pm = []
loc_dist = []
loc_rv = []
loc_pos = []

print len(comparisons)/span,len(comparisons)
j=0
for i in numpy.arange(0,len(comparisons),span):

    sort_final = numpy.argsort(numpy.ma.asarray(comparisons[i:i+span]['Probability']))[::-1]
    sort_pm =  numpy.argsort(numpy.ma.asarray(comparisons[i:i+span]['sig_pm']))
    sort_dist = numpy.argsort(numpy.ma.asarray(comparisons[i:i+span]['sig_dist']))
    sort_rv =  numpy.argsort(numpy.ma.asarray(comparisons[i:i+span]['sig_rv']))
    sort_pos = numpy.argsort(numpy.ma.asarray(comparisons[i:i+span]['sig_pos']))

    mgp_real.append(comparisons[i]['Note'].split(";")[0])
    mgp_note.append(comparisons[i]['Note'].split(";")[1])


    sig,loc,val = starsort(comparisons,sort_final,'Probability',i,span)
    sig_final.append(sig)
    loc_final.append(loc)
    val_final.append(val)
    sig,loc,val = starsort(comparisons,sort_pm,'sig_pm',i,span)
    sig_pm.append(sig)
    loc_pm.append(loc)
    sig,loc,val = starsort(comparisons,sort_dist,'sig_dist',i,span)
    sig_dist.append(sig)
    loc_dist.append(loc)
    sig,loc,val = starsort(comparisons,sort_rv,'sig_rv',i,span)
    sig_rv.append(sig)
    loc_rv.append(loc)
    sig,loc,val = starsort(comparisons,sort_pos,'sig_pos',i,span)
    sig_pos.append(sig)
    loc_pos.append(loc)



print loc_dist
mgp_real = numpy.asarray(mgp_real)
mgp_note = numpy.asarray(mgp_note)
sig_final = numpy.asarray(sig_final)
sig_pm = numpy.asarray(sig_pm)
sig_dist = numpy.asarray(sig_dist)
sig_rv = numpy.asarray(sig_rv)
sig_pos = numpy.asarray(sig_pos)
loc_final = numpy.asarray(loc_final)
loc_pm = numpy.asarray(loc_pm)
loc_dist = numpy.asarray(loc_dist)
loc_rv = numpy.asarray(loc_rv)
loc_pos = numpy.asarray(loc_pos)
val_final = numpy.asarray(val_final)
print val_final



outfile = open('{0:}/output.forward'.format(foldername),'wb')

for i in range(len(groups)):
    print groups[i]
    num_real = len(numpy.where((numpy.asarray(mgp_real) == groups[i]))[0])

    matched_final = len(numpy.where((mgp_real == groups[i]) & (sig_final == i))[0])
    good_final = len(numpy.where((mgp_real == groups[i]) & (sig_final == i) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
    print numpy.where((mgp_real ==groups[i]))[0]#,mgp_real,groups[i]
    forplot_final = loc_final[numpy.where((mgp_real == groups[i]))[0]]
    num_final = len(numpy.where((mgp_real == groups[i]) & (sig_final != -1))[0]) 
    num_good_final = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_final != -1))[0]) 

    matched_pm = len(numpy.where((mgp_real == groups[i]) & (sig_pm == i))[0])
    good_pm = len(numpy.where((mgp_real == groups[i]) & (sig_pm == i) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
    forplot_pm = loc_pm[numpy.where((mgp_real == groups[i]))]
    num_pm = len(numpy.where((mgp_real == groups[i]) & (sig_pm != -1))[0]) 
    num_good_pm = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_pm != -1))[0]) 

    matched_dist = len(numpy.where((mgp_real ==groups[i]) & (sig_dist == i))[0])
    good_dist = len(numpy.where((mgp_real ==groups[i]) & (sig_dist == i) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
    forplot_dist = loc_dist[numpy.where((mgp_real ==groups[i]))]
    num_dist = len(numpy.where((mgp_real ==groups[i]) & (sig_dist != -1))[0]) 
    num_good_dist = len(numpy.where((mgp_real ==groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_dist != -1))[0]) 

    matched_rv = len(numpy.where((mgp_real ==groups[i]) & (sig_rv == i))[0])
    good_rv = len(numpy.where((mgp_real ==groups[i]) & (sig_rv == i) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
    forplot_rv = loc_rv[numpy.where((mgp_real ==groups[i]))]
    num_rv = len(numpy.where((mgp_real ==groups[i]) & (sig_rv != -1))[0]) 
    num_good_rv = len(numpy.where((mgp_real ==groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_rv != -1))[0]) 

    matched_pos = len(numpy.where((mgp_real ==groups[i]) & (sig_pos == i))[0])
    good_pos = len(numpy.where((mgp_real ==groups[i]) & (sig_pos == i) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
    forplot_pos = loc_pos[numpy.where((mgp_real ==groups[i]))]
    num_pos = len(numpy.where((mgp_real ==groups[i]) & (sig_pos != -1))[0]) 
    num_good_pos = len(numpy.where((mgp_real ==groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_pos != -1))[0]) 


    fig = plt.figure(figsize=(10,18))
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(322)
    ax3 = fig.add_subplot(323)
    ax4 = fig.add_subplot(324)
    ax5 = fig.add_subplot(325)

#    print forplot_final
    ax1.plot(forplot_final,'b.')
    ax2.plot(forplot_pm,'b.')
    ax3.plot(forplot_dist,'b.')
    ax4.plot(forplot_rv,'b.')
    ax5.plot(forplot_pos,'b.')

    ax1.set_ylim((-1,span))
    ax2.set_ylim((-1,span))
    ax3.set_ylim((-1,span))
    ax4.set_ylim((-1,span))
    ax5.set_ylim((-1,span))

    ax1.set_title('FINAL: {0:3d} {1:3d} / {2:3d} {3:3d}'.format(matched_final,good_final,num_final,num_good_final))
    ax2.set_title('PM: {0:3d} {1:3d} / {2:3d} {3:3d}'.format(matched_pm,good_pm,num_pm,num_good_pm))
    ax3.set_title('DIST: {0:3d} {1:3d} / {2:3d} {3:3d}'.format(matched_dist,good_dist,num_dist,num_good_dist))
    ax4.set_title('RV: {0:3d} {1:3d} / {2:3d} {3:3d}'.format(matched_rv,good_rv,num_rv,num_good_rv))
    ax5.set_title('POS: {0:3d} {1:3d} / {2:3d} {3:3d}'.format(matched_pos,good_pos,num_pos,num_good_pos))

    outfile.write('\n{0:},{1:},{2:}\n'.format(i,groups[i],num_real))
    outfile.write('FINAL:,{0:3d},{1:3d},/,{2:3d},{3:3d}\n'.format(matched_final,good_final,num_final,num_good_final))
    outfile.write('PM:,{0:3d},{1:3d},/,{2:3d},{3:3d}\n'.format(matched_pm,good_pm,num_pm,num_good_pm))
    outfile.write('DIST:,{0:3d},{1:3d},/,{2:3d},{3:3d}\n'.format(matched_dist,good_dist,num_dist,num_good_dist))
    outfile.write('RV:,{0:3d},{1:3d},/,{2:3d},{3:3d}\n'.format(matched_rv,good_rv,num_rv,num_good_rv))
    outfile.write('POS:,{0:3d},{1:3d},/,{2:3d},{3:3d}\n'.format(matched_pos,good_pos,num_pos,num_good_pos))

    plt.savefig('{0:}/forward_{1:}.png'.format(foldername,groups[i].replace(' ','_')))
    plt.clf()
    plt.close()

outfile.close()


outfile = open('{0:}/output.reverse'.format(foldername),'wb')

outfile.write('Real,  Matched,  Matched_Members, Good_Members,  Above_Threshold, Threshold_Members,  Threshold_Good_Members\n')

for i in range(len(groups)):
    num_real = len(numpy.where((numpy.asarray(mgp_real) == groups[i]))[0])
    num_good = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
    outfile.write('\n{0:},{1:},{2:},\n'.format(groups[i],num_real,num_good))

    matched_final = len(numpy.where((sig_final == i))[0])
    matched_real_final = len(numpy.where((mgp_real == groups[i]) & (sig_final == i))[0])
    matched_good_final = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_final == i))[0])
    matched_threshold_final = len(numpy.where((sig_final == i) & (val_final > 20.0))[0])
    matched_threshold_real_final = len(numpy.where((mgp_real == groups[i]) & (sig_final == i) & (val_final > 20.0))[0])
    matched_threshold_good_final = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_final == i) & (val_final > 20.0))[0])

    matched_pm = len(numpy.where((sig_pm == i))[0])
    matched_real_pm = len(numpy.where((mgp_real == groups[i]) & (sig_pm == i))[0])
    matched_good_pm = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_pm == i))[0])
    #matched_threshold_pm = len(numpy.where((sig_pm == i) & (sig_pm > 20.0))[0])
    #matched_threshold_real_pm = len(numpy.where((mgp_real == groups[i]) & (sig_pm == i) & (sig_pm > 20.0))[0])
    #matched_threshold_good_pm = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_pm == i) & (sig_pm > 20.0))[0])

    matched_dist = len(numpy.where((sig_dist == i))[0])
    matched_real_dist = len(numpy.where((mgp_real == groups[i]) & (sig_dist == i))[0])
    matched_good_dist = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_dist == i))[0])
    #matched_threshold_dist = len(numpy.where((sig_dist == i) & (sig_dist > 20.0))[0])
    #matched_threshold_real_dist = len(numpy.where((mgp_real == groups[i]) & (sig_dist == i) & (sig_dist > 20.0))[0])
    #matched_threshold_good_dist = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_dist == i) & (sig_dist > 20.0))[0])

    matched_rv = len(numpy.where((sig_rv == i))[0])
    matched_real_rv = len(numpy.where((mgp_real == groups[i]) & (sig_rv == i))[0])
    matched_good_rv = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_rv == i))[0])
    #matched_threshold_rv = len(numpy.where((sig_rv == i) & (sig_rv > 20.0))[0])
    #matched_threshold_real_rv = len(numpy.where((mgp_real == groups[i]) & (sig_rv == i) & (sig_rv > 20.0))[0])
    #matched_threshold_good_rv = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_rv == i) & (sig_rv > 20.0))[0])

    matched_pos = len(numpy.where((sig_pos == i))[0])
    matched_real_pos = len(numpy.where((mgp_real == groups[i]) & (sig_pos == i))[0])
    matched_good_pos = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_pos == i))[0])
    #matched_threshold_pos = len(numpy.where((sig_pos == i) & (sig_pos > 20.0))[0])
    #matched_threshold_real_pos = len(numpy.where((mgp_real == groups[i]) & (sig_pos == i) & (sig_pos > 20.0))[0])
    #matched_threshold_good_pos = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_pos == i) & (sig_pos > 20.0))[0])




    outfile.write('{0:},{1:},{2:},{3:},{4:},{5:},\n'.format(matched_final,matched_real_final,matched_good_final,matched_threshold_final,matched_threshold_real_final,matched_threshold_good_final))
    outfile.write('{0:},{1:},{2:},\n'.format(matched_pm,matched_real_pm,matched_good_pm))
    outfile.write('{0:},{1:},{2:},\n'.format(matched_dist,matched_real_dist,matched_good_dist))
    outfile.write('{0:},{1:},{2:},\n'.format(matched_rv,matched_real_rv,matched_good_rv))
    outfile.write('{0:},{1:},{2:},\n'.format(matched_pos,matched_real_pos,matched_good_pos))

outfile.close()

# now, contamination levels:


for i in range(len(groups)):
    print groups[i]
    outfile = open('{0:}/output.contamination.{1:}'.format(foldername,groups[i].replace(' ','_')),'wb')
    num_real = len(numpy.where((numpy.asarray(mgp_real) == groups[i]))[0])
    for j in range(len(groups)):

        matched_final = len(numpy.where((mgp_real == groups[i]) & (sig_final == j))[0])
        good_final = len(numpy.where((mgp_real == groups[i]) & (sig_final == j) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
        #    print numpy.where((mgp_real ==groups[i]))[0],mgp_real,groups[i]
        forplot_final = loc_final[numpy.where((mgp_real ==groups[i]))[0]]
        num_final = len(numpy.where((mgp_real == groups[i]) & (sig_final != -1))[0]) 
        num_good_final = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_final != -1))[0]) 

        matched_pm = len(numpy.where((mgp_real == groups[i]) & (sig_pm == j))[0])
        good_pm = len(numpy.where((mgp_real == groups[i]) & (sig_pm == j) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
        forplot_pm = loc_pm[numpy.where((mgp_real == groups[i]))]
        num_pm = len(numpy.where((mgp_real == groups[i]) & (sig_pm != -1))[0]) 
        num_good_pm = len(numpy.where((mgp_real == groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_pm != -1))[0]) 
        
        matched_dist = len(numpy.where((mgp_real ==groups[i]) & (sig_dist == j))[0])
        good_dist = len(numpy.where((mgp_real ==groups[i]) & (sig_dist == j) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
        forplot_dist = loc_dist[numpy.where((mgp_real ==groups[i]))]
        num_dist = len(numpy.where((mgp_real ==groups[i]) & (sig_dist != -1))[0]) 
        num_good_dist = len(numpy.where((mgp_real ==groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_dist != -1))[0]) 
        
        matched_rv = len(numpy.where((mgp_real ==groups[i]) & (sig_rv == j))[0])
        good_rv = len(numpy.where((mgp_real ==groups[i]) & (sig_rv == j) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
        forplot_rv = loc_rv[numpy.where((mgp_real ==groups[i]))]
        num_rv = len(numpy.where((mgp_real ==groups[i]) & (sig_rv != -1))[0]) 
        num_good_rv = len(numpy.where((mgp_real ==groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_rv != -1))[0]) 
        
        matched_pos = len(numpy.where((mgp_real ==groups[i]) & (sig_pos == j))[0])
        good_pos = len(numpy.where((mgp_real ==groups[i]) & (sig_pos == j) & ((mgp_note == "Good") | (mgp_note == "(good)")))[0])
        forplot_pos = loc_pos[numpy.where((mgp_real ==groups[i]))]
        num_pos = len(numpy.where((mgp_real ==groups[i]) & (sig_pos != -1))[0]) 
        num_good_pos = len(numpy.where((mgp_real ==groups[i]) & ((mgp_note == "Good") | (mgp_note == "(good)")) & (sig_pos != -1))[0]) 

        outfile.write('\n{0:},{1:},{2:}\n'.format(j,groups[j],num_real))
        outfile.write('FINAL:,{0:3d},{1:3d},/,{2:3d},{3:3d}\n'.format(matched_final,good_final,num_final,num_good_final))
        outfile.write('PM:,{0:3d},{1:3d},/,{2:3d},{3:3d}\n'.format(matched_pm,good_pm,num_pm,num_good_pm))
        outfile.write('DIST:,{0:3d},{1:3d},/,{2:3d},{3:3d}\n'.format(matched_dist,good_dist,num_dist,num_good_dist))
        outfile.write('RV:,{0:3d},{1:3d},/,{2:3d},{3:3d}\n'.format(matched_rv,good_rv,num_rv,num_good_rv))
        outfile.write('POS:,{0:3d},{1:3d},/,{2:3d},{3:3d}\n'.format(matched_pos,good_pos,num_pos,num_good_pos))

    outfile.close()
