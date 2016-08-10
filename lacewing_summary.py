from astropy.io import ascii
import numpy as np
from sys import argv
import lacewing


################################################
### LACEwING Summary 1.5
### ARR 2016-08-10
### 1.5 Now uses LACEwING's own routines to handle
###     any potential number/order of groups
################################################

infile = argv[1]
lines = ascii.read(infile)

moving_groups = lacewing.moving_group_loader()

moving_groups = [moving_groups[x] for x in range(len(moving_groups)) if moving_groups[x].name != "Field"]

numgroups = len(moving_groups)

outfile = open('{0:}.summary.csv'.format(infile),'wb')

# regular output is a one line per entry summary, also in .csv form
outfile.write('Name,Note,Group,Probability,Predicted Dist,Predicted Dist uncertainty,Predicted RV,Predicted RV uncertainty,')
for i in range(len(moving_groups)):
    outfile.write('{0:},'.format(moving_groups[i].name))
outfile.write('\n')


for i in np.arange(0,len(lines),numgroups):
    #print lines[i]
    name = lines[i]["Name"]
    note = lines[i]["Note"]
    grouporder = np.argsort(lines[i:i+numgroups]['Probability'])[::-1]
    #print grouporder
    if lines[i+grouporder[0]]['Probability'] < 20:
        outfile.write('{0:},{1:},(None),,,,,,'.format(name,note))
    else:
        matchgroup = lines[i+grouporder[0]]['Group']
        matchsig   = lines[i+grouporder[0]]['Probability']
        matchdist  = lines[i+grouporder[0]]['kin_dist']
        matchedist = lines[i+grouporder[0]]['kin_edist']
        matchrv   = lines[i+grouporder[0]]['kin_rv']
        matcherv  = lines[i+grouporder[0]]['kin_erv']
        outfile.write('{0:},{1:},{2:},{3: 5.0f},{4: 6.2f},{5: 5.2f},{6:+6.2f},{7: 5.2f},'.format(name,note,matchgroup,matchsig,matchdist,matchedist,matchrv,matcherv))

    for j in np.arange(i,i+numgroups):
        outfile.write('{0: 5.0f},'.format(lines[j]['Probability']))
    outfile.write(',\n')

outfile.close()

    
