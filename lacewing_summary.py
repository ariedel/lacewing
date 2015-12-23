from astropy.io import ascii
import numpy as np
from sys import argv

infile = argv[1]
file = open(infile,'rb')
readtable = ascii.get_reader(Reader=ascii.Basic)
readtable.header.splitter.delimiter = ','
readtable.data.splitter.delimiter = ','
readtable.header.start_line = 0
readtable.data.start_line = 1
lines = readtable.read(file)
file.close()

numgroups = 14

outfile = open('{0:}.summary.csv'.format(infile),'wb')

outfile.write('Name,Note,Group,Probability,Predicted Dist,Predicted Dist uncertainty,Predicted RV,Predicted RV uncertainty,eps Cha,eta Cha,TW Hya,beta Pic,Octans,Tuc-Hor,Columba,Argus,AB Dor,Pleiades,Her-Lyr,Coma Ber,Ursa Major,Hyades,\n')
#outfile.write('Name,Note,Group,Probability,Predicted Dist,Predicted Dist uncertainty,Predicted RV,Predicted RV uncertainty,TW_Hya,beta_Pic,Tuc-Hor,Columba,Carina,Argus,AB_Dor\n')
#outfile.write('Name,Note,Group,Probability,Predicted Dist,Predicted Dist uncertainty,Predicted RV,Predicted RV uncertainty,beta_Pic,AB_Dor,Tuc-Hor,TW_Hya,Cha-Near,Columba\n')

for i in np.arange(0,len(lines),numgroups):
    #print lines[i]
    name = lines[i]["Name"]
    note = lines[i]["note"]
    grouporder = np.argsort(lines[i:i+numgroups]['probability'])[::-1]
    #print grouporder
    if lines[i+grouporder[0]]['probability'] < 20:
        outfile.write('{0:},{1:},(None),    ,    ,    ,    ,    ,'.format(name,note))
    else:
        matchgroup = lines[i+grouporder[0]]['Group']
        matchsig   = lines[i+grouporder[0]]['probability']
        matchdist  = lines[i+grouporder[0]]['kin_dist']
        matchedist = lines[i+grouporder[0]]['kin_edist']
        matchrv   = lines[i+grouporder[0]]['kin_rv']
        matcherv  = lines[i+grouporder[0]]['kin_erv']
        outfile.write('{0:},{1:},{2:},{3: 5.0f},{4: 6.2f},{5: 5.2f},{6:+6.2f},{7: 5.2f},'.format(name,note,matchgroup,matchsig,matchdist,matchedist,matchrv,matcherv))

    for j in np.arange(i,i+numgroups):
        outfile.write('{0: 5.0f},'.format(lines[j]['probability']))
    outfile.write(',\n')

outfile.close()

    
