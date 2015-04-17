# lacewing
Repository for the LACEwING moving group identification code (Riedel et al. in prep)

LocAting Constituent mEmbers In Nearby Groups (LACEwING) is an astrophysics code that uses the kinematics (positions and motions) of stars to determine if they are members of one of the 10 nearby young moving groups or 4 nearby open clusters. It considers membership in:
epsilon Chameleon
eta Chameleon (open cluster)
TW Hydra
beta Pic
Octans
Tuc-Hor
Columba
Argus
AB Doradus
The Pleiades (open cluster)
Hercules-Lyra
Coma Berenices (open cluster)
Ursa Major
The Hyades (open cluster)

It is written for Python 2.7 and depends upon Numpy, Scipy, and Astropy modules.

The basic operation is simple: 

1. Prepare a .csv file with a header containing the name, RA, eRA, DEC, eDEC, pi, epi, rv, erv, pmra, epmra, pmdec, epmdec (the order does not matter) and fill in the columns with information about the stars of interest. Empty columns will be either ignored or errors estimated by the code.

2. Run LACEwING:
python lacewing.py <filename> [<young?> <outputfilename> <verboseoutput?>]

 * If you do not specify any arguments, you will get a help message.

 * If you input a filename and nothing else, LACEwING will use field star parameters and regular output (in .csv format) will be sent to lacewing_output.csv

 * If you specify "young" as the third keyword, LACEwING will treat all the stars as if they are unlikely to be ordinary field objects. This increases the estimated probabilities that the star is a member of one of the nearby young moving groups.

 * If you specify "verbose" as the fifth keyword, LACEwING will print out details information about how well the stars matched to each of the 14 moving groups, with information about how the individual proper motion, distance, RV, and spatial position matches turned out... which is useful for in-depth studies of why LACEwING predicted what it did. Verbose output is 14 lines per star, rather than all the information on one line.



To run the LACEwING code itself requires four helper modules and one .csv file:
lacewing.py
converge.py
kinematics.py
ellipse.py
astrometry.py
Moving_Groups_all.csv



---------------------------------------------------------------

TRACEwING is an epicyclic traceback program that's designed to determine the distance between a star and the center of a moving group as a function of time going back into the past.

The theory is that all the stars in a moving group formed in the same place at roughly the same time, so running the clock back should put a true member very close to the center of the moving group, at the time the moving group formed.

Operation:
python tracewing.py <inputfile> <group_to_fit> <method> <ending_timestep> <age_of_group> <number_of_monte_carlo_iterations>

Inputfile is a .csv file of the same format that LACEwING requires, which should contain kinematic information for all the stars you're interested in.

"Group_to_Fit" is the name of a moving group. which should have a file with a name like: Moving_group_<name>_epicyclic.dat in the same folder. That file should contain the ellipse parameters of the nearby young moving group as a function of time, and can be generated with tracewing_mgp.py

"Method" must currently be "epicyclic" to use the epicyclic traceback technique.

"Ending Timestep" must be the time in Myr ago you want the calculations to stop. (for instance, -40 for 40 Myr ago). It must be no farther back than the time you ran the moving group itself back to.

"Age of Group" is for the vertical line that plots at the age of the moving group (for instance, for beta Pic use -25).

"Number of Monte Carlo Iterations" governs the accuracy of the traceback. The more traces, the more memory it takes. Unless you have a lot of memory, set this to 10000 or so.

TRACEwING requires one module and a set of .dat files containing the ellipse parameters of the moving groups as a function of time.
tracewing.py
kinematics.py
