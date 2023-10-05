'''
Mark Hartigan
February 21, 2023

Description: defines generateScripts() and runScripts(), which generate and run
(respectively) GMAT scripts for a suite of training data to feed into a
surrogate model. This uses a predefined template script and modifies it using
parameters that are defined globally.

NOTE: This script only works on python 3.7, which is an issue with the GMAT API
'''

from load_gmat import *    # GMAT API, need to download

import os                       # to read/write files

# Orbital elements over which to iterate
sma = [1858, 1861, 1864]
ecc = [0.023, 0.043, 0.063]
inc = [89, 90, 91]
raan= [-1, 0, 1]
aop = [269, 270, 271]
ta  = [-1, 0, 1]

# Name of GMAT variable corresponding to each orbital element
lines = ['GMAT Orbiter.SMA', 'GMAT Orbiter.ECC', 'GMAT Orbiter.INC',
         'GMAT Orbiter.RAAN', 'GMAT Orbiter.AOP', 'GMAT Orbiter.TA']

# function to generate new file path + name for GMAT report file
path = lambda x: f"GMAT ReportFile.Filename = 'D:\Documents\Georgia Tech\_PNT\cislunar-nav\data\\training\{x}_MoonOrb.txt'\n"

def generateScripts():
    '''
    Opens the template script and iteratively modifies to create and save new
    scenarios that span the parameter space.
    '''

    # Store script in array
    temp = []
    # open template file and store all the lines in program memory
    with open('gmat/training/temp_MoonOrb.script', 'r') as template:
        for line in template:
            temp.append(line)

    idx = 0

    for i in range(len(sma)):    # iterate over orbital elements to change
        for j in range(len(ecc)):
            for k in range(len(inc)):
                for l in range(len(raan)):
                    for m in range(len(aop)):
                        for n in range(len(ta)):
                            # path + name of new script to make
                            name = f'gmat/training/{idx}_MoonOrb.script'

                            # writes new file line by line, replacing ones in template as needed
                            with open(name, 'w') as file:
                                # replace each line that starts with one of our variables with new value
                                for line in temp:
                                    if line.startswith(lines[0]):
                                        file.write(lines[0] + ' = ' + str(sma[i]) + ';\n')
                                    elif line.startswith(lines[1]):
                                        file.write(lines[1] + ' = ' + str(ecc[j]) + ';\n')
                                    elif line.startswith(lines[2]):
                                        file.write(lines[2] + ' = ' + str(inc[k]) + ';\n')
                                    elif line.startswith(lines[3]):
                                        file.write(lines[3] + ' = ' + str(raan[l]) + ';\n')
                                    elif line.startswith(lines[4]):
                                        file.write(lines[4] + ' = ' + str(aop[m]) + ';\n')
                                    elif line.startswith(lines[5]):
                                        file.write(lines[5] + ' = ' + str(ta[n]) + ';\n')
                                    elif line.startswith('GMAT ReportFile.Filename'):
                                        file.write(path(idx))
                                    else:   # if not a line we want to change, just write it to new file
                                        file.write(line)

                            idx += 1

def runScripts():
    '''
    Runs all GMAT scripts in a given directory

    NOTE: Just tries to run anything that is a file, so make sure there are
    only GMAT .script files in the directory (folders are OK)
    '''
    dir = 'gmat/training/'  # relative path to location of script files
    files = []
    for file in os.listdir(dir):    # create list of files in directory
        f = os.path.join(dir, file)

        if os.path.isfile(f):
            files.append(f)

    for f in files:                 # run each file that was found
        print('Running ' + f + '...')
        gmat.LoadScript(f)
        gmat.RunScript()

if __name__ == "__main__":          # executes iff this python script is run
    generateScripts()
    runScripts()