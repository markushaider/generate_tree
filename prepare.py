#!/usr/bin/python

import os
import sys

if(len(sys.argv)!=2):
    print 'Give the snapshot directory as argument'
    sys.exit()
if(os.access(sys.argv[1],os.R_OK)==False):
    print 'Snapshot directory not found or wrong access rights'
    sys.exit()
snapshot_directory = sys.argv[1]

# MASS_PARTICLE in solar masses
particle_mass = 1.269E9 # important for calculating mass_res_ok
boxsize = 32 # in megaparsec

# generate snapshot_names.dat
cwd = os.getcwd()
os.chdir(snapshot_directory)
os.system("ls snapshot_* > "+cwd+"/snapshot_names.dat")
os.chdir(cwd)

INBASE = "INBASE = "+ snapshot_directory + "\n"
OUTBASE = "OUTBASE = " + os.getcwd() + "/rockstar_output/\n"
SNAPSHOT_NAMES = "SNAPSHOT_NAMES = "+ os.getcwd()+"/snapshot_names.dat\n"

filename = 'rockstar-server.cfg'
file = open(filename,'r')
lines = file.readlines()
file.close()
for i in range(len(lines)):
    if(lines[i].find('INBASE')==0):
        lines[i] = INBASE
        print 'set '+ lines[i]
    if(lines[i].find('OUTBASE')==0):
        lines[i] = OUTBASE
        print 'set '+ lines[i]
    if(lines[i].find('SNAPSHOT_NAMES')==0):
        lines[i] = SNAPSHOT_NAMES
        print 'set '+ lines[i]
file = open(filename,'w')
file.writelines(lines)
file.close()

filename = 'rockstar-client.cfg'
file = open(filename,'r')
lines = file.readlines()
file.close()
for i in range(len(lines)):
    if(lines[i].find('INBASE')==0):
        lines[i] = INBASE
        #print 'set '+ lines[i]
    if(lines[i].find('OUTBASE')==0):
        lines[i] = OUTBASE
        #print 'set '+ lines[i]
    if(lines[i].find('SNAPSHOT_NAMES')==0):
        lines[i] = SNAPSHOT_NAMES
        #print 'set '+ lines[i]
file = open(filename,'w')
file.writelines(lines)
file.close()

BOX_WIDTH = 'BOX_WIDTH=' + str(boxsize)
MASS_RES_OK = 'MASS_RES_OK='+str(particle_mass*1000)+' #Halo mass above which there are probably'

filename = 'ctree.cfg'
file = open(filename,'r')
lines = file.readlines()
file.close()
for i in range(len(lines)):
    if(lines[i].find('BOX_WIDTH')==0):
        lines[i] = BOX_WIDTH
        print 'set '+ lines[i]
    if(lines[i].find('MASS_RES_OK')==0):
        lines[i] = MASS_RES_OK
        print 'set '+ lines[i]
file = open(filename,'w')
file.writelines(lines)
file.close()
