#!/usr/bin/python

import os
import sys

# generate snapshot_names.dat
cwd = os.getcwd()
os.chdir('./rockstar_output/')
os.system("perl extract_scales.pl out_*.list > DescScales.txt")
os.chdir(cwd)
print 'Wrote DescScales in ./rockstar_output/\n'
