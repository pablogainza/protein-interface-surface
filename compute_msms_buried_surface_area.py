#!/usr/local/bin/python
import sys
import pdb
import os
from subprocess import Popen, PIPE
from utils.mergechains import *
from utils.xyzrn import *
from config import env
from Bio.PDB import * 

# Pablo Gainza 2016-2017. LPDI IBI STI EPFL
if len(sys.argv) > 3:
  print "Compute and output the interface surfaces of a set of PDB chains."
  print "Currently this program receives each chain for a PDB as a separate input."
  print "Usage: "+sys.argv[0]+" {chain1} {chain2} [chain3]... "
  sys.exit(1)
pdbID = sys.argv[1][1:4]
complex_file = env.tmpDirectory+"complex.pdb"

mergePDBs(sys.argv[1:4], complex_file)

# Now go through every PDB file and generate an xyzrn file. 
allPDB_files = list(sys.argv[1:4])
allPDB_files.append(complex_file)

count = 0
xyzrn_filenames =[]
for pdbfile in allPDB_files:
  xyzrn_filenames.append(env.tmpDirectory+`count`+".xyzrn")
  output_pdb_as_xyzrn(pdbfile, xyzrn_filenames[count])
  count += 1

# Now run MSMS for every xyzrn file
for xyzrn_filename in xyzrn_filenames:
  FNULL = open(os.devnull, 'w')
  args= [env.msms_bin, "-if",xyzrn_filename,"-of",xyzrn_filename, "-af", xyzrn_filename]
  print env.msms_bin+" "+`args`
  p2 = Popen(args, stdout=PIPE, stderr=PIPE)
  stdout, stderr = p2.communicate()
  print stdout 
  print stderr
