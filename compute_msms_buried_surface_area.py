#!/usr/local/bin/python
import sys
import pdb
import os
from utils.mergechains import *
from Bio.PDB import * 

# Pablo Gainza 2016-2017. LPDI IBI STI EPFL
if len(sys.argv) > 3:
  print "Compute and output the interface surfaces of a set of PDB chains."
  print "Currently this program receives each chain for a PDB as a separate input."
  print "Usage: "+sys.argv[0]+" {chain1} {chain2} [chain3]... "
  sys.exit(1)
msms_bin = ""
if 'MSMS_BIN' in os.environ:
   msms_bin = os.environ['MSMS_BIN']
else:
  print "ERROR: MSMS_BIN not set. Variable should point to MSMS program."
  sys.exit(1)
tmpDirectory = "/tmp/"
pdbID = sys.argv[1][1:4]
complex_file = tmpDirectory+"complex.pdb"

mergePDBs(sys.argv[1:4], complex_file)

