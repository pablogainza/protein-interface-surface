import os
tmpDirectory = "tmp/"
outputDirectory = "output/"
outputDatabaseJsonFile = outputDirectory+"db.json"
msms_bin= ""
if 'MSMS_BIN' in os.environ:
   msms_bin = os.environ['MSMS_BIN']
else:
  print "ERROR: MSMS_BIN not set. Variable should point to MSMS program."
  sys.exit(1)
# epsilon_area_change: for an atom to be considered an interface atom, its SES must change 
#       this much between bound and unbound.
epsilon_area_change = 0.1
# minimum_ppi_area: Any interface that has a size less than minimum_ppi_area will be 
#       discarded: it is too small to be relevant.
minimum_ppi_area = 100
# Require that all three vertices of a triangle be in the interface to include said triangle. 
require_all_three_vertices = True
# Ignore hydrogens; color nitrogens blue.
ignore_hydrogens = True
