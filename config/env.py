import os
tmpDirectory = "tmp/"
outputDirectory = "output/"
msms_bin= ""
if 'MSMS_BIN' in os.environ:
   msms_bin = os.environ['MSMS_BIN']
else:
  print "ERROR: MSMS_BIN not set. Variable should point to MSMS program."
  sys.exit(1)
# epsilon_area_change: for an atom to be considered an interface atom, its SES must change 
#       this much between bound and unbound.
epsilon_area_change = 0.01 

