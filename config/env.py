import os
tmpDirectory = "tmp/"
msms_bin= ""
if 'MSMS_BIN' in os.environ:
   msms_bin = os.environ['MSMS_BIN']
else:
  print "ERROR: MSMS_BIN not set. Variable should point to MSMS program."
  sys.exit(1)
