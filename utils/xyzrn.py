from Bio.PDB import * 
from config import env 
radii = {}
radii['N'] = "1.540000"
radii['O'] = "1.400000"
radii['C'] = "1.740000"
radii['H'] = "1.200000"
radii['S'] = "1.800000"
polarHydrogens = {}
polarHydrogens['ALA'] = ['H']
polarHydrogens['GLY'] = ['H']
polarHydrogens['SER'] = ['H', 'HG']
polarHydrogens['THR'] = ['H', 'HG1']
polarHydrogens['LEU'] = ['H']
polarHydrogens['ILE'] = ['H']
polarHydrogens['VAL'] = ['H']
polarHydrogens['ASN'] = ['H', 'HD21', 'HD22']
polarHydrogens['GLN'] = ['H', 'HE21', 'HE22']
polarHydrogens['ARG'] = ['H', 'HH11', 'HH12', 'HH21', 'HH22', 'HE']
polarHydrogens['HIS'] = ['H', 'HD1', 'HD2', 'HE1', 'HE2']
polarHydrogens['TRP'] = ['H', 'HD1', 'HD1']
polarHydrogens['PHE'] = ['H']
polarHydrogens['TYR'] = ['H', 'HH']
polarHydrogens['GLU'] = ['H', ]
polarHydrogens['ASP'] = ['H']
polarHydrogens['LYS'] = ['H', 'HZ1', 'HZ2', 'HZ3']
polarHydrogens['PRO'] = []
polarHydrogens['CYS'] = ['H']
polarHydrogens['MET'] = ['H']
def output_pdb_as_xyzrn(pdbfilename, xyzrnfilename):
  parser = PDBParser()
  struct = parser.get_structure(pdbfilename, pdbfilename)
  outfile = open(xyzrnfilename, 'w')
  for atom in struct.get_atoms():
    name = atom.get_name()
    residue = atom.get_parent()
    resname = residue.get_resname()
    atomtype = name[0]
    color = 'Green'
    if atomtype in radii and resname in polarHydrogens:
      if atomtype == 'O':
        color = 'Red'
      if atomtype == 'N':
        color = 'Blue'
      if atomtype == 'H' and not env.ignore_hydrogens :
        if name in polarHydrogens[resname]:
          color = 'Blue' # Polar hydrogens
      coords = "{:.06f} {:.06f} {:.06f}".format(atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2])
      full_id = "{}_{:d}_{}_{}_{}".format(atom.get_full_id()[2],\
                atom.get_full_id()[3][1], resname, name, color)
  
      outfile.write(coords+" "+radii[atomtype]+" 1 "+full_id+"\n")

