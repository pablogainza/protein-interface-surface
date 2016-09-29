from Bio.PDB import * 
radii = {}
radii['N'] = "1.540000"
radii['O'] = "1.400000"
radii['C'] = "1.740000"
radii['H'] = "1.200000"
radii['S'] = "1.800000"
def output_pdb_as_xyzrn(pdbfilename, xyzrnfilename):
  parser = PDBParser()
  struct = parser.get_structure(pdbfilename, pdbfilename)
  outfile = open(xyzrnfilename, 'w')
  for atom in struct.get_atoms():
    name = atom.get_name()
    atomtype = name[0]
    coords = "{:.06f} {:.06f} {:.06f}".format(atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2])
    residue = atom.get_parent()
    full_id = "{}_{:d}_{}_{}".format(atom.get_full_id()[2],\
              atom.get_full_id()[3][1], residue.get_resname(), name)

    outfile.write(coords+" "+radii[atomtype]+" 1 "+full_id+"\n")

