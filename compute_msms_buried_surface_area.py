#!/usr/local/bin/python
import sys
import pdb
import os
from subprocess import Popen, PIPE
from utils.mergechains import *
from utils.xyzrn import *
from config import env
from sets import Set
from Bio.PDB import * 

# Pablo Gainza 2016-2017. LPDI IBI STI EPFL
if len(sys.argv) > 3:
  print "Compute and output the interface surfaces of a set of PDB chains."
  print "Currently this program receives each chain for a PDB as a separate input."
  print "Usage: "+sys.argv[0]+" {chain1} {chain2} [chain3]... "
  sys.exit(1)
pdbID = sys.argv[1][1:4]
full_input_filename = sys.argv[1].split("/")[-1]
basename = full_input_filename[0:4]
complex_file = env.tmpDirectory+basename+"_complex.pdb"

mergePDBs(sys.argv[1:4], complex_file)

# Now go through every PDB file and generate an xyzrn file. 
allPDB_files = list(sys.argv[1:4])
allPDB_files.append(complex_file)

filename_roots = []
xyzrn_filename_template = env.tmpDirectory+"{0}.xyzrn"
for pdbfile in allPDB_files:
  chain = pdbfile.split("/")[-1][4]
  filename_roots.append(basename+chain)
  output_pdb_as_xyzrn(pdbfile, xyzrn_filename_template.format(basename+chain))

# Now run MSMS for every xyzrn file
for root in filename_roots:
  FNULL = open(os.devnull, 'w')
  full_filename = xyzrn_filename_template.format(root)
  args= [env.msms_bin, "-if",full_filename,"-of",full_filename, "-af", full_filename]
  print env.msms_bin+" "+`args`
  p2 = Popen(args, stdout=PIPE, stderr=PIPE)
  stdout, stderr = p2.communicate()
  print stdout 
  print stderr
# Now compute the buried interface 
# We compute the buried interface by reading first the area file of unbound and bound complexes. 
#  We will have a dictionary with each atom name and its Solvent excluded surface (SES)
bound_ses = {}
# First read the SES for the bound 
area_filename_template = env.tmpDirectory+"{0}.area"
# The bound complex is indexd by the last number
complex_area_filename = area_filename_template.format(basename+"_")
bound_ses_file = open(complex_area_filename)
next(bound_ses_file) # ignore header line
for line in bound_ses_file:
  fields = line.split()
  bound_ses[fields[3]] = fields[1]

# Now read each of the unbound files and output interface and triangulation files for 
# each unbound file.
for root in filename_roots[0:-1]:
  chain_ses_file = open(area_filename_template.format(root))
  next(chain_ses_file) # ignore header
  iface_atoms = Set()
  for line in chain_ses_file:
    fields = line.split()
    atomName = fields[3]
    ses_value = float(fields[1])
    # If this atom is not in the complex surface, or its SES has changed more than epsilon
    if atomName not in bound_ses.keys() or (ses_value - float(bound_ses[atomName]) > env.epsilon_area_change):
      iface_atoms.add(atomName)
  # Now we will create new output files for the vertices and triangulations
  # But first we need to read in the vertices
  vert_filename = env.tmpDirectory+"{}.xyzrn.vert".format(root)
  vert_file = open(vert_filename)
  # Ignore first three lines
  next(vert_file)
  next(vert_file)
  next(vert_file)
  # Go through the rest of the file
  old_vertices = []
  old_vertices.append("IGNORE") # ignore first element.
  vertices_in_interface_old_indexing = Set()
  ix_counter = 1
  for line in vert_file:
    fields = line.split()
    atomName = fields[9]
    if atomName in iface_atoms: 
      vertices_in_interface_old_indexing.add(ix_counter)
    old_vertices.append(fields[0:6]+[fields[9]])
    ix_counter += 1

  # Now read the lines that correspond to the face triangles.
  face_filename = env.tmpDirectory+"{}.xyzrn.face".format(root)
  face_file = open(face_filename)
  # Ignore first three lines
  next(face_file)
  next(face_file)
  next(face_file)
  tris_old_ix_dict = {}
  for line in face_file:
    fields = line.split()
    # the tree is defined by the three coordinates
    T = map(lambda field: int(field), fields[0:3])
    # Is any of the three vertices in T in the interface? 
    is_vertex_in_interface = len(filter(lambda x: x in vertices_in_interface_old_indexing, T)) > 0
    if is_vertex_in_interface: 
      for v in T:
        if v not in tris_old_ix_dict:
          tris_old_ix_dict[v] = []
      tris_old_ix_dict[v].append(T) # Add this triangle to the triangles for this vertex.
  # Now go through all vertices again and select those in the interface
  map_from_old_to_new = {}
  map_from_new_to_old = {}
  new_vertex_ix = 1
  for old_vertex_ix in range(1, len(old_vertices)):
    if old_vertex_ix in tris_old_ix_dict:
      map_from_old_to_new[old_vertex_ix] = new_vertex_ix
      map_from_new_to_old[new_vertex_ix] = old_vertex_ix
      new_vertex_ix += 1
  # We now have the all the new indices and a map between them
  # Create the new vertices
  out_file_vertices = open(env.outputDirectory+root+".vert", 'w')
  for vertex_ix in range(1,new_vertex_ix):
    new_vertex = [`vertex_ix`]+old_vertices[map_from_new_to_old[vertex_ix]]
    out_file_vertices.write(' '.join(new_vertex)+'\n')
  # Now create a new set with the new triangles
  tris_new_ix_set = Set()
  for key in tris_old_ix_dict.keys():
    old_tris = tris_old_ix_dict[key]
    for old_tri in old_tris:
      new_tri = tuple(map(lambda v: map_from_old_to_new[v], old_tri))
      tris_new_ix_set.add(new_tri)
  # Finally, write the new triangles to the output file.
  out_file_face = open(env.outputDirectory+root+".face", 'w')
  for tri in tris_new_ix_set:
    out_file_face.write(' '.join(str(x) for x in tri)+'\n')

