#!/usr/local/bin/python
# Pablo Gainza 2016-2017. LPDI IBI STI EPFL
import sys
import pdb
import os
from subprocess import Popen, PIPE
from utils.mergechains import *
from utils.xyzrn import *
from config import env
from sets import Set
from Bio.PDB import * 
import numpy as np
import shutil
from glob import glob
from scipy.spatial import distance

# Recursive function to find all connected components in a surface (to separate interfaces). 
#   In order to keep the stack small, the tree works breadth-first
def get_connected_components(cur_vertex, disconnected_components, graph):
  my_connected_components = set()
  my_connected_components.add(cur_vertex)
  neighbors = graph[cur_vertex]
  # Breadth first search tree
  neighbors_to_visit = set()
  for n in neighbors:
    my_connected_components.add(n)
    if n in disconnected_components:
      disconnected_components.discard(n)
      neighbors_to_visit.add(n)
  # Once all neighbors have been added, we go down the tree
  for n in neighbors_to_visit:
    new_set = get_connected_components(n, disconnected_components, graph)
    my_connected_components = my_connected_components.union(new_set)
  return my_connected_components

## Area of a triangle from its sides
def tri_area(a, b, c):
  s = (a + b + c) / 2
  area = (s*(s-a)*(s-b)*(s-c)) ** 0.5
  return area

def get_area_triangle_set(triangles, vertices):
  total_area = 0
  for tri in triangles:
    p0 = vertices[tri[0]]
    p1 = vertices[tri[1]]
    p2 = vertices[tri[2]]
    a = distance.euclidean(p0, p1)
    b = distance.euclidean(p1, p2)
    c = distance.euclidean(p2, p0)
    total_area += tri_area(a,b,c)
  return total_area

if "-includeH" in sys.argv:
  sys.argv.remove("-includeH")
  env.ignore_hydrogens = False
allPDB_files = list(sys.argv[1:])
if len(allPDB_files) < 2:
  print "Compute and output the interface surfaces of a set of PDB chains."
  print "Currently this program receives each chain for a PDB as a separate input."
  print "It is assumed that all chains belong to the same PDB file with the same PDB id."
  print "[-includeH]: hydrogens are included and will be colored depending on whether they are polar (blue) or not(green).\
            If this option is not set then nitrogens are painted blue and hydrogens are ignored."
  print "Usage: "+sys.argv[0]+" [-includeH] {chain1} {chain2} [chain3]... "
  sys.exit(1)

# The name of the PDB is taken only from the first file. 
full_input_filename = allPDB_files[0].split("/")[-1]
basename = full_input_filename[0:4]
complex_file = env.tmpDirectory+basename+"_complex.pdb"

## output directory 
hierarchOutputDirectory = env.outputDirectory+"/"+(basename[1:3]).lower()+"/"+basename+"/"
if not os.path.exists(hierarchOutputDirectory):
  os.makedirs(hierarchOutputDirectory)

# Copy all pdb files to destination directory
for pdb_file in allPDB_files:
  shutil.copy(pdb_file, hierarchOutputDirectory)

mergePDBs(allPDB_files, complex_file)

# Now go through every PDB file and generate an xyzrn file. 
allPDB_files.append(complex_file)

filename_roots = []
xyzrn_filename_template = env.tmpDirectory+"{0}.xyzrn"
for pdbfile in allPDB_files:
  chain = pdbfile.split("/")[-1][4]
  filename_roots.append(basename+chain)
  output_pdb_as_xyzrn(pdbfile, xyzrn_filename_template.format(basename+chain))

# Now run MSMS for every xyzrn file
for filename_root in filename_roots:
  FNULL = open(os.devnull, 'w')
  full_filename = xyzrn_filename_template.format(filename_root)
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
for filename_root in filename_roots[0:-1]:
  chain_ses_file = open(area_filename_template.format(filename_root))
  next(chain_ses_file) # ignore header
  iface_atoms = Set()
  for line in chain_ses_file:
    fields = line.split()
    atomName = fields[3]
    ses_value = float(fields[1])
    # If this atom is not in the complex surface, or its SES has changed more than epsilon
    if atomName not in bound_ses.keys() or abs(ses_value - float(bound_ses[atomName])) > env.epsilon_area_change:
      iface_atoms.add(atomName)
  # Now we will create new output files for the vertices and triangulations
  # But first we need to read in the vertices
  vert_filename = env.tmpDirectory+"{}.xyzrn.vert".format(filename_root)
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
    old_vertices.append(fields)
    ix_counter += 1

  # Now read the lines that correspond to the face triangles.
  face_filename = env.tmpDirectory+"{}.xyzrn.face".format(filename_root)
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
    is_vertex_in_interface = False
    if not env.require_all_three_vertices:
      # Is any of the three vertices in T in the interface? 
      is_vertex_in_interface = len(filter(lambda x: x in vertices_in_interface_old_indexing, T)) > 0
    else: 
      # Are all three vertices in T in the interface: 
      is_vertex_in_interface = len(filter(lambda x: x in vertices_in_interface_old_indexing, T)) > 2
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
  count_vertices_in_full_interface = new_vertex_ix
  # We now have the all the new indices and a map between them
  # Create the new vertices
  out_file_vertices = open(hierarchOutputDirectory+filename_root+".vert", 'w')
  for vertex_ix in range(1,new_vertex_ix):
    new_vertex = old_vertices[map_from_new_to_old[vertex_ix]]
    out_file_vertices.write(' '.join(new_vertex)+'\n')
  # Now create a new set with the new triangles
  tris_new_ix_set = Set()
  for key in tris_old_ix_dict.keys():
    old_tris = tris_old_ix_dict[key]
    for old_tri in old_tris:
      new_tri = tuple(map(lambda v: map_from_old_to_new[v], old_tri))
      tris_new_ix_set.add(new_tri)
  # Finally, write the new triangles to the output file.
  out_file_face = open(hierarchOutputDirectory+filename_root+".face", 'w')
  for tri in tris_new_ix_set:
    out_file_face.write(' '.join(str(x) for x in tri)+" 0 0\n")
    
  ####### Separate the vertices into sets of connected components
  # We use the information in the triangulation.
  # First convert the triangulation into a graph. 
  graph = {}
  for tri in tris_new_ix_set: 
    x, y, z = tri
    if x not in graph:
      graph[x] = set()
    graph[x].add(y)
    graph[x].add(z)
    if y not in graph:
      graph[y] = set()
    graph[y].add(x)
    graph[y].add(z)
    if z not in graph:
      graph[z] = set()
    graph[z].add(x)
    graph[z].add(y)

  # Find the connected components. 
  # Add all elements to disconnected components. 
  disconnected_components = set()
  for v in graph.keys():
    disconnected_components.add(v)
  forest = []
  # Go through the disconnected components until they all belong to a tree. 
  tree_count = 0
  while len(disconnected_components) > 0: 
    # Remove first element.
    rootv = disconnected_components.pop()
    new_tree = get_connected_components(rootv, disconnected_components, graph)
    # Create a vertex list for this tree. 
    vertex_coords = {}
    for vertex_ix in new_tree: 
      vertex_coords[vertex_ix] = (float(old_vertices[map_from_new_to_old[vertex_ix]][0]),
                                  float(old_vertices[map_from_new_to_old[vertex_ix]][1]),
                                  float(old_vertices[map_from_new_to_old[vertex_ix]][2]))
    # Find the triangles that correspond to this connected component.
    triangles_for_this_component = set()
    for triangle in tris_new_ix_set:
      x, y, z = triangle
      if x in new_tree or y in new_tree or z in new_tree: 
        triangles_for_this_component.add(triangle)
    area = get_area_triangle_set(triangles_for_this_component, vertex_coords)

    if area > env.minimum_ppi_area:
      print "Area of this interface (Angs^2): " + `area`
      # Write a vert file and a face file for the interface in this tree
      # Create the new vertices
      output_vertfilename = "{}{}{:02d}.vert".format(hierarchOutputDirectory,filename_root,tree_count)
      out_file_vertices = open(output_vertfilename, 'w')
      sorted_tree = sorted(new_tree)
      # The triangles of this connected component will have a new indexing, which means we need a new map. 
      map_from_new_ix_to_this_connected_component_ix = {}
      ix_this_connected_component = 1
      for ix_this_connected_component in range(1, len(sorted_tree)+1):
        map_from_new_ix_to_this_connected_component_ix[sorted_tree[ix_this_connected_component-1]] = ix_this_connected_component
      
      for vertex_ix in sorted_tree:
        new_vertex = old_vertices[map_from_new_to_old[vertex_ix]]
        out_file_vertices.write(' '.join(new_vertex)+'\n')
      # Now write out the faces.
      output_facefilename = "{}{}{:02d}.face".format(hierarchOutputDirectory,filename_root,tree_count)
      out_file_face = open(output_facefilename, 'w')
      for tri in triangles_for_this_component:
        tri_ix_for_this_component = []
        tri_ix_for_this_component.append(map_from_new_ix_to_this_connected_component_ix[tri[0]])
        tri_ix_for_this_component.append(map_from_new_ix_to_this_connected_component_ix[tri[1]])
        tri_ix_for_this_component.append(map_from_new_ix_to_this_connected_component_ix[tri[2]])
        outline = "{:d} {:d} {:d} 0 0".format(tri_ix_for_this_component[0], tri_ix_for_this_component[1], tri_ix_for_this_component[2])
        out_file_face.write(outline+'\n')
      tree_count += 1

