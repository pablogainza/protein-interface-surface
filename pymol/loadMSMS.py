# Pablo Gainza Cirauqui 2016 
# This pymol function loads msms files into pymol. 
from pymol import cmd
import sys, urllib, zlib
import subprocess
import os,math,re
import string
from pymol.cgo import *
import Queue
import threading
import os.path

colorDict = {'sky': [COLOR, 0.0, 0.76, 1.0 ],
             'sea': [COLOR, 0.0, 0.90, 0.5 ],
             'yellowtint': [COLOR, 0.88, 0.97, 0.02 ],
             'hotpink': [COLOR, 0.90, 0.40, 0.70 ],
             'greentint': [COLOR, 0.50, 0.90, 0.40 ],
             'blue': [COLOR, 0.0, 0.0, 1.0 ],
             'green': [COLOR, 0.0, 1.0, 0.0 ],
             'yellow': [COLOR, 1.0, 1.0, 0.0 ],
             'orange': [COLOR, 1.0, 0.5, 0.0],
             'red': [COLOR, 1.0, 0.0, 0.0],
             'black': [COLOR, 0.0, 0.0, 0.0],
             'gray': [COLOR, 0.9, 0.9, 0.9] }

def msms(fileroot, color="green", name='msms', dotSize=0.2, lineSize = 1):
  vertfilename = fileroot+".vert"
  facefilename = fileroot+".face"
  # Read vertices first
  vertfile = open(fileroot+".vert")
  vertdata = (vertfile.read().rstrip()).split('\n')
  vertfile.close()
  ## Data structures
  verts = {}
  for line in vertdata:
    fields = line.split()
    # Line is formatted: i x y z n1 n2 n3 x x x name
    xyz = fields[1:4]
    normals = fields [4:7]
    name = [fields[-1]]
    verts[fields[0]] = xyz+normals+name
  # Read faces (triangles)
  face_file_exists = False
  if os.path.isfile(fileroot+".face"):
    face_file_exists = True
  faces = []
  if face_file_exists:
    facefile = open(fileroot+".face")
    facedata = (facefile.read().rstrip()).split('\n')
    facefile.close()
    for line in facedata:
      faces.append(line.split())
  # Draw vertices and normals
  obj = []
  for vertkey in verts:
    vert = verts[vertkey]
    name = vert[-1]
    colorToAdd = colorDict[color]
    if "Blue" in name:
      colorToAdd = colorDict["blue"]
    if "Red" in name:
      colorToAdd = colorDict["red"]
    # Vertices
    obj.extend(colorToAdd)
    obj.extend([SPHERE, float(vert[0]), float(vert[1]), float(vert[2]), dotSize])
    obj.extend([LINEWIDTH, lineSize])
    # Normals 
    obj.extend([BEGIN, LINES])
    obj.extend(colorToAdd)
    obj.extend([VERTEX, float(vert[0]), float(vert[1]), float(vert[2])])
    obj.extend([VERTEX, float(vert[0])+float(vert[3]), \
              float(vert[1])+float(vert[4]), float(vert[2])+float(vert[5])])
    obj.append(END)
  cmd.load_cgo(obj,"vert_"+fileroot, 1.0)
  obj =[]
  colorToAdd = colorDict['gray']
  # Draw triangles (faces)
  if face_file_exists:
    for tri in faces: 
      pairs = [[tri[0],tri[1]], [tri[0],tri[2]], [tri[1],tri[2]]]
      for pair in pairs: 
        vert1 = verts[pair[0]]
        vert2 = verts[pair[1]]
        obj.extend([BEGIN, LINES])
        obj.extend(colorToAdd)
        obj.extend([VERTEX, float(vert1[0]), float(vert1[1]), float(vert1[2])])
        obj.extend([VERTEX, float(vert2[0]), float(vert2[1]), float(vert2[2])])
        obj.append(END)
    cmd.load_cgo(obj,"mesh_"+fileroot, 1.0)

def myfunc(resi,resn,name):
  print '%s`%s/%s' % (resn ,resi, name)
def get_obj_name(model):
  return cmd.get_names("objects", 1, "sele")

def msms_surface(selection, file_name="default.vert"):
  myspace={'myfunc':myfunc}
  full_object_name = cmd.get_names("objects", 1, "sele")[0]
  cmd.save(full_object_name, full_object_name)
  msms_bin = os.environ['MSMS_BIN']
  args= [env.msms_bin, "-if",full_object_name,"-of",full_object_name, "-af", full_object_name]
  p2 = Popen(args, stdout=PIPE, stderr=PIPE)
  stdout, stderr = p2.communicate()

  print stdout
  print stderr
  cmd.iterate(selection,'myfunc(resi,resn,name)', space=myspace)
  
  
