import sys
import pdb
import os
from Bio.PDB import * 

# Merge a set of PDBs into a single one and save it.
def mergePDBs(input_pdb_filenames, output_pdb_filename):
  parser = PDBParser()
  # Read over each chain
  allChainFilenames = input_pdb_filenames
  allChains = []
  for filename in allChainFilenames:
    ## Each structure should have exactly one chain and we only care about model 0.
    struct = parser.get_structure(filename, filename)
    models = Selection.unfold_entities(struct, 'M')
    assert len(models) == 1
    chains = Selection.unfold_entities(struct, 'C')
    assert len(chains) == 1
    for chain in chains:
      chain.detach_parent()
      allChains.append(chain)

  pdbio = PDBIO()
  structBuild = StructureBuilder.StructureBuilder()
  structBuild.init_structure("output")
  structBuild.init_seg(" ")
  structBuild.init_model(0)
  outputStruct = structBuild.get_structure()
  for chain in allChains: 
    outputStruct[0].add(chain)
  
  io = PDBIO()
  io.set_structure(outputStruct)
  io.save(output_pdb_filename)
