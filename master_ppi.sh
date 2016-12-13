#!/bin/bash 
# Pablo Gainza LPDI EPFL 2016
# This script depends on the  compute_msms_buried_surface_area.py script. 
# For any four letter PDB id, "PDBX", this script computes the buried surface areas of that PDB and
# outputs a set of files PDBX001, PDBX002, PDBX003 which contain the PPIs in PDBX. Currently, 
# it is configured to work on Jaume Bonet's PDB database at EPFL.
## Global config options 
# include hydrogens. 
INCLUDE_H=1
if [ $# -ne 1 ]
then
  echo "Usage: $0 [PDBid]"
  exit 1
fi
if [ $HOSTNAME == "macbook" ]
then
  mkdir -p pdb_database
  mounted=$(mount | grep cleanstr | wc -l | xargs)
  if [ $mounted == 0 ]
  then
    echo "mounting remote volume" 
    sshfs -o reconnect gainza@castor.epfl.ch:/work/upcorreia/databases/pdb/cleanstr pdb_database
  fi
elif [ $HOSTNAME == "castor" ]
then
  ln -s /work/upcorreia/databases/pdb/cleanstr pdb_database
fi
pdb=$1
directory=$(echo $pdb | cut -b2,3 | tr '[:upper:]' '[:lower:]')
mkdir -p tmp_master
if [ -e pdb_database/$directory/$pdb*00*.pdb.gz ]
then 
  cp pdb_database/$directory/$pdb*001.pdb.gz tmp_master
else
  cp pdb_database/$directory/$pdb* tmp_master
fi
cd tmp_master
for file in *
do
  gunzip $file
  if (( $INCLUDE_H == 1 ))
  then
    fileroot="${file%.pdb.gz}"
    reduce -Trim $fileroot.pdb > $fileroot\_noH.pdb
    reduce $fileroot\_noH.pdb > $fileroot\_H.pdb
    echo "Fileroot = $fileroot"
    rm $fileroot.pdb
    rm $fileroot\_noH.pdb
  else
    # Remove hydrogens in case structure has tehm 
    fileroot="${file%.pdb.gz}"
    reduce -Trim $fileroot.pdb > $fileroot\_noH.pdb
    mv $fileroot\_noH.pdb $fileroot.pdb
  fi

done
cd ..
if (( $INCLUDE_H == 1 ))
then
  OPTIONS="-includeH"
else
  OPTIONS=""
fi
echo ./compute_msms_buried_surface_area.py $OPTIONS tmp_master/*
./compute_msms_buried_surface_area.py $OPTIONS tmp_master/*
#rm tmp_master/*
