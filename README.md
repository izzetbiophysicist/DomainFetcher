# cut_protein_domains
Fetches structures from the PDB and cuts them into single domain structures according to PFAM annotation. The script demands a PDB-pfam mapping that can be found at http://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/  and a list of PDBs in csv  format containing the PDB codes in a single column. 

# Usage
Rscript cut_domains.R pdb_list.csv

The pdb_list.csv file should follow the format of the example


A "domains" directory will be created in the current directory, and within it, a folder for each domain. The domains extracted from each PDB structure will be in their correspondent folders. The script will also create files listing PDBs that couldn't be downloades and PDBs lacking author's domain mapping in the pdb-pfam mapping file. 

