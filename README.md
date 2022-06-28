# fetch_and_write
Code to fetch a list of PDBID structure and write it out as a PDB file. 

## Usage example: 
python3 fetch_and_write.py --pdbid-list-file tests/test.txt --out-file-directory tests/ --max-threads 2

### Caveats/Concerns
MDAnalysis confuses segment IDs and chain IDs, relative to the columns in a RCSB .pdb file. This results in a change in placement of that structural information (chainID) in .pdb files written by MDAnalysis. A bit unfortunate and potentially fixable. More work to be done.  


