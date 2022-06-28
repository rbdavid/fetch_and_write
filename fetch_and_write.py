
"""
Fetch and save to file PDB structures
"""

import argparse
import mmtf
import MDAnalysis
import warnings
from pathlib import Path
from joblib import Parallel, delayed


#######################################
# FUNCTIONS
#######################################

def fetch_write_single_structure(pdb_meta_lst, output_directory = './'):
    ''' grabs PDBIDs and optionally a chainID to fetch the structure data and write out to file
    INPUT:
        pdb_meta_lst: list, len of 1 or 2 (elements beyond will be ignored); 0th element is the PDBID, the optional 1st element is chainID
    OUTPUT:
        returns the pdb_meta_lst appended with the filename (if successfully fetched and wrote to file) or "failed" (if fetching the PDBID or writing the .pdb file fails)
    '''
    output_directory = Path(output_directory)
    pdbID = pdb_meta_lst[0]
    outFileName = pdbID
    # check if chainID is included in the input list
    if len(pdb_meta_lst) > 1:
        chainID = pdb_meta_lst[1]
        outFileName += '_' + chainID
    else:
        chainID = None

    # Pull the structure from the RCSB database. Store as mmtf object
    try:
        pdb = mmtf.fetch(pdbID)
    except Exception as e:
        print(f"The PDBID ({pdbID}) was not able to be fetched. Error: {e}")
        pdb_meta_lst.append(f"The PDBID ({pdbID}) was not able to be fetched. Error: {e}")
        return pdb_meta_lst
    
    # check if chainID is not None AND the chainID is in the list of chains stored in the mmtf object
    if chainID and chainID not in pdb.chain_id_list:
        print(f"The expected chain ID ({chainID}) is not found in the PDB's chain list")
        pdb_meta_lst.append(f"The expected chain ID ({chainID}) is not found in the PDB's chain list")
        return pdb_meta_lst
    
    # loading the mmtf structure object into an MDAnalysis universe object
    u = MDAnalysis.Universe(pdb)
    # make chain selection if chainID is not None
    if chainID:
        try:
            # a bit disappointing that RCSB and MDAnalysis use different naming conventions for the same thing, chainID
            sel = u.select_atoms(f'segid {chainID}')
        except Exception as e:
            print(f"The expected chain ID ({chainID}) does not return an atom selection, check formatting of the chainID. Error: {e}")
            pdb_meta_lst.append(f"The expected chain ID ({chainID}) does not return an atom selection, check formatting of the chainID. Error: {e}")
            return pdb_meta_lst
    # make all selection if no chainID used
    else:
        sel = u.select_atoms('all')
    
    # write atom selection to file, if possible
    try:
        with warnings.catch_warnings():
            # ignore some annoying warnings from sel.write line due to missing information (chainIDs, elements, and record_types). 
            warnings.simplefilter('ignore',UserWarning)
            # test for the potential for writing a pdb file failing due to size limitations
            # return 0 (denoting failure) 
            if sel.n_atoms > 99999:
                print(f'Number of atoms ({sel.n_atoms}) is too large for the pdb file format; need to store in some other file format.')
                pdb_meta_lst.append(f'Number of atoms ({sel.n_atoms}) is too large for the pdb file format; need to store in some other file format.')
                return pdb_meta_lst
            # write the file out to a name and pass the name back to main
            else:
                sel.write(output_directory / f'{outFileName}.pdb')
                pdb_meta_lst.append(outFileName)
                return pdb_meta_lst
    except Exception as e:
        print(f"Writing to file ({outFileName}) failed. Error: {e}")
        pdb_meta_lst.append(f"Writing to file ({outFileName}) failed. Error: {e}")
        return pdb_meta_lst


#######################################
# MAIN
#######################################
if __name__ == '__main__':
    # read arguments
    parser = argparse.ArgumentParser(description='Read in a list of PDB IDs and subsequent chain IDs (optional) that will be fetched and written to file.')
    parser.add_argument('--pdbid-list-file','-pdbIDs',required=True,type=str,help='string, path to the file within which the PDB IDs to be gathered and saved are stored; expected format: "PDBID   [chainID]\n" where the chainID is optional.')
    parser.add_argument('--out-file-directory','-outdir',required=True,type=str,help='string, path to the directory within which the .pdb files will be saved.')
    parser.add_argument('--max-threads','-c',required=True,type=int,help='number of threads to be made available to run this code')
    args = parser.parse_args()

    args.out_file_directory = Path(args.out_file_directory)

    # read the pdbid list file, gathering lines
    with open(args.pdbid_list_file,'r') as pdbid_file:
        lines = pdbid_file.readlines()
    # split the lines, create a list of lists with each sublist having length 1 or 2; element 0 is the PDBID, element 1 (if present) is the chainID; don't worry about empty lines here
    lines = [line.split() for line in lines if line.strip() != '']

    # run the gather script across the number of cores denoted in the arguments
    results = Parallel(n_jobs=args.max_threads,prefer='threads')(
            delayed(fetch_write_single_structure)(pdbid_list,output_directory = args.out_file_directory) for pdbid_list in lines)

    # write the results list out to a log file
    with open(args.out_file_directory / 'fetching.log','w') as log_file:
        for result in results:
            string = ' '.join(result) + '\n'
            log_file.write(string)


