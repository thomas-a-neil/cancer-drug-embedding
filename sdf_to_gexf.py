from multiprocessing import Pool
import os
import pickle

import networkx as nx
from rdkit import Chem
from rdkit.Chem import Draw

sdf_filepath = 'data/Chem2D_Jun2016.sdf'
pickle_filepath = 'data/all_mols.pkl'
div5_sdf_filepath = 'data/Div5_2DStructures_Oct2014.sdf'
training_data_dir = 'training_data'

div5_symbols = {'As', 'Br', 'C', 'Cl', 'F', 'I', 'N', 'O', 'P', 'S'}

symbol_to_index = {'As': 0,
                   'Br': 1,
                   'C': 2,
                   'Cl': 3,
                   'F': 4,
                   'I': 5,
                   'N': 6,
                   'O': 7,
                   'P': 8,
                   'S': 9}

index_to_symbol = {0: 'As',
                   1: 'Br',
                   2: 'C',
                   3: 'Cl',
                   4: 'F',
                   5: 'I',
                   6: 'N',
                   7: 'O',
                   8: 'P',
                   9: 'S'}


# for reconstructing dataset later, much faster from pickle, supposedly
# -rw-rw-r-- 1 nthomas nthomas 166M Dec  5 20:44 all_mols.pkl
# -rwx------ 1 nthomas nthomas 678M Feb  9  2017 Chem2D_Jun2016.sdf*
def sdf_to_pickle(sdf_filepath):
    suppl = Chem.SDMolSupplier(sdf_filepath)
    all_mols = [mol for mol in suppl]
    with open(pickle_filepath, 'wb') as f:
        pickle.dump(all_mols, f)


def load_all_mols_from_pickle():
    with open(pickle_filepath, 'rb') as f:
        all_mol = pickle.load(f)
    return all_mol


def get_all_chemical_symbols(all_mols):
    symbols = set()
    for mol in all_mols:
        if mol is None:
            continue
        for j in range(mol.GetNumAtoms()):
            symbols.add(mol.GetAtomWithIdx(j).GetSymbol())
    print(symbols)
    return symbols


def get_div_symbols():
    all_mols = Chem.SDMolSupplier(div5_sdf_filepath)
    symbols = get_all_chemical_symbols(all_mols)
    return symbols


def get_all_molecules_with_allowed_symbols(mols, nsc_ids, allowed_symbols):
    """
    Filtering

    TODO nthomas: move this into pandas, probably.
    """
    filtered_molecules = []
    filtered_nsc_ids = []
    for i, mol in enumerate(mols):
        if mol is None:
            continue
        mol_symbols = set()
        for j in range(mol.GetNumAtoms()):
            mol_symbols.add(mol.GetAtomWithIdx(j).GetSymbol())
        if mol_symbols <= allowed_symbols and '.' not in Chem.MolToSmiles(mol):
            filtered_molecules.append(mol)
            filtered_nsc_ids.append(nsc_ids[i])
    return filtered_molecules, filtered_nsc_ids


def convert_molecule_to_graph(mol):
    G = nx.Graph()
    for i in range(mol.GetNumAtoms()):
        symbol = mol.GetAtomWithIdx(i).GetSymbol()
        G.add_node(i, value=symbol_to_index[symbol])
    for i in range(mol.GetNumBonds()):
        start_idx = mol.GetBondWithIdx(i).GetBeginAtomIdx()
        end_idx = mol.GetBondWithIdx(i).GetEndAtomIdx()
        # TODO nthomas: bond types?
        # mol.GetBondBetweenAtoms(start_idx, end_idx).GetBondType()
        G.add_edge(start_idx, end_idx)
    return G


def write_molecule_to_gexf(mol, path):
    graph = convert_molecule_to_graph(mol)
    nx.write_gexf(graph, path)


def get_nsc_ids_from_sdf(filepath, delimiter='$$$$\n'):
    nsc_ids = []
    # some bytes in the .sdf file cannot be converted to unicode
    # errors='replace' replaces malformed bytes with a ?
    # we hope that the malformed bytes are not characters we need for NSC id
    with open(filepath, 'r', errors='replace') as f:
        is_nsc = True
        for line in f:
            if is_nsc:
                nsc_ids.append(line.strip())
                is_nsc = False
            if line == delimiter:
                is_nsc = True
    return nsc_ids


def generate_gexf_training_data(directory):
    # sdf_to_pickle(sdf_filepath)
    all_mols = load_all_mols_from_pickle()
    all_nsc_ids = get_nsc_ids_from_sdf(sdf_filepath)
    # TODO nthomas: switch to generator so that all mols don't
    filtered_mols, filtered_nsc_ids = get_all_molecules_with_allowed_symbols(all_mols, all_nsc_ids, div5_symbols)

    # a wrapper so I can use pool.map
    def pool_wrapper_write_molecule_to_dir(directory):

        def pool_wrapper_write_molecule(nsc_id_mol_tuple):
            nsc_id, mol = nsc_id_mol_tuple
            path = os.path.join(directory, nsc_id + '.gexf')
            write_molecule_to_gexf(mol, path)

        return pool_wrapper_write_molecule

    pool = Pool()
    write_molecule = pool_wrapper_write_molecule_to_dir(directory)
    pool.map(write_molecule, zip(filtered_nsc_ids, filtered_mols))


def write_mol_png(nsc_id, mol, directory):
        path = os.path.join(directory, nsc_id + '.png')
        Draw.MolToFile(mol, path)


def get_all_pngs():
    all_mols = load_all_mols_from_pickle()
    all_nsc_ids = get_nsc_ids_from_sdf(sdf_filepath)
    filtered_mols, filtered_nsc_ids = get_all_molecules_with_allowed_symbols(all_mols, all_nsc_ids, div5_symbols)

    molecule_image_dir = '/data/molecule_images'

    pool = Pool()
    pool.starmap(write_mol_png, zip(filtered_nsc_ids, filtered_mols, [molecule_image_dir]*len(filtered_nsc_ids)))


def get_num_atoms_per_molecule(all_mols):
    all_nsc_ids = get_nsc_ids_from_sdf(sdf_filepath)
    filtered_mols, filtered_nsc_ids = get_all_molecules_with_allowed_symbols(all_mols, all_nsc_ids, div5_symbols)
    num_atoms = []
    for mol in filtered_mols:
        num_atoms.append(mol.GetNumAtoms())

    with open('num_atoms_to_nsc.csv', 'w') as f:
        for num_atoms, nsc_id in zip(num_atoms, filtered_nsc_ids):
            f.write('{}, {}\n'.format(num_atoms, nsc_id))
