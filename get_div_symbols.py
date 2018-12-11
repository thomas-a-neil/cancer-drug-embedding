from rdkit import Chem


div5_symbols = {'As', 'Br', 'C', 'Cl', 'F', 'I', 'N', 'O', 'P', 'S'}


def get_all_molecules_with_allowed_symbols(all_mols, allowed_symbols):
    filtered_molecules = []
    for mol in all_mols:
        if mol is None:
            continue
        mol_symbols = set()
        for j in range(mol.GetNumAtoms()):
            mol_symbols.add(mol.GetAtomWithIdx(j).GetSymbol())
        if mol_symbols <= allowed_symbols:
            filtered_molecules.append(mol)
    return filtered_molecules


def get_all_chemical_symbols(all_mols):
    symbols = set()
    for i in range(len(all_mols)):
        m = all_mols[i]
        if m is None:
            continue
        for j in range(m.GetNumAtoms()):
            if '*' == m.GetAtomWithIdx(j).GetSymbol():
                print(Chem.MolToSmiles(m))
                print(i+1)
                break
            symbols.add(m.GetAtomWithIdx(j).GetSymbol())
    print(symbols)
    return symbols


