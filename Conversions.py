# Imports
from rdkit import Chem
import numpy as np
import selfies as sf
import networkx as nx
import matplotlib.pyplot as plt
from urllib.request import urlopen
from urllib.parse import quote
import nltk as nk
nk.download('punkt')


## Molecular Name to .mol
def name_to_inchi(name):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(name) + '/stdinchi'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Did not work'

def inchi_to_mol(inchi):
    # Create the .mol file from the InChI code
    mol = Chem.inchi.MolFromInchi(inchi, sanitize=True, treatWarningAsError=False)
    # Explicitly state the hydrogen atoms
    mol = Chem.AddHs(mol)
    return mol


## .mol and SMILES and SELFIES
def mol_to_smiles(mol):
    # Remove the hydrogen atoms for SMILES
    mol = Chem.RemoveHs(mol)
    # Convert the molecule to SMILES code
    smiles = Chem.MolToSmiles(mol)
    return smiles

def smiles_to_mol(smiles):
    # Create the molecular file
    mol = Chem.MolFromSmiles(smiles, sanitize=True, treatWarningAsError=False)
    # Explicitly state the hydrogen atoms
    mol = Chem.AddHs(mol)
    return mol

def smiles_to_selfies(smiles)
    selfies = sf.encoder(smiles)
    return selfies

def selfies_to_smiles(selfies)
    smiles = sf.decoder(smiles)
    return smiles


## Create graphics
def generate_adjacency_matrix(mol):
    # Get number of atoms from the .mol file
    num_atoms = mol.GetNumAtoms()
    adjacency_matrix = np.zeros((num_atoms, num_atoms), dtype=float)

    # Iterate through each pair of atom to find bond type
    for bond in mol.GetBonds():
        atom1_index = bond.GetBeginAtomIdx()
        atom2_index = bond.GetEndAtomIdx()
        bond_type = bond.GetBondTypeAsDouble()

        adjacency_matrix[atom1_index][atom2_index] = bond_type
        adjacency_matrix[atom2_index][atom1_index] = bond_type

    return adjacency_matrix

def generate_vertices(mol):
    # Extract the list of atoms
    atoms = mol.GetAtoms()

    # Get the element symbols from the atoms
    element_list = [atom.GetSymbol() for atom in atoms]

    return element_list

def build_graph(mol):
    adj_matrix = generate_adjacency_matrix(mol)
    g = from_numpy_matrix(adj_matrix)
    
    return g


## Storing More Properties
def mol_to_vertices_list(mol):
    vertices_list = []

    for atom in mol.GetAtoms():
        vertex_data = {
            'atomic_num': atom.GetAtomicNum(),
            'formal_charge': atom.GetFormalCharge(),
            'chiral_tag': atom.GetChiralTag(),
            'hybridization': atom.GetHybridization(),
            'num_explicit_hs': atom.GetNumExplicitHs(),
            'is_aromatic': atom.GetIsAromatic()
        }
        vertices_list.append((atom.GetIdx(), vertex_data))

    return vertices_list

def mol_to_matrix(mol):
    num_atoms = mol.GetNumAtoms()
    adjacency_matrix = np.zeros((num_atoms, num_atoms), dtype=int)

    for bond in mol.GetBonds():
        begin_atom_idx = bond.GetBeginAtomIdx()
        end_atom_idx = bond.GetEndAtomIdx()
        adjacency_matrix[begin_atom_idx, end_atom_idx] = 1
        adjacency_matrix[end_atom_idx, begin_atom_idx] = 1

    return adjacency_matrix


## Getting SMILES From Graph
def nx_to_mol(vertices, matrix):
    mol = Chem.RWMol()
    for id in range(0, len(vertices)):
        symbol = vertices[id]
        atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(symbol)
        a=Chem.Atom(atomic_num)
        molIdx = mol.AddAtom(a)
        
    num_rows, num_cols = matrix.shape
    
    for row in range(num_rows):
        for column in range(num_cols):
            if column < row:
                value = matrix[row, column]
                if value == 1.0:
                    bond = Chem.BondType.SINGLE
                    mol.AddBond(row, column, bond)
                elif value == 1.5:
                    bond = Chem.BondType.AROMATIC
                    mol.AddBond(row, column, bond)
                elif value == 2.0:
                    bond = Chem.BondType.DOUBLE
                    mol.AddBond(row, column, bond)
                elif value == 3.0:
                    bond = Chem.BondType.TRIPLE
                    mol.AddBond(row, column, bond)
                    
    Chem.SanitizeMol(mol)
    return mol


## Getting Additional Properties
def nx_to_mol(vertices, matrix):
    mol = Chem.RWMol()
    atomic_nums = []
    chiral_tags = []
    formal_charges = []
    aromaticity = []
    hybridizations = []
    explicitHs = []
    for node in vertices:
        a=Chem.Atom(atomic_nums[node])
        a.SetChiralTag(chiral_tags[node])
        a.SetFormalCharge(formal_charges[node])
        a.SetIsAromatic(aromaticity[node])
        a.SetHybridization(hybridizations[node])
        a.SetNumExplicitHs(explicitHs[node])
        molIdx = mol.AddAtom(a)

    # FIX THE BOND STUFF: EITHER STORE SEPERATE LIST OR FIGURE OUT RDKIT STUFF
    # bond_types = nx.get_edge_attributes(G, 'bond_type')
    # num_rows, num_cols = matrix.shape
    
    for row in range(num_rows):
        for column in range(num_cols):
            if column < row:
                value = matrix[row, column]
                if value == 1.0:
                    bond = Chem.BondType.SINGLE
                    mol.AddBond(row, column, bond)
                elif value == 1.5:
                    bond = Chem.BondType.AROMATIC
                    mol.AddBond(row, column, bond)
                elif value == 2.0:
                    bond = Chem.BondType.DOUBLE
                    mol.AddBond(row, column, bond)
                elif value == 3.0:
                    bond = Chem.BondType.TRIPLE
                    mol.AddBond(row, column, bond)
     
    Chem.SanitizeMol(mol)
    return mol



