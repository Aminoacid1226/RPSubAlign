# 2024.01.10 Random smiles are generated according to cano
import random
from rdkit import Chem

def randomize_smiles(smiles, random_type="unrestricted"):
    """
    Given the SMILES of a molecule, a random SMILES is returned.
    :param smiles: The entered SMILES string
    :param random_type: Type of randomization ("unrestricted" or "restricted")
    :return: Random SMILES string for the same molecule or None if the molecule is invalid
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    if random_type == "unrestricted":
        return Chem.MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=False)
    elif random_type == "restricted":
        new_atom_order = list(range(mol.GetNumAtoms()))
        random.shuffle(new_atom_order)
        random_mol = Chem.RenumberAtoms(mol, newOrder=new_atom_order)
        return Chem.MolToSmiles(random_mol, canonical=False, isomericSmiles=False)
    else:
        raise ValueError(f"type '{random_type}' invalid")

def read_smiles_from_file(file_type,file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()
    if file_type == "tgt":
        # tgt files use this sentence
        smiles_list=[''.join(line.strip().split()) for line in lines]
    if file_type == "src":
        # USPTO-50k-src files use this sentence
        # smiles_list = [''.join(line.strip().split()[1:]) for line in lines]
        smiles_list = [''.join(line.strip().split()) for line in lines]
    classes = [line.strip().split()[0] for line in lines]
    return classes, smiles_list

def write_randomized_smiles_to_file(output_file_path, classes, smiles_list):
    with open(output_file_path, "w") as output_file:
        for cls, smiles in zip(classes, smiles_list):
            random_smiles = randomize_smiles(smiles)
            if random_smiles is not None:
                spaced_random_smiles = " ".join(random_smiles)
                output_file.write(f"{cls} {spaced_random_smiles}\n")

file_type = "src"
input_file_path = "../data/SMILES representation/USPTO-MIT/cano/src-test.txt"
output_file_path = f"../data/SMILES representation/USPTO-MIT/Random/{file_type}-test.txt"

classes, smiles_list = read_smiles_from_file(file_type,input_file_path)

write_randomized_smiles_to_file(output_file_path, classes, smiles_list)
