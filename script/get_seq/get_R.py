# 2024.01.10 
import numpy as np
import re
import random

#R_SMILES
def get_root_id(mol,root_map_number):
    root = -1
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomMapNum() == root_map_number:
            root = i
            break
    return root
def get_cano_map_number(smi,root=-1):
    atommap_mol = Chem.MolFromSmiles(smi)
    canonical_mol = Chem.MolFromSmiles(clear_map_canonical_smiles(smi,root=root))
    cano2atommapIdx = atommap_mol.GetSubstructMatch(canonical_mol)
    correct_mapped = [canonical_mol.GetAtomWithIdx(i).GetSymbol() == atommap_mol.GetAtomWithIdx(index).GetSymbol() for i,index in enumerate(cano2atommapIdx)]
    atom_number = len(canonical_mol.GetAtoms())
    if np.sum(correct_mapped) < atom_number or len(cano2atommapIdx) < atom_number:
        cano2atommapIdx = [0] * atom_number
        atommap2canoIdx = canonical_mol.GetSubstructMatch(atommap_mol)
        if len(atommap2canoIdx) != atom_number:
            return None
        for i, index in enumerate(atommap2canoIdx):
            cano2atommapIdx[index] = i
    id2atommap = [atom.GetAtomMapNum() for atom in atommap_mol.GetAtoms()]

    return [id2atommap[cano2atommapIdx[i]] for i in range(atom_number)]
def smi_tokenizer(smi):
    pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return ' '.join(tokens)
def clear_map_canonical_smiles(smi, canonical=True, root=-1):
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        for atom in mol.GetAtoms():
            if atom.HasProp('molAtomMapNumber'):
                atom.ClearProp('molAtomMapNumber')
        return Chem.MolToSmiles(mol, isomericSmiles=True, rootedAtAtom=root, canonical=canonical)
    else:
        return smi
def get_retro_rsmiles(data):
    
    pt = re.compile(r':(\d+)]') 
    product = data['product']
    reactant = data['reactant']
    augmentation = data['augmentation']
    pro_mol = Chem.MolFromSmiles(product)
    rea_mol = Chem.MolFromSmiles(reactant)
    """checking data quality"""
    rids = sorted(re.findall(pt, reactant))
    pids = sorted(re.findall(pt, product))
    return_status = {
        "status":0,
        "src_data":[],
        "tgt_data":[],
        "edit_distance":0,
    }
    pro_atom_map_numbers = list(map(int, re.findall(r"(?<=:)\d+", product))) 

    reactant = reactant.split(".")
    reversable = False  # no shuffle
    # augmentation = 100
    if augmentation == 999:
        product_roots = pro_atom_map_numbers
        times = len(product_roots)
    else:
        product_roots = [-1]
        # reversable = len(reactant) > 1

        max_times = len(pro_atom_map_numbers)
        times = min(augmentation, max_times)
        if times < augmentation:  # times = max_times
            product_roots.extend(pro_atom_map_numbers)
            product_roots.extend(random.choices(product_roots, k=augmentation - len(product_roots)))
        else:  # times = augmentation
            while len(product_roots) < times:
                product_roots.append(random.sample(pro_atom_map_numbers, 1)[0])
                # pro_atom_map_numbers.remove(product_roots[-1])
                if product_roots[-1] in product_roots[:-1]:
                    product_roots.pop()
        times = len(product_roots)

        assert times == augmentation
        if reversable:
            times = int(times / 2)
    # candidates = []
    for k in range(times):
        pro_root_atom_map = product_roots[k]
        pro_root = get_root_id(pro_mol, root_map_number=pro_root_atom_map)
        cano_atom_map = get_cano_map_number(product, root=pro_root)
        if cano_atom_map is None:
            return_status["status"] = "error_mapping"
            return return_status
        pro_smi = clear_map_canonical_smiles(product, canonical=True, root=pro_root)
        aligned_reactants = []
        aligned_reactants_order = []
        rea_atom_map_numbers = [list(map(int, re.findall(r"(?<=:)\d+", rea))) for rea in reactant]
        used_indices = []
        for i, rea_map_number in enumerate(rea_atom_map_numbers):
            for j, map_number in enumerate(cano_atom_map):
                # select mapping reactans
                if map_number in rea_map_number:
                    rea_root = get_root_id(Chem.MolFromSmiles(reactant[i]), root_map_number=map_number)
                    rea_smi = clear_map_canonical_smiles(reactant[i], canonical=True, root=rea_root)
                    aligned_reactants.append(rea_smi)
                    aligned_reactants_order.append(j)
                    used_indices.append(i)
                    break
        sorted_reactants = sorted(list(zip(aligned_reactants, aligned_reactants_order)), key=lambda x: x[1])
        aligned_reactants = [item[0] for item in sorted_reactants]
        reactant_smi = ".".join(aligned_reactants)
        product_tokens = smi_tokenizer(pro_smi)
        reactant_tokens = smi_tokenizer(reactant_smi)

        return_status['src_data'].append(product_tokens)
        return_status['tgt_data'].append(reactant_tokens)

        if reversable:
            aligned_reactants.reverse()
            reactant_smi = ".".join(aligned_reactants)
            product_tokens = smi_tokenizer(pro_smi)
            reactant_tokens = smi_tokenizer(reactant_smi)
            return_status['src_data'].append(product_tokens)
            return_status['tgt_data'].append(reactant_tokens)
    assert len(return_status['src_data']) == data['augmentation']
    edit_distances = []
    for src,tgt in zip(return_status['src_data'],return_status['tgt_data']):
        edit_distances.append(textdistance.levenshtein.distance(src.split(),tgt.split()))
    return_status['edit_distance'] = np.mean(edit_distances)
    return return_status


def generate_retro_rsmiles(src_file, tgt_file, output_src_file, output_tgt_file):
    src = open(src_file)
    src_lines = src.readlines()
    src_smiles_list = [''.join(line.strip().split()[1:]) for line in src_lines]
    classes = [line.strip().split()[0] for line in src_lines]
    src_mol_list = [Chem.MolFromSmiles(smiles) for smiles in src_smiles_list]

    src_map_list = []
    for m in src_mol_list:
       
        [a.SetProp('Molatommapnumber', str(i + 1)) for (i, a) in enumerate(m.GetAtoms())]
        [a.SetAtomMapNum(i) for (i, a) in enumerate(m.GetAtoms(), start=1)]
        src_map_list.append(Chem.MolToSmiles(m))

    tgt = open(tgt_file)
    tgt_lines = tgt.readlines()
    tgt_smiles_list = [''.join(line.strip().split()) for line in tgt_lines]
    tgt_mol_list = [Chem.MolFromSmiles(smiles) for smiles in tgt_smiles_list]

    tgt_map_list = []
    for m in tgt_mol_list:
   
        [a.SetProp('Molatommapnumber', str(i + 1)) for (i, a) in enumerate(m.GetAtoms())]
        [a.SetAtomMapNum(i) for (i, a) in enumerate(m.GetAtoms(), start=1)]
    
        tgt_map_list.append(Chem.MolToSmiles(m))

    with open(output_src_file, "w") as f, open(output_tgt_file, "w") as fb:
        for i in range(len(src_mol_list)):
            data = {
                'product': src_map_list[i],
                'reactant': tgt_map_list[i],
                'augmentation': 5
            }

            R_Smiles = get_retro_rsmiles(data)

            f.write(classes[i] + " " + " ".join(R_Smiles["src_data"]) + '\n')
            fb.write(" ".join(R_Smiles["tgt_data"]) + '\n')


generate_retro_rsmiles(
    src_file=r"../data/SMILES representation/USPTO-50K/cano/src-test.txt",
    tgt_file=r"../data/SMILES representation/USPTO-50K/cano/tgt-test.txt",
    output_src_file=r"../SMILES representation/data/USPTO-50K/R/src_test.txt",
    output_tgt_file=r"../SMILES representation/data/USPTO-50K/R/tgt_test.txt"
)