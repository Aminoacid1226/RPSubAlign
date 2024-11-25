# -*- coding: UTF-8 -*-
"""
@Date    ：2024-05-31 20:39 
@Software: PyCharm
@Project ：SA-20240403-matric 
@File    ：matrics_new_3.0.py
@Author  ：AminoAcid
@Description: 0531 Iterate over all TXTS in the folder
"""
import os
from rdkit import Chem
import numpy
from numpy import mean, std
import datetime

folder_path = "pred" # where save model output result

file_output_template = "metrics/test-{}.txt"  # metrics output file path

file_test = "SMILES representation/data/UPSTO-50k/cano/tgt-test.txt"  # standard target

current_datetime = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def canonicalize_smiles(s):
    mol = Chem.MolFromSmiles(s)
    if mol is not None:
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    else:
        return None


def get_test_smiles(file_test, cano=True):
    with open(file_test, 'r') as f:
        test_smiles = [''.join(line.strip().split()) for line in f.readlines()]
    if cano:
        test_smiles = [canonicalize_smiles(s) for s in test_smiles]
    return test_smiles


def get_pred_smiles(file_pred, cano=True):
    with open(file_pred, 'r') as f:
        pred_smiles = [''.join(line.strip().split()) for line in f.readlines()]
    if cano:
        pred_smiles = [canonicalize_smiles(s) for s in pred_smiles]
    return pred_smiles


def get_pred_dict(smiles, best_num):
    pred = {i: [] for i in range(1, best_num + 1)}
    for i, smile in enumerate(smiles):
        pred[(i % best_num) + 1].append(smile)
    return pred


def get_rank(pred_dict, test_smiles, best_num):
    rank = [0] * len(test_smiles)
    for i in range(len(test_smiles)):
        target = test_smiles[i]
        for j in range(1, best_num + 1):
            if target == pred_dict[j][i]:
                rank[i] = j
                break
    return rank


def get_rank_maxfrag(pred_dict, tsmiles, best_num):
    rank = [0] * len(tsmiles)

    for i in range(len(tsmiles)):
        target_list = tsmiles[i].split(".")
        if len(target_list) > 1:
            target = max(target_list, key=len)
        else:
            target = target_list

        for j in range(1, best_num + 1):
            pred_list = pred_dict[j][i]
            if pred_list is None:
                continue

            if not isinstance(pred_list, list):
                pred_list = pred_list.split(".")

            if len(pred_list) > 1:
                pred = max(pred_list, key=len)
            else:
                pred = pred_list
            # print(pred)
            if target == pred:
                rank[i] = j
                break

    return rank


best_num = 10
aug = 10

num_pre_smiles = best_num * aug

test_smiles = get_test_smiles(file_test, cano=True)

for file_name in os.listdir(folder_path):
    if file_name.endswith(".txt"):
        file_pred = os.path.join(folder_path, file_name)
        pred_smiles = get_pred_smiles(file_pred, cano=True)

        total = len(test_smiles)

        pred_smiles2 = numpy.array(pred_smiles).reshape(total, aug, best_num)

        pred_smiles3 = [[] for _ in range(aug)]
        for i in range(aug):
            for j in range(total):
                pred_smiles3[i].extend(list(pred_smiles2[j][i]))

        dict_recovery = {i: [] for i in range(1, best_num + 1)}
        dict_valid = {i: [] for i in range(1, best_num + 1)}
        dict_maxfrag = {i: [] for i in range(1, best_num + 1)}

        for k in range(1, aug + 1):
            smiles = pred_smiles3[k - 1]

            pred_dict = get_pred_dict(smiles, best_num=best_num)
            rank = get_rank(pred_dict, test_smiles, best_num=best_num)
            rank_maxfrag = get_rank_maxfrag(pred_dict, test_smiles, best_num=best_num)

            correct = 0
            correct_maxfrag = 0
            pred = []
            for i in range(1, best_num + 1):
                correct += rank.count(i)
                pred.extend(pred_dict[i])
                valid = len(pred) - pred.count(None)

                dict_recovery[i].append(correct / total)
                dict_valid[i].append(valid / (total * i))
                correct_maxfrag += rank_maxfrag.count(i)
                dict_maxfrag[i].append(correct_maxfrag / total)

        output_file = file_output_template.format(file_name.split('.')[0])
        with open(output_file, "w") as f:
        
            f.write(f"Output File: {output_file}\n")
            f.write(f"Script Start Time: {current_datetime}\n")
            f.write(f"file_pred: {file_pred}\n\n")

      
            f.write(f"----------------------------- recovery ---------------------------------\n")
            for i in range(1, aug + 1):
                f.write("recovery_test " + str(i) + "\n")
                for m in range(1, best_num + 1):
                    f.write("top " + str(m) + " : " + str(dict_recovery[m][i - 1]) + "\n")

            f.write(f"----------------------------- maxfrag ---------------------------------\n")
            for i in range(1, aug + 1):
                f.write("maxfrag_test " + str(i) + "\n")
                for m in range(1, best_num + 1):
                    f.write("top " + str(m) + " : " + str(dict_maxfrag[m][i - 1]) + "\n")

            f.write(f"----------------------------- valid ----------------------------------\n")
            for i in range(1, aug + 1):
                f.write("valid_test " + str(i) + "\n")
                for m in range(1, best_num + 1):
                    f.write("top " + str(m) + " : " + str(dict_valid[m][i - 1]) + "\n")

            f.write(f"----------------------------- total ---------------------------------\n")
            for i in range(1, best_num + 1):
                output_line = 'top_{} : recovery {:.2f}±{:.2f}% || maxfrag {:.2f}±{:.2f}% || valid {:.2f}±{:.2f}%\n' \
                    .format(i, mean(dict_recovery[i]) * 100, std(dict_recovery[i]) * 100,
                            mean(dict_maxfrag[i]) * 100, std(dict_maxfrag[i]) * 100,
                            mean(dict_valid[i]) * 100, std(dict_valid[i]) * 100)
                f.write(output_line)

print("Done!")
