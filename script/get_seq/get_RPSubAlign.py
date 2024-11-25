# -*- coding: utf-8 -*-
# SA
"""
Created on Tue Aug 30 16:53:44 2022

@author: user
"""

import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS
import random

aug = 40
tfile=open('../data/SMILES representation/UPSTO-50K/cano/tgt-train.txt')
sfile=open('../data/SMILES representation/UPSTO-50K/cano/src-train.txt')

tfw=open('../data/SMILES representation/UPSTO-50k/SA_40x/tgt-train.txt','w')
sfw=open('../data/SMILES representation/UPSTO-50k/SA_40x/src-train.txt','w')

lines1=tfile.readlines()
lines2=sfile.readlines()

tsmiles=[''.join(line.strip().split()) for line in lines1]
ssmiles=[''.join(line.strip().split()[1:]) for line in lines2]
classes=[line.strip().split()[0] for line in lines2]

for i in range(len(tsmiles)):
    
    smol=Chem.MolFromSmiles(ssmiles[i])
    tsmiles_list=tsmiles[i].split('.')
    tmols_list=[Chem.MolFromSmiles(s) for s in tsmiles_list]
    cl=classes[i]
    res=[]
    for mol in tmols_list:
        res.append(rdFMCS.FindMCS([smol,mol]))
    anums=[r.numAtoms for r in res]
    index=anums.index(max(anums))
    patt=res[index].queryMol
    
    a1=smol.GetSubstructMatch(patt)
    a2=tmols_list[index].GetSubstructMatch(patt)
    
    index_list=list(range(len(a1)))
    
    for j in range(aug):
        
        random.shuffle(index_list)    
        a1=[a1[n] for n in index_list]
        a2=[a2[n] for n in index_list]
    
        order_s=a1+list(set(list(range(smol.GetNumAtoms()))).difference(set(a1)))
        random_mol = Chem.RenumberAtoms(smol, newOrder=order_s) 
        ssmiles_new=Chem.MolToSmiles(random_mol, canonical=False)
        sfw.write(cl+' '+' '.join(ssmiles_new)+'\n')
    
        order_t1=a2+list(set(list(range(tmols_list[index].GetNumAtoms()))).difference(set(a2)))
        random_mol = Chem.RenumberAtoms(tmols_list[index], newOrder=order_t1) 
        tsmiles_new1=Chem.MolToSmiles(random_mol, canonical=False)
        tsmiles_new_other=''
        
        if len(tmols_list)>1:
            tmols_other=[tmols_list[n] for n in list(range(len(tmols_list))) if n !=index]            
            for m in tmols_other:
                tsmiles_new_other+='.'
                tsmiles_new_other+=Chem.MolToSmiles(m,doRandom=True)
        
        tsmiles_new=tsmiles_new1+tsmiles_new_other
        tfw.write(' '.join(tsmiles_new)+'\n')
        

tfw.close()
sfw.close()    
  