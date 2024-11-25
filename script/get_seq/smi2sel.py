import selfies as sf
import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS
import random
import re

'''
tfile=open('new_data/smiles/random_1x/tgt_train.txt')
sfile=open('new_data/smiles/random_1x/src_train.txt')

tfw=open('new_data/random_1x/tgt_train.txt','w')
sfw=open('new_data/random_1x/src_train.txt','w')

lines1=tfile.readlines()
lines2=sfile.readlines()
str_all=''
tsmiles=[''.join(line.strip().split()) for line in lines1]
ssmiles=[''.join(line.strip().split()[1:]) for line in lines2]
classes=[line.strip().split()[0] for line in lines2]
for i in range(len(tsmiles)):

    str_sf=sf.encoder(ssmiles[i])
    tgt_sf=sf.encoder(tsmiles[i])
    cl = classes[i]
    str_tokens = list(sf.split_selfies(str_sf))
    tgt_tokens = list(sf.split_selfies(tgt_sf))
    sfw.write(cl+' '+' '.join(str_tokens)+'\n')
    tfw.write(' '.join(tgt_tokens) + '\n')
tfw.close()
sfw.close()


'''
def selfies_tokenizer(selfies):
    tokens = list(sf.split_selfies(selfies))
    return ' '.join(tokens)#,  sf.len_selfies(selfies)
    # return tokens
smiles='CC(C)CC(=O)c1ccc(O)nc1'
str_sf = sf.encoder(smiles)
print(str_sf)
print(selfies_tokenizer(str_sf))
# '''