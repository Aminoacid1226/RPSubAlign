import os
import selfies as sf


folder_path = 'pred/smiles'


for filename in os.listdir(folder_path):
    if filename.endswith(".txt"):
       
        file_path = os.path.join(folder_path, filename)
        
    
        with open(file_path, 'r') as sfile:
            lines = sfile.readlines()
        
  
        pred_sf = [''.join(line.strip().split()) for line in lines]
        

        output_file_path = os.path.join(folder_path, filename.replace('.txt', '_smi.txt'))
        
        with open(output_file_path, 'w') as sfw:
            for i in range(len(pred_sf)):
                pred_smi = sf.decoder(pred_sf[i])
                sfw.write(pred_smi + '\n')


print("Done!")
