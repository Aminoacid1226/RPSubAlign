# Read me

The RPSubAlign method offers a novel molecular characterization for retrosynthesis by aligning reactant and product structures, thereby enhancing the comprehension of structural grammar by the model.

This repository contains all the data and source code for the RPSubAlign method. You can follow the steps below to implement RPSubAlign for inverse synthesis prediction.

Full data enhancement data and trained models can be found at [Google-drive](https://drive.google.com/drive/folders/1Mbx0m5HQhaZ2Wzv9jp35vxbiOpkb1Ahx?usp=drive_link). Or you can download our model directly from [hugging face](https://huggingface.co/Aminoacid1226/RPSA/tree/main). If you want to reproduce the results of the article, we recommend that you download the full data and model and organize the file as follows.

-raw data : raw data csv files for the USPO-50K and USPO-MIT datasets.

-script：

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--get_seq：get different molecular characterization scripts.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--Config_yml：retrosynthesis model preparation, training and use profile.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--metrics：result evaluation script.

-SMILES representation:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--data: data that can be used directly for model training and testing. Contains three different SMILES representation data for the USPTO-50K and USPTO-MIT datasets, respectively.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--model: The trained model, which can be directly used to reproduce the results of the article or to make inverse synthesis prediction. 

-SELFIES  representation: similar as SMILES representation.

# Environmental preparation

All code is written in Python3.6 and modeled using the pytorch version of openmnt-2.2.0. The following are the packages required for the method.

```python
·python 3.6
·opennmt-py 2.2.0
·rdkit 2020.09.1.0
·numpy 1.19.5
·pandas 1.1.5
·selfies 2.0.0
·torch 1.10.1
```

Or you can create the desired environment with one click using the yaml files provided.

```bash
# conda env create -f RPSubAlign.yaml -n your_env_name
conda env create -f RPSubAlign.yaml -n sa
```

# Substructure alignment

Substructure alignment of reactant and product data is required before retrosynthesis model training, which can be done by running the following script.

```bash
python3 get_RPSubAlign.py
```

You can change the data input path and processing result output path in the script.

# Retrosynthetic prediction

The model framework provided by opennmt is used to train and predict the retroynthesis.The configuration file can be found in `script/Config_yml`.

### step1：tokenization

The commands provided by opennmt can be used for quick word segmentation, such as:

```bash
onmt_build_vocab -config preprocess.yml
```

Running this command produces a vocabulary file in the specified folder.

### step2：model training

Similar to the previous step, you can train the model using the following command:

```bash
onmt_train -config train.yml
```

Running this command generates the chekpoint of the model in the specified folder.

### step3：model prediction

Similar to the previous step, you can use the following command to make predictions with the trained model, and don't forget to change the model path before using the command：

```bash
onmt_translate -config translate.yml
```

Running this command generates predictions in the specified folder.

# Metrics

We provide evaluation scripts for three indicators, Top-N, Maxfrag and Validity, which can be directly run to obtain results. You can find it in `script/metrcis` and run it.
