# Identification of the Core Chemical Structure in SureChEMBL Patents
This repository contains scripts for the extraction of Core Chemical Structures (CCS) representing the patent claim in SureChEMBL patents, and its use to filter out non-pharmacologically relevant compounds associated to the patent (e.g. common reactants, starting materials) while retaining molecules expemplyfing the patent claim.

# Overview
This repository contains scripts designed for the extraction of Core Chemical Structures (CCS) from SureChEMBL patents. These scripts are used to represent the key patent claims and filter out non-pharmacologically relevant compounds, such as common reactants or starting materials, while retaining molecules that exemplify the patent claims. More details can be found in [Falaguera and Mestres "Identification of the Core Chemical Structure in SureChEMBL Patents" (JCIM, 2021)](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00151). Please cite this article if you use this code.

# System Requirements

## Hardware requirements
The scripts require only a standard computer with enough RAM to support the in-memory operations.

## Software requirements

### OS Requirements
This package is supported for macOS. The package has been tested on the following systems:

- macOS: Sequoia (15.0.1)

### Python Dependencies
```
numpy # 1.26.4
rdkit # 2017.09.1
```

# Execution Guide

The following code takes a dictionary containing the SMILES strings of molecules associated with a patent. It processes the data to extract the Core Chemical Structures (CCS) from these molecules. The result is a dictionary where each CCS is a key, and the molecules containing the CCS are stored in the `filteredMolecules` field.

```
import ccs
idx_smiles = {'CHEMBL293613': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(C(=O)O)c3)cc21',
 'CHEMBL429540': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(CNCCCc4ccccc4)c3)cc21',
 'CHEMBL63861': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(CN4CCN(Cc5ccccc5)CC4)c3)cc21',
 'CHEMBL12465': 'NCc1ccn(-c2cc3c(cc2[N+](=O)[O-])[nH]c(=O)c(=O)n3CC(=O)O)c1',
 'CHEMBL67828': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(CN4CCN(CCc5ccccc5)CC4)c3)cc21',
 'CHEMBL416685': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(CN4CCC(c5ccccc5)CC4)c3)cc21',
 'CHEMBL268284': 'O=C(O)Cn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3cccc3)cc21'}
CCS = getCCSforMolecules(idx_smiles=idx_smiles)
```
