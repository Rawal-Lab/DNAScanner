# DNAScanner
DNAScanner is a tool wirtten in Python that scans genomic DNA for a number of different physicochemical properties by incorporating biophysical, thermodynamic, protein interactions and sequence based features.

### Requirements

Miniconda/Anaconda

Git

### Installation
Clone the repo

```bash
git clone https://github.com/Rawal-Lab/DNAScanner.git
cd DNAScanner/
```

Set up environment 

```bash
conda env create -f DNAScanner_env.yml
```
Done.

### Procedure 
Run via terminal
Syntax : 

```bash
python3 DNAScanner.py 'input-filename'.fasta
```

### Output
Nucleotide_Concentration : These outputs comprise of csv files mentioning the concentration of the nucleotide group in the given block.

Parameters : These outputs comprise of csv files mentioning the parameters of the nucleotide group wrt. their positions in the sequence.

Plots : HTML generated plots with sliders.  
