# README

## Environment
python3.6+

## Files

### Baseline
human1, human2, human3, mouse1, mouse2, mouse3

### rawdata
rawdata: 3 batches of human and mouse gene list known, we used human1 to process further modification.

### modified
fine-tune

simplified

### database
Two dictionaries for gene name and ID transform. Human scRNA expression data.


## Commands
`Please keep the structure of the files`

### Raw gene list to proper gene pairs 
Add the list needed to be transformed, according to its species.
```python
python3 human_data_process.py human1.csv
```

### Gene pairs to NEPDF data
Add the index, corresponding pair file and a number indicating whether it has labels.
```python
python3 NEPDF_human.py index_train.csv train.csv 1
```
 
### To use mouse data　& human data
Download its database to `database` folder.
```bash
https://s3.amazonaws.com/mousescexpression/rank_total_gene_rpkm.h5
```
Human database is in the folder database.

### To run svm.py 
```
python3 svm.py
```
 



# cs-433-project-2-refrain
# cs-433-project-2-refrain
# cs-433-project-2-refrain
