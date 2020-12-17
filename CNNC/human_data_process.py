'''=================================================================================
This file transform raw human gene list to gene pairs.
sys.argv[1]: "human1.csv"
====================================================================================
'''
import pandas as pd
import scanpy as sc
import numpy as np
import csv
import json, re,os, sys
from matplotlib.image import NonUniformImage
import matplotlib.pyplot as plt
import h5py
import random

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
random.seed(2048)
save_dir = os.path.join(os.getcwd(),'generatedata')
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

## load gene expression lib
store=h5py.File("database/hey.hdf5","r+")

counts=store['counts']

keys=[key for key in counts.keys()]
keystmp=keys

## load known gene
with open("rawdata/"+sys.argv[1], newline='') as csvfile:
    reader = csv.reader(csvfile)
    rec = [row for row in reader]
    
for i in rec:
    i[0]=i[0].replace(u'\xa0', u'')


## remove gene lib does not contain from known gene

nokey=[]
for i in rec:
    try:
        keystmp.remove(i[0])
    except ValueError:
        print (i[0])
        nokey.append(i[0])
    else:
        print ("getkey")
      
for i in nokey:
    rec.remove([i])

length=len(rec)
print (length)

## generate known list that can be found in the library
with open(save_dir+"/rawgenelist.csv",'w') as f:
    csv_write = csv.writer(f)
    for i in range(length):
        csv_head=rec[i]
        csv_write.writerow(csv_head)  

## randomly select different gene from lib 
rand31 = random.sample(keystmp, length) 
with open(save_dir+"/randomlist.csv",'w') as f:
    csv_write = csv.writer(f)
    for i in range(length):
        csv_head=rand31[i]
        csv_write.writerow([csv_head])

## generate 1/3,2/3 gene pair
pair31=[]
for i in range(length):
	pair31.append([rec[i][0],rand31[i]])
with open(save_dir+"/generatepair.csv",'w') as f:
	csv_write= csv.writer(f)
	for i in range(length):
		csv_head=pair31[i]
		csv_write.writerow(csv_head)

print (len(pair31))
sl=int(length/4)

pair11 = pair31[0:sl]
pair20 = pair31[sl:length]

# with open(save_dir+"/testpair.csv",'w') as f:
# 	csv_write= csv.writer(f)
# 	for i in range(len(pair11)):
# 		csv_head=pair11[i]
# 		csv_write.writerow(csv_head)

# with open(save_dir+"/trainpair.csv",'w') as f:
# 	csv_write= csv.writer(f)
# 	for i in range(len(pair20)):
# 		csv_head=pair20[i]
# 		csv_write.writerow(csv_head)

##=============================================================================
# Generate train pair and test pair
#==============================================================================
# gene_known_train: a,b -> a*b*2

gene_known_train= [w[0] for w in pair20]
gene_unknown_train= [w[1] for w in pair20]

train1=[]
for i in range(length-sl):
    for j in range(i+1,length-sl):
        train1.append([gene_known_train[i],gene_known_train[j],'1'])
        train1.append([gene_known_train[j],gene_known_train[i],'1'])
for i in range(length-sl):
    train1.append([gene_known_train[i],gene_known_train[i],'1'])

train2=[]
for i in range(length-sl):
    for j in range(length-sl):
        train2.append([gene_known_train[i],gene_unknown_train[j],'0'])

train=[]
# (a,b) (b,a)
for i in range(0,(length-sl)*(length-sl-1),2): 
    train.append(train1[i])
    train.append(train1[i+1])
    train.append(train2[i])
    train.append(train2[i+1])
# (a,a)
for i in range(length-sl):
    train.append(train1[length-sl+i])
    train.append(train2[length-sl+i])


with open(save_dir+"/train.csv",'w') as f:
    csv_write = csv.writer(f)
    for i in train:
        csv_write.writerow(i)   

# gene_known_test: a,c -> a*c*2

gene_known_test= [w[0] for w in pair11]
gene_unknown_test= [w[1] for w in pair11]

test1=[]
for i in range(length-sl):
    for j in range(sl):
        test1.append([gene_known_train[i],gene_known_test[j],'1'])
test2=[]
for i in range(length-sl):
    for j in range(sl):
        test2.append([gene_known_train[i],gene_unknown_test[j],'0'])

test=test1+test2
random.shuffle(test)
with open(save_dir+"/test.csv",'w') as f:
    csv_write = csv.writer(f)
    for i in test:
        csv_write.writerow(i) 

whole = train+test
with open(save_dir+"/whole.csv",'w') as f:
    csv_write = csv.writer(f)
    for i in whole:
        csv_write.writerow(i) 

with open(save_dir+"/index_train.csv",'w') as f:
    csv_write = csv.writer(f)
    for i in range(0,len(train),4):
        csv_write.writerow([i])
    csv_write.writerow([(length-sl)**2*2])
with open(save_dir+"/index_test.csv",'w') as f:
    csv_write = csv.writer(f)
    for i in range(0,len(test),4):
        csv_write.writerow([i])     
    csv_write.writerow([(length-sl)*sl*2])
with open(save_dir+"/index_whole.csv",'w') as f:
    csv_write = csv.writer(f)
    for i in range(0,len(whole),4):
        csv_write.writerow([i])     
    csv_write.writerow([(length-sl)*length*2])


store.close()
