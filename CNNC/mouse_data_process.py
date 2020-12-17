'''=================================================================================
This file transform raw gene list to gene pairs.
sys.argv[1]: "mouse1.csv"
====================================================================================
'''
import pandas as pd
import numpy as np
import csv
import json, re,os, sys
from matplotlib.image import NonUniformImage
import matplotlib.pyplot as plt
import h5py
import random
import math

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
random.seed(2048)
save_dir = os.path.join(os.getcwd(),'generatedata')
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

def get_gene_list(file_name):
	with open(file_name, newline='') as csvfile:
		reader = csv.reader(csvfile)
		ensembl = [row for row in reader]
	names=np.array(ensembl)[:,1]
	length=np.array(ensembl)[:,0]
	dic = dict(zip(names,length))
	return dic


def get_sc_list(file_name):
	with open(file_name, newline='') as csvfile:
		reader = csv.reader(csvfile)
		ensembl = [row for row in reader]
	names=np.array(ensembl)[:,0]
	length=np.array(ensembl)[:,1]
	dic = dict(zip(names,length))
	return dic


store=pd.HDFStore("database/rank_total_gene_rpkm.h5")
counts=store['rpkm']
store.close()

keys=[key for key in counts.keys()]
keystmp=keys

with open("rawdata/"+sys.argv[1], newline='') as csvfile:
    reader = csv.reader(csvfile)
    rec0 = [row for row in reader]
    
for i in rec0:
    i[0]=i[0].replace(u'\xa0', u'')

dic =get_gene_list("name-ensembl.csv")
scdic=get_sc_list("name-number.csv")

## remove gene that cannot convert to symbol name & ID and that lib does not contain

nokey=[]
s=0
for i in rec0:
    try:
        counts[int(scdic[dic[i[0]]])]
    except KeyError:
        print (i)
        nokey.append(i)
        s=s+1
    else:
        print (i,dic[i[0]])
print (s)
for i in nokey:
    rec0.remove(i)


length=len(rec0)
print (length)
## randomly select different gene from lib
rand31 = random.sample(keystmp, length) 
        
        
pair31 = []
for i in range(length):
    pair31.append([int(scdic[dic[rec0[i][0]]]),rand31[i]])




##=============================================================================
# Generate train pair and test pair
#==============================================================================

length=len(pair31)
print (length)

sl=math.floor(length/4)
pair11 = pair31[0:sl]
pair20 = pair31[sl:length]


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

# gene_known_test: 20,11 -> 20*11*2

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



