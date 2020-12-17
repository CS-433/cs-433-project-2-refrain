'''=================================================================================
This file transform gene pairs to NEPDF matrix which will be fed to CNN.
sys.argv[1]: "index_train.csv"
sys.argv[2]: "train.csv"
sys.argv[3]: '1'
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

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
datapath="generatedata/"
num_class=sys.argv[3]

if sys.argv[1]=='index_train.csv':
    save_dir = os.path.join(os.getcwd(),'NEPDF_data_train')
elif sys.argv[1]=='index_test.csv':
    save_dir = os.path.join(os.getcwd(),'NEPDF_data_test')
else:
    save_dir = os.path.join(os.getcwd(),'NEPDF_data_whole')

if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

def get_sepration_index(file_name):
    import numpy as np
    index_list = []
    s = open(file_name, 'r')
    for line in s:
        index_list.append(int(line))
    return (np.array(index_list))

def get_sc_list(file_name):# name to number
	with open(file_name, newline='') as csvfile:
		reader = csv.reader(csvfile)
		ensembl = [row for row in reader]
	names=np.array(ensembl)[:,0]
	length=np.array(ensembl)[:,1]
	dic = dict(zip(names,length))
	return dic


h_gene_list =get_sc_list("database/name-number.csv")

store = pd.HDFStore("database/rank_total_gene_rpkm.h5")
counts = store['rpkm']
store.close()
#store_bulk = pd.HDFStore("data/mouse_bulk.h5")
#counts_bulk = store_bulk['rpkm']
#store_bulk.close()
print('read lib')

gene_pair_label = []
gene_pair_index = get_sepration_index(datapath+sys.argv[1])

with open(datapath+sys.argv[2],'r') as f: #train/test/whole -> corresponding index changes
    reader = csv.reader(f)
    gene_pair_label = [row for row in reader]

for i in range(len(gene_pair_index)-1):   #### many sperations
    print (i)
    start_index = gene_pair_index[i]
    end_index = gene_pair_index[i+1]
    x = []
    y = []
    z = []
    for gene_pair in gene_pair_label[start_index:end_index]: ## each speration
        if num_class == '1': 
            x_gene_name,y_gene_name,label = gene_pair[0],gene_pair[1],gene_pair[2]
            y.append(label)
        else:
            x_gene_name, y_gene_name = gene_pair[0], gene_pair[1]
        z.append(x_gene_name+'\t'+y_gene_name)
 
        x_gene = np.log10(counts[int(x_gene_name)]*10 + 10 ** -4) 
        y_gene = np.log10(counts[int(y_gene_name)]*10 + 10 ** -4) 
        H_T = np.histogram2d(x_gene, y_gene, bins = 32)
        H = H_T[0].T
        HT = (np.log10(H / x_gene.shape[0] + 10 ** -4) + 4) / 4       

        x.append(HT)

    if (len(x)>0):
        xx = np.array(x)[:, :, :, np.newaxis]
    else:
        xx=np.array(x)
    np.save(save_dir+'/Nxdata_tf' + str(i) + '.npy', xx)
    if num_class=='1':
        np.save(save_dir+'/ydata_tf' + str(i) + '.npy', np.array(y))
    np.save(save_dir+'/zdata_tf' + str(i) + '.npy', np.array(z))

