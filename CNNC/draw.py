'''=================================================================================
This file explores the difference between label 1 group and label 0 group.
sys.argv[1]: "human1/human2/human3/mouse1/mouse2/mouse3"
====================================================================================
'''
from __future__ import print_function
import scipy.stats
import numpy as np
import h5py
import csv
import sys,os
import matplotlib.pyplot as plt
import pandas as pd
from pandas import Series,DataFrame
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
sns.set(font_scale=1.5)
import palettable
from math import*

def load_data_TF2(indel_list,data_path): # cell type specific  ## random samples for reactome is not enough, need borrow some from keggp
    import random
    import numpy as np
    xxdata_list = []
    yydata = []
    count_set = [0]
    count_setx = 0
    for i in indel_list:#len(h_tf_sc)):
        xdata = np.load(data_path+'/Nxdata_tf' + str(i) + '.npy')
        ydata = np.load(data_path+'/ydata_tf' + str(i) + '.npy')
        for k in range(len(ydata)):
            xxdata_list.append(xdata[k,:,:,:])
            yydata.append(ydata[k])
        count_setx = count_setx + len(ydata)
        count_set.append(count_setx)
        #print (i,len(ydata))
    yydata_array = np.array(yydata)
    yydata_x = yydata_array.astype('int')
    #print (np.array(xxdata_list).shape)
    return((np.array(xxdata_list),yydata_x,count_set))

csv_reader = csv.reader(open(sys.argv[1]+'/generatedata/index_train.csv','r'))
trainindex=np.array(list(csv_reader)).shape[0]
csv_reader = csv.reader(open(sys.argv[1]+'/generatedata/index_test.csv','r'))
testindex=np.array(list(csv_reader)).shape[0]

length_TF_train =int(trainindex-1)
length_TF_test =int(testindex-1)

data_path_train = sys.argv[1]+"/NEPDF_data_train"
data_path_test = sys.argv[1]+"/NEPDF_data_test"

num_classes = int(2)

train_TF = [i for i in range (length_TF_train)]
test_TF = [i for i in range (length_TF_test)]

(x_train, y_train,count_set_train) = load_data_TF2(train_TF,data_path_train)
(x_test, y_test,count_set) = load_data_TF2(test_TF,data_path_test)
print(x_train.shape, 'x_train samples')
print(x_test.shape, 'x_test samples')


data_train=np.mean(x_train,axis=0).squeeze()
data_test=np.mean(x_test,axis=0).squeeze()


length=x_train.shape[0]/2
size=int(sqrt(length)*(sqrt(length)-1))
train_1=np.zeros((32,32,1))
train_2=np.zeros((32,32,1))
for i in range(0, size*2,4):
    train_1=train_1+(x_train[i])
    train_1=train_1+(x_train[i+1])
    train_2=train_2+(x_train[i+2])
    train_2=train_2+(x_train[i+3])

train1_foo=train_1.squeeze()
train2_foo=train_2.squeeze()

x=np.arange(32)
y=np.arange(32)
x,y=np.meshgrid(x,y)
fig = plt.figure()
ax=Axes3D(fig)
ax.plot_surface(x, y, train2_foo, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_xlabel("gene a",fontsize=14, rotation=-15)
ax.set_ylabel("gene b",fontsize=14, rotation=60)
ax.set_zlabel("Coexpression",fontsize=14)
ax.zaxis.set_ticks_position('none')
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

ax.zaxis.set_tick_params(labelsize=10)
plt.savefig("label1.png",dpi=300)


