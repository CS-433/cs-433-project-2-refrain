import csv
import os
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import joblib
from sklearn import decomposition, svm, tree, model_selection

# performance visualization
def performance_visualization(labels, scores):
    p = Performance(labels, scores)
    print("confusion matrix")
    print(p.get_confusion_matrix())
    acc = p.accuracy()
    pre = p.precision()
    rec = p.recall()
    print('accuracy: %.2f' % acc)
    print('precision: %.2f' % pre)
    print('recall: %.2f' % rec)
    p.roc_plot()

# for classifier's performance evaluation
class Performance:

    def __init__(self, labels, scores, threshold=0.5):
        # true value of y
        self.labels = labels
        # predicted value of y
        self.scores = scores
        # classification threshold
        self.threshold = threshold
        self.db = self.get_db()
        self.TP, self.FP, self.FN, self.TN = self.get_confusion_matrix()

    def accuracy(self):
        return (self.TP + self.TN) / (self.TP + self.FN + self.FP + self.TN)

    def precision(self):
        return self.TP / (self.TP + self.FP)

    def recall(self):
        return self.TP / (self.TP + self.FN)

    def auc(self):
        auc = 0.
        prev_x = 0
        xy_arr = self.roc_coord()
        for x, y in xy_arr:
            if x != prev_x:
                auc += (x - prev_x) * y
                prev_x = x
        return auc

    # return coordinate of roc
    def roc_coord(self):
        xy_arr = []
        tp, fp = 0., 0.
        neg = self.TN + self.FP
        pos = self.TP + self.FN
        for i in range(len(self.db)):
            tp += self.db[i][0]
            fp += 1 - self.db[i][0]
            xy_arr.append([fp / neg, tp / pos])
        return xy_arr

    # plot roc
    def roc_plot(self):
        auc = self.auc()
        xy_arr = self.roc_coord()
        x = [_v[0] for _v in xy_arr]
        y = [_v[1] for _v in xy_arr]
        plt.title("ROC curve (AUC = %.4f)" % auc)
        plt.ylabel("True Positive Rate")
        plt.xlabel("False Positive Rate")
        plt.plot(x, y)
        plt.grid()
        plt.savefig("ROC.pdf")
        plt.show()

    def get_db(self):
        db = []
        for i in range(len(self.labels)):
            db.append([self.labels[i], self.scores[i]])
        db = sorted(db, key=lambda x: x[1], reverse=True)
        return db

    # calculate confusion matrix
    def get_confusion_matrix(self):
        tp, fp, fn, tn = 0., 0., 0., 0.
        for i in range(len(self.labels)):
            if self.labels[i] == 1 and self.scores[i] >= self.threshold:
                tp += 1
            elif self.labels[i] == 0 and self.scores[i] >= self.threshold:
                fp += 1
            elif self.labels[i] == 1 and self.scores[i] < self.threshold:
                fn += 1
            else:
                tn += 1
        return [tp, fp, fn, tn]


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

csv_reader = csv.reader(open('../human1/generatedata/index_train.csv','r'))
trainindex=np.array(list(csv_reader)).shape[0]
csv_reader = csv.reader(open('../human1/generatedata/index_test.csv','r'))
testindex=np.array(list(csv_reader)).shape[0]

length_TF_train =int(trainindex-1)
length_TF_test =int(testindex-1)

data_path_train = "../human1/NEPDF_data_train"
data_path_test = "../human1/NEPDF_data_test"

num_classes = int(2)

train_TF = [i for i in range (length_TF_train)]
test_TF = [i for i in range (length_TF_test)]

(xtrain, ytrain,count_set_train) = load_data_TF2(train_TF,data_path_train)
(xtest, ytest,count_set) = load_data_TF2(test_TF,data_path_test)
print(ytrain.shape, 'x_train samples')
print(ytest.shape, 'x_test samples')


def get_cluster(set_, path):
    store = h5py.File(path, "r")
    rpkm = store["counts"]
    cluster = []
    for item in set_:
        cluster.append(np.log10(np.array(rpkm[item]) + 10 ** -2))
    store.close()
    return np.array(cluster)

with open("../human1/generatedata/generatepair.csv", "r") as f:
    reader = csv.reader(f)
    whole_pair = [row for row in reader]

set_1 = set()
set_2 = set()
for item in whole_pair:
    set_1.add(item[0])
    set_2.add(item[1])
length=xtrain.shape[0]
xtrain=xtrain.reshape(length,1024)
length2=xtest.shape[0]
xtest=xtest.reshape(length2,1024)
cluster_1 = get_cluster(set_1, "../database/hey.hdf5")
cluster_2 = get_cluster(set_2, "../database/hey.hdf5")
#all_data = np.r_[cluster_1, cluster_2]
all_data=np.r_[xtrain, xtest]
#labels = np.array([1] * 31 + [0] * 31)
labels=np.r_[ytrain, ytest]
print (xtrain[:,0])
#print (all_data)
clf = svm.SVC()
#x_train, x_test, y_train, y_test = model_selection.train_test_split(all_data, labels, test_size=0.25, random_state=48)
clf.fit(xtrain, ytrain)

score_train=clf.score(xtrain, ytrain)
score_test=clf.score(xtest, ytest)
y_predict=clf.predict(xtest)
for i in range(y_predict.shape[0]):
    if y_predict[i]>0.5:
        y_predict[i]=1
    else:
        y_predict[i]=0 


print(score_train)
print(score_test)


labels = ytest
scores = y_predict
performance_visualization(labels, scores)


joblib.dump(clf, "svmmodel.m")
