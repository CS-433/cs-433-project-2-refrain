import numpy as np
import keras
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping,ModelCheckpoint

def get_data(x_path, y_path):
    x = np.load(x_path)
    y = np.load(y_path)
    y = y.astype(np.float)
    return x, y

'''x_train, y_train = get_data("NEPDF_data/train_Nxdata_tf.npy", "NEPDF_data/train_ydata_tf.npy")
x_test, y_test = get_data("NEPDF_data/test_Nxdata_tf.npy", "NEPDF_data/test_ydata_tf.npy")'''
x_train, y_train = get_data("NEPDF_data/Nxdata_tf.npy", "NEPDF_data/ydata_tf.npy")

model = Sequential()
model.add(Conv2D(32, (3, 3), padding='same',input_shape=x_train.shape[1:]))
model.add(Activation('relu'))
model.add(Conv2D(32, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Conv2D(64, (3, 3), padding='same'))
model.add(Activation('relu'))
model.add(Conv2D(64, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Conv2D(128, (3, 3), padding='same'))
model.add(Activation('relu'))
model.add(Conv2D(128, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Flatten())
model.add(Dense(512))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(3))
model.add(Activation('softmax'))
model.load_weights("trained_models/KEGG_keras_cnn_trained_model_shallow.h5")

new_model = Sequential()
for layer in model.layers[:-2]:
    #layer.trainable = False                       # freeze param or not
    new_model.add(layer)

new_model.add(Dense(1))
new_model.add(Activation("sigmoid"))
sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
new_model.compile(optimizer=sgd,loss='binary_crossentropy',metrics=['accuracy'])

'''new_model.add(Dense(2))
new_model.add(Activation("softmax"))
sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
new_model.compile(optimizer=sgd,loss='hinge',metrics=['accuracy'])'''

early_stop = EarlyStopping(monitor="val_acc", patience=150, verbose=0, mode="auto")
checkpoint = ModelCheckpoint(filepath="", monitor="val_acc", verbose=1, save_best_only=True, mode="auto", period=1)
#history = new_model.fit(x=x_train, y=y_train, batch_size=32, epochs=200, validation_data=(x_test, y_test), shuffle=True)#, callbacks=[checkpoint, early_stop])
history = new_model.fit(x=x_train, y=y_train, batch_size=32, epochs=200, validation_split=0.2, shuffle=True)#, callbacks=[checkpoint, early_stop])
print(history.history["val_acc"])