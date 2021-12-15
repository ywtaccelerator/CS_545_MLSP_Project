import numpy as np
import matplotlib.pyplot as plt
import h5py
import gc
import matplotlib
matplotlib.use('Agg')
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM
from tensorflow.keras.layers import Layer
from tensorflow.keras.layers import Bidirectional
import tensorflow.keras.backend as K
import tensorflow as tf
from tensorflow.keras.layers import Dense

import math
import matplotlib.pyplot as plt


num_features_after_pca = 60
my_num_epochs = 2000
my_batch_size = 2000


#read data
data_set_folder = '../../data/CS545_data_set_PCA_60_20211213_intact/'
x_path = data_set_folder + 'x_matrix.mat'
x_matrix_mat = h5py.File(x_path, mode='r')
x_matrix = np.transpose(x_matrix_mat['x_matrix'][:].astype(np.float32))
num_features, num_samples_of_each_subject = x_matrix.shape

TN_subjectNum = 22
subject_num_list_for_test = [21, 22] #change every time

# construct the training set
are_training_sets_not_initialized = True
for subject_num in range(1, TN_subjectNum + 1):
    if subject_num in subject_num_list_for_test:
            continue
    subject_num_str = str(subject_num + 100000)[-2:]
    y_path = data_set_folder + 'sj_' + subject_num_str + '_Y.mat'
    y_list_mat = h5py.File(y_path, mode='r')
    y_list = np.transpose(y_list_mat['y_list'][:].astype(np.float32))

    # balance the samples of each subject
    count0 = 0
    count1 = 0
    for i in range(len(y_list)):
        if y_list[i][0] == 0:
            count0 = count0 + 1
        else:
            count1 = count1 + 1

    x_matrix_y0 = np.zeros((num_features, count0))
    x_matrix_y1 = np.zeros((num_features, count1))
    x_matrix_y0_pointer = 0
    x_matrix_y1_pointer = 0

    for column_num in range(num_samples_of_each_subject):
        x_sample = x_matrix[:, column_num]
        if y_list[column_num][0] == 0:
            x_matrix_y0[:, x_matrix_y0_pointer] = x_sample
            x_matrix_y0_pointer = x_matrix_y0_pointer + 1
        else:
            x_matrix_y1[:, x_matrix_y1_pointer] = x_sample
            x_matrix_y1_pointer = x_matrix_y1_pointer + 1

    # shuffle all columns of x_matrix_y0
    perm = np.arange(count0)
    np.random.shuffle(perm)
    x_matrix_y0 = x_matrix_y0.T[perm]
    x_matrix_y0 = x_matrix_y0.T

    # shuffle all columns of x_matrix_y1
    perm = np.arange(count1)
    np.random.shuffle(perm)
    x_matrix_y1 = x_matrix_y1.T[perm]
    x_matrix_y1 = x_matrix_y1.T

    if count0 < count1:
        x_matrix_y1 = x_matrix_y1[:, :count0]
    else:
        x_matrix_y0 = x_matrix_y0[:, :count1]

    if are_training_sets_not_initialized:
        training_set_y0 = x_matrix_y0
        training_set_y1 = x_matrix_y1
        are_training_sets_not_initialized = False
    else:
        training_set_y0 = np.hstack((training_set_y0, x_matrix_y0))
        training_set_y1 = np.hstack((training_set_y1, x_matrix_y1))


y0_feature_list = []
for sample_n in range(len(training_set_y0[0])):
    features = np.reshape(training_set_y0[:, sample_n], (num_features_after_pca, -1), 'F')
    y0_feature_list.append(features)
del training_set_y0
gc.collect()

y1_feature_list = []
for sample_n in range(len(training_set_y1[0])):
    features = np.reshape(training_set_y1[:, sample_n], (num_features_after_pca, -1), 'F')
    y1_feature_list.append(features)
del training_set_y1
gc.collect()


time_steps = len(y0_feature_list[0][0]) #300
num_features = len(y0_feature_list[0]) #320

#input_array:total_sample * time_steps * features
#nead to transpose every sample.
x_input = []
y_input = []
for y0_feature in y0_feature_list:
    x_input.append(np.transpose(y0_feature).tolist())
    y_input.append(0)
del y0_feature_list
gc.collect()

for y1_feature in y1_feature_list:
    x_input.append(np.transpose(y1_feature).tolist())
    y_input.append(1)
del y1_feature_list
gc.collect()


# shuffle all X and Y with the same order
seed_num = np.random.randint(0, 99999)
np.random.seed(seed_num)
np.random.shuffle(x_input)
np.random.seed(seed_num)
np.random.shuffle(y_input)


#test data
#x_matrix
#reshape feature, move right part of the matrix to under left
x_feature_list = []

for sample_n in range(len(x_matrix[0])): #len(x_matrix[0])
    features = np.reshape(x_matrix[:, sample_n], (num_features_after_pca, -1), 'F')
    x_feature_list.append(features)
test_input = []
for x_feature in x_feature_list:
    test_input.append(np.transpose(x_feature).tolist())

del x_matrix
gc.collect()

test_validation_x = test_input
test_validation_x.extend(test_input)


test_validation_y = []
# subject # 21
subject_num_str = str(subject_num_list_for_test[0] + 100000)[-2:]
y_path = data_set_folder + 'sj_' + subject_num_str + '_Y.mat'
y_list_mat = h5py.File(y_path, mode='r')
y_list = np.transpose(y_list_mat['y_list'][:].astype(np.float32))

for i in range(len(y_list)):
    if np.int(y_list[i][0]) == 0:
        test_validation_y.append(0)
    elif np.int(y_list[i][0]) == 1:
        test_validation_y.append(1)

# subject # 22
subject_num_str = str(subject_num_list_for_test[1] + 100000)[-2:]
y_path = data_set_folder + 'sj_' + subject_num_str + '_Y.mat'
y_list_mat = h5py.File(y_path, mode='r')
y_list = np.transpose(y_list_mat['y_list'][:].astype(np.float32))

for i in range(len(y_list)):
    if np.int(y_list[i][0]) == 0:
        test_validation_y.append(0)
    elif np.int(y_list[i][0]) == 1:
        test_validation_y.append(1)

# Visualization
def display_loss_curve(loss, val_loss):
    plt.clf()
    plt.rcParams['figure.dpi'] = 600
    plt.rcParams['savefig.dpi'] = 600
    plt.plot([i + 1 for i in range(len(loss))], loss, linewidth=1, alpha=0.8)
    plt.plot([i + 1 for i in range(len(val_loss))], val_loss, linewidth=1, alpha=0.8)
    # Setup the legend
    plt.legend(["Loss", "Validation loss"])
    # Labels for the x and y axis
    plt.ylabel("Value")
    plt.xlabel("Epoch")
    plt.title("Loss curve")
    plt.grid()
    # plt.show()
    plt.savefig('loss_curve.png',
                format='png')

def display_accuracy_curve(accuracy, val_accuracy):
    plt.clf()
    plt.rcParams['figure.dpi'] = 600
    plt.rcParams['savefig.dpi'] = 600
    plt.plot([i + 1 for i in range(len(accuracy))], accuracy, linewidth=1, alpha=0.8)
    plt.plot([i + 1 for i in range(len(val_accuracy))], val_accuracy, linewidth=1, alpha=0.8)
    # Setup the legend
    plt.legend(["Accuracy", "Validation accuracy"])
    # Labels for the x and y axis
    plt.ylabel("Value")
    plt.xlabel("Epoch")
    plt.title("Accuracy curve")
    plt.grid()
    # plt.show()
    plt.savefig('accuracy_curve.png',
                format='png')


#Add attention
#https://machinelearningmastery.com/adding-a-custom-attention-layer-to-recurrent-neural-network-in-keras/

class attention(Layer):
    def __init__(self,**kwargs):
        super(attention,self).__init__(**kwargs)

    def build(self,input_shape):
        self.W=self.add_weight(name='attention_weight', shape=(input_shape[-1],1),
                               initializer='random_normal', trainable=True)
        self.b=self.add_weight(name='attention_bias', shape=(input_shape[1],1),
                               initializer='zeros', trainable=True)
        super(attention, self).build(input_shape)

    def call(self,x):
        # Alignment scores. Pass them through tanh function
        e = K.tanh(K.dot(x,self.W)+self.b)
        # Remove dimension of size 1
        e = K.squeeze(e, axis=-1)
        # Compute the weights
        alpha = K.softmax(e)
        # Reshape to tensorFlow format
        alpha = K.expand_dims(alpha, axis=-1)
        # Compute the context vector
        context = x * alpha
        context = K.sum(context, axis=1)
        return context

#LSTM
from tensorflow.keras.layers import RNN, LSTMCell
model = Sequential()
two_lstm = RNN([LSTMCell(80), LSTMCell(80), LSTMCell(80)], input_shape=(time_steps, num_features),  return_sequences=True)
model.add(Bidirectional(two_lstm, input_shape=(time_steps, num_features)))
model.add(attention())
model.add(Dense(1, activation="sigmoid"))
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

history = model.fit(x_input, y_input, validation_data=(test_validation_x, test_validation_y), epochs=my_num_epochs, batch_size=my_batch_size, verbose=2)
model.summary()
model.save("trained_LSTM_model_with_multiplayer.h5")



# Prediction
print("\n\nFitting finished!\n\n")
print("model.evaluate(test_validation_x, test_validation_y):\n")
print(model.evaluate(test_validation_x, test_validation_y))
print('\n\n')


from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score as sklearn_f1_score
from sklearn.metrics import confusion_matrix, classification_report

predicted_y = model.predict(test_validation_x)
model_accuracy = accuracy_score(test_validation_y,
                                np.around(predicted_y))
macro_f1_score = sklearn_f1_score(test_validation_y,
                                np.around(predicted_y),
                                average='macro')
micro_f1_score = sklearn_f1_score(test_validation_y,
                                np.around(predicted_y),
                                average='micro')


# write out the macro_f1_score
macro_f1_score_txtFile = open('macro_f1_score.txt', 'w')
macro_f1_score_txtFile.write(str(macro_f1_score))
macro_f1_score_txtFile.close()

# write out the model_accuracy
model_accuracy_txtFile = open('model_accuracy.txt', 'w')
model_accuracy_txtFile.write(str(model_accuracy))
model_accuracy_txtFile.close()

# write out the micro_f1_score
micro_f1_score_txtFile = open('micro_f1_score.txt', 'w')
micro_f1_score_txtFile.write(str(micro_f1_score))
micro_f1_score_txtFile.close()







lstm_output = model.predict(test_input)
# print("========================output======================")
# print(lstm_output)

# res=0
# for o in lstm_output:
#     res += o
# print("average of output:")
# print(res / len(lstm_output))


transformed_lstm_output = []
ones = 0
zeros = 0
for i in range(len(lstm_output)):
    if lstm_output[i] > 0.5:
        transformed_lstm_output.append(1)
        ones+=1
    else:
        transformed_lstm_output.append(0)
        zeros+=1

# print("===============# 1s: ================")
# print(ones)

# print("===============#  0s: ================")
# print(zeros)

#test y21
subject_num_str = str(subject_num_list_for_test[0] + 100000)[-2:]
y_path = data_set_folder + 'sj_' + subject_num_str + '_Y.mat'
y_list_mat = h5py.File(y_path, mode='r')
y_list = np.transpose(y_list_mat['y_list'][:].astype(np.float32))

match = 0
match_ones = 0
match_zeros = 0
for i in range(len(y_list)):
    if y_list[i][0] == transformed_lstm_output[i]:
        match += 1
        if y_list[i][0] == 0:
                match_zeros += 1
        else:
                match_ones += 1

print("match:")
print(match)
print("match ones:")
print(match_ones)
print("match zeros:")
print(match_zeros)

accuracy = match / len(y_list)
print("=====================accuracy y21: =====================")
print(accuracy)
acc_21 = accuracy

#test y22
subject_num_str = str(subject_num_list_for_test[1] + 100000)[-2:]
y_path = data_set_folder + 'sj_' + subject_num_str + '_Y.mat'
y_list_mat = h5py.File(y_path, mode='r')
y_list = np.transpose(y_list_mat['y_list'][:].astype(np.float32))

match = 0
match_ones = 0
match_zeros = 0
for i in range(len(y_list)):
    if y_list[i][0] == transformed_lstm_output[i]:
        match += 1
        if y_list[i][0] == 0:
                match_zeros += 1
        else:
                match_ones += 1

print("match:")
print(match)
print("match ones:")
print(match_ones)
print("match zeros:")
print(match_zeros)

accuracy = match / len(y_list)
print("=====================accuracy y22: =====================")
print(accuracy)
acc_22 = accuracy


display_loss_curve(history.history['loss'], history.history['val_loss'])
display_accuracy_curve(history.history['accuracy'], history.history['val_accuracy'])

np.savez('predictions_on_intact.npz', prediction_on_intact = np.array(transformed_lstm_output))

# write out the model_accuracy
model_accuracy_matching = (acc_21 + acc_22) / 2
model_accuracy_txtFile = open('model_accuracy_matching.txt', 'w')
model_accuracy_txtFile.write(str(model_accuracy_matching))
model_accuracy_txtFile.close()

np.savez('training_history.npz', history = history)


# write out the confusion matrices
txtFile = open('Confusion matrices.txt', 'w')
txtFile.write(
    '\n\nTraining confusion matrix and classification report:' + '\n')
training_pred = model.predict(x_input)
txtFile.write(str(confusion_matrix(y_input, np.around(training_pred))) + '\n\n')
txtFile.write(str(classification_report(y_input, np.around(training_pred))) + '\n\n')

txtFile.write(
    '\nTest confusion matrix and classification report:' + '\n')
txtFile.write(str(confusion_matrix(test_validation_y, np.around(predicted_y))) + '\n\n')
txtFile.write(str(classification_report(test_validation_y, np.around(predicted_y))) + '\n\n')

txtFile.close()

