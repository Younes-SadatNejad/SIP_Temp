# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 18:45:05 2021.

@author: mahdi
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestCentroid
import statistics
import math
from scipy import stats
from scipy.stats import linregress
import pandas as pd
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.metrics import hinge_loss

# %% Functions


def unit_vector(vector):
    """
    Compute the unit vector.

    Parameters
    ----------
    vector : numpy array
        The input vector.

    Returns
    -------
    TYPE : numpy array
        The unit vector of the input.

    """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """
    Calculate the angle between two vectors.

    Parameters
    ----------
    v1 : numpy array
        vector 1.
    v2 : numpu array
        vector 2.

    Returns
    -------
    TYPE :
        The angle between two vectors in raidan.

    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def projection_on_line(c_center_1, c_center_2, original_data):
    """
    Calculate the projection of one data points on the line going through \
        bothcluster centers.

    Parameters
    ----------
    c_center_1 : numpy 1 by 2 array
        first center coordinates.
    c_center_2 : numpy 1 by 2 array
        scond center coordinates.
    original_data : numpy n by 2 array
        data points.

    Returns
    -------
    projection : numpy array
        the coordinates of the points projected on to the line going through\
            the line which connects the two centers.

    """
    vector_data = original_data - c_center_1
    projection_line = c_center_1 - c_center_2
    projection = c_center_1 + np.dot(vector_data, projection_line) /\
        np.dot(projection_line, projection_line) * projection_line
    return projection


def calculate_center(original_data):
    """
    Calculate the center of data points for the label.

    Parameters
    ----------
    original_data : numpy array
        The data points.

    Returns
    -------
    center_co : numpy array
        The coordinates of the center point.

    """
    avr_vec = np.sum(original_data, axis=0)
    center_co = avr_vec/original_data.shape[0]
    return center_co


def calculate_pvar(pdata):
    """
    Calculate the variance of the data projected on to the line.

    Parameters
    ----------
    pdata : numpy array
        the coordinates of the data projected on the line

    Returns
    -------
    data_var : numpy array
        the variance of the projected data points on the line.

    """
    c_center = calculate_center(pdata)
    mean_vec = np.full(pdata.shape, c_center)
    temp_disvec = pdata - mean_vec
    temp_vec = []
    for i in range(pdata.shape[0]):
        sign_v = np.dot(unit_vector(temp_disvec[1, :]),
                        unit_vector(temp_disvec[i, :]))
        temp_valu = np.sign(sign_v) * np.linalg.norm(temp_disvec[i, :])
        temp_vec.append(temp_valu)
    # temp_vec = np.linalg.norm(temp_disvec, axis=1)
    temp_vec = np.array(temp_vec)
    data_var = np.var(temp_vec)
    return data_var


def calculate_dvar(pdata):
    """
    Calculate the variance of the data based on the distance from central\
        point.

    Parameters
    ----------
    pdata : numpy array
        the coordinates of the data projected on the line

    Returns
    -------
    data_var : numpy array
        the variance of the projected data points on the line.

    """
    c_center = calculate_center(pdata)
    mean_vec = np.full(pdata.shape, c_center)
    temp_disvec = pdata - mean_vec
    temp_vec = np.linalg.norm(temp_disvec, axis=1)
    temp_pvec = np.power(temp_vec, 2)
    temp_sum = np.sum(temp_pvec)
    data_var = temp_sum / pdata.shape[0]
    return data_var


def rotate_data(X_data, y):
    """
    Do the rotation to make variance calculation easier.

    Parameters
    ----------
    X_data : numpy array
        The data points that we want to rotata.
    y : numpy array
        Labels for X_data.

    Returns
    -------
    X_rotated : numpy array
        Rotated numpy array.

    """
    X_datap = X_data[y == 1]
    X_datan = X_data[y == -1]
    center_p = calculate_center(X_datap)
    center_n = calculate_center(X_datan)
    slope = (center_p[1] - center_n[1])/(center_p[0] - center_n[0])
    # slope = (X_data[0, 1] - X_data[1, 1])/(X_data[0, 0] - X_data[1, 0])
    angle = (math.atan(slope))
    theta = -angle
    c, s = np.cos(theta), np.sin(theta)
    rotation_mat = np.array(((c, -s), (s, c)))
    X_rotated = []
    for i in range(X_data.shape[0]):
        X_rot = rotation_mat.dot(X_data[i])
        X_rotated.append(X_rot)
    X_rotated = np.array(X_rotated)
    return X_rotated


# %% Generating the data
n_samples_1 = 2000
n_samples_2 = 2000
centers = [[-2, 0.0], [2, 2.0]]  # cluster centers
clusters_std = [0.7, 0.7]  # cluster std_dev
X, y = make_blobs(n_samples=[n_samples_1, n_samples_2],
                  centers=centers,
                  cluster_std=clusters_std,
                  random_state=0, shuffle=False)

y = np.where(y == 1, 1, -1)


# %% Preprocessing step
scaler = StandardScaler()
# X_s = scaler.fit_transform(X)
X_s = X
X_pos = X_s[y == 1]
X_neg = X_s[y == -1]
center_1 = NearestCentroid()
center_1.fit(X_s, y)
data_centers = center_1.centroids_
c_y = np.array([[1], [-1]])
pos_center = calculate_center(X_pos)
neg_center = calculate_center(X_neg)
print(f'The cluster centers are: {center_1.centroids_}')

# %% calculating S&S for clusters

# Calulate the distance of the centers
distance = np.linalg.norm(data_centers[0, :] - data_centers[1, :])

# First projecting the data on to the line which go through the cetners
X_pro = []
for i in range(X_s.shape[0]):
    projected_data = projection_on_line(data_centers[0, :], data_centers[1, :],
                                        X_s[i])
    X_pro.append(projected_data)

X_pro = np.array(X_pro)

X_pro_pos = X_pro[y == 1]
X_pro_neg = X_pro[y == -1]
var_x_pos = calculate_pvar(X_pro_pos)
var_x_neg = calculate_pvar(X_pro_neg)
total_var = ((X_pro_pos.shape[0] * var_x_pos) +
             (X_pro_neg.shape[0] * var_x_neg)) / (X_pro_pos.shape[0] +
                                                  X_pro_neg.shape[0])
sigma = np.sqrt(total_var)
SandS = 20 * np.log10(distance / (6 * sigma))
# Projection of the data on to the X axis
X_rota = rotate_data(X_pro, y)
X_rota_pos = X_rota[y == 1]
X_rota_neg = X_rota[y == -1]

# %% Plotting the data and centeral points
fig, ax = plt.subplots()
ax.scatter(X_s[:, 0], X_s[:, 1], marker="o", s=20,
           color=["coral" if y == -1 else "cyan" for y in y])
ax.scatter(data_centers[:, 0], data_centers[:, 1],
           color=["lime" if y == 1 else "r" for y in c_y])

# %% plotting the projection on to the line going throught two centers
fig, ax = plt.subplots()
# xmin, xmax = -10, 10
# ax.set_xlim([xmin, xmax])
# ax.set_ylim([xmin, xmax])

# Move left y-axis and bottim x-axis to centre, passing through (0,0)
# ax.spines['left'].set_position('zero')
# ax.spines['bottom'].set_position('zero')
# Eliminate upper and right axes
# ax.spines['right'].set_color('none')
# ax.spines['top'].set_color('none')

# Show ticks in the left and lower axes only
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# make the box square shape
ax.set_aspect('equal')

ax.scatter(X_pro[:, 0], X_pro[:, 1], marker="o", s=20,
           color=["r" if y == -1 else "b" for y in y],  alpha=0.5)
ax.scatter(X_s[:, 0], X_s[:, 1], alpha=0.5)
start, end = ax.get_xlim()
ax.xaxis.set_ticks(np.arange(start, end, 3.0))
ax.set_title('Projected and datas')

# %% Plotting the rotated data

fig, ax = plt.subplots()
# xmin, xmax = -5, 0
# ax.set_xlim([xmin, xmax])
# ax.set_ylim([xmin, xmax])

# Move left y-axis and bottim x-axis to centre, passing through (0,0)
# ax.spines['left'].set_position('zero')`
# ax.spines['bottom'].set_position('zero')
# Eliminate upper and right axes
# ax.spines['right'].set_color('none')
# ax.spines['top'].set_color('none')

# Show ticks in the left and lower axes only
# ax.xaxis.set_ticks_position('bottom')
# ax.yaxis.set_ticks_position('left')

# make the box square shape
# ax.set_aspect('equal')

ax.scatter(X_rota[:, 0], X_rota[:, 1], marker="o", s=20,
           color=["r" if y == -1 else "b" for y in y])
start, end = ax.get_xlim()
ax.xaxis.set_ticks(np.arange(start, end, 3.0))
# %% Ishtiaque approch
# make a dataframe with following columns
cols = ['iteration', 'C', 'Margin', 'Train_hinge_loss', 'cost_training',
        'Test_hinge_loss', 'cost_testing']
lst = []
iteration_num = 10
for i in range(1, iteration_num):
    X_train, X_test, y_train, y_test = train_test_split(X_s, y, test_size=0.40,
                                                        random_state=1)
    i = i
    Cs = np.logspace(-1, 2, 1000).tolist()
    Cs = np.array(Cs)
    clf = svm.SVC(kernel='linear', C=Cs)
    C = []
    Margin = []
    train_errors = []
    test_errors = []
    number_of_misclassified_train_points = []
    number_of_misclassified_test_points = []
    Train_hinge_loss = []
    cost_training = []
    Test_hinge_loss = []
    cost_testing = []

    for C in Cs:
        clf.set_params(C=C)
        clf.fit(X_train, y_train)
        i = i
        w = clf.coef_[0]
        y_train_predict = clf.predict(X_train)
        train_error = metrics.mean_squared_error(y_train, y_train_predict)
        train_errors.append(train_error)

        misclassified_train = np.where(y_train != y_train_predict)
        number_of_misclassified_train_points.append(misclassified_train)

        pred_decision_train = clf.decision_function(X_train)
        hinge_loss_train = hinge_loss(y_train, pred_decision_train)
        Train_hinge_loss.append(hinge_loss_train)

        pred_decision_test = clf.decision_function(X_test)
        hinge_loss_test = hinge_loss(y_test, pred_decision_test)
        Test_hinge_loss.append(hinge_loss_test)

        cost_train = 1/2 * np.dot(w, w) + C * hinge_loss_train
        cost_training.append(cost_train)

        cost_test = 1/2 * np.dot(w, w) + C * hinge_loss_test
        cost_testing.append(cost_test)

        # alpha=clf.dual_coef_
        # alphas.append(alpha)
        # ξ=y_train*clf.decision_function(X_train)
        # ξs.append(ξ)
        a = -w[0] / w[1]
        M = 2 / np.sqrt(np.sum(w ** 2))
        Margin.append(M)

        lst.append([i, C, M, hinge_loss_train, cost_train, hinge_loss_test,
                    cost_test])

comp_list = []
df = pd.DataFrame(lst, columns=cols)

for i in range(iteration_num):
    temp_df = df[df['iteration'] == i]
    temp_ar = temp_df.to_numpy()
    comp_list.append(temp_ar)

del comp_list[0]

array_sum = comp_list[0] + comp_list[1]
for i in range(len(comp_list)-2):
    array_sum = array_sum + comp_list[i+2]

averaged_data = array_sum/len(comp_list)

# plotting the average
fig, ax = plt.subplots()
ax.plot(averaged_data[:, 2], averaged_data[:, 5])

ax.set(xlabel='C values', ylabel='test cost',
       title='test')
ax.grid()
df.to_excel(r'dataset_one.xlsx', index=False, header=True)
# %%
# fit the model and get the separating hyperplane
clf = svm.SVC(kernel='linear', C=1.0)
clf.fit(X_s, y)

# fit the model and get the separating hyperplane using weighted classes
wclf = svm.SVC(kernel='linear', class_weight={1: 10})
wclf.fit(X_s, y)

fig, ax = plt.subplots()
# plot the samples
ax.scatter(X_s[:, 0], X_s[:, 1], c=y, cmap=plt.cm.Paired, edgecolors='k')

# plot the decision functions for both classifiers
ax = plt.gca()
xlim = ax.get_xlim()
ylim = ax.get_ylim()

# create grid to evaluate model
xx = np.linspace(xlim[0], xlim[1], 30)
yy = np.linspace(ylim[0], ylim[1], 30)
YY, XX = np.meshgrid(yy, xx)
xy = np.vstack([XX.ravel(), YY.ravel()]).T

# get the separating hyperplane
Z = clf.decision_function(xy).reshape(XX.shape)

# plot decision boundary and margins
a = ax.contour(XX, YY, Z, colors='k', levels=[0], alpha=0.5, linestyles=['-'])

# get the separating hyperplane for weighted classes
Z = wclf.decision_function(xy).reshape(XX.shape)

# plot decision boundary and margins for weighted classes
b = ax.contour(XX, YY, Z, colors='r', levels=[0], alpha=0.5, linestyles=['-'])

plt.legend([a.collections[0], b.collections[0]], ["non weighted", "weighted"],
           loc="upper right")
plt.show()
