# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 10:36:38 2021

@author: mahdi
"""

import numpy as np
from scipy.linalg import toeplitz
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
from matplotlib.pyplot import figure
import matplotlib.colors as mcolors
import matplotlib as mpl
from numpy import linalg as LA

# %% Figure settings

# figure(num=None, figsize=(8, 7), dpi=100, facecolor='w', edgecolor='k')
# plt.rcParams['figure.figsize'] = (13, 9)
plt.style.use(['default'])
# plt.style.use('dracula.mplstyle')
rc('font', **{'family': 'serif', 'serif': ['Times']})
font = {'size': 9}
mpl.rc('font', **font)
plt.rcParams['font.size'] = '9'
plt.rcParams["font.family"] = "Times New Roman"

# %% Functions


def generate_A(n):
    """
    Generate the A toeplitz matrix of input.

    Parameters
    ----------
    n : int
        Length of the input data.

    Returns
    -------
    A: numpy arra
       The toeplitz input matrix.

    """
    # Bernouli sequence as input
    U = np.random.binomial(size=n, n=1, p=0.5)
    # U = np.arange(1, 6)
    for i in range(len(U)):
        if U[i] == 0:
            U[i] = -1

    A = toeplitz(U)
    n_row = A.shape[1]

    for i in range(n_row):
        A[i+1:, i] = 0
    A = np.transpose(A)
    return A * 10


def parameters_t(m):
    """
    Generate the paramters vector.

    Parameters
    ----------
    m : int
        length of the parameter.

    Returns
    -------
    None.

    """
    param_vec = np.zeros(m)
    for i in range(m-1):
        param_vec[i+1] = 0.3 * np.power(0.5, i) + 3 * i * np.power(0.8, i)
    return param_vec


# %%
data_length = 100
A = generate_A(data_length)
theta_vec = parameters_t(data_length)

fig, ax = plt.subplots()
ax.stem(theta_vec)
db_r = 15  # SRN in dB
y_bar = A @ theta_vec
sigma_2 = ((np.sum(np.power(y_bar, 2)))/len(y_bar))/np.power(10, db_r/10)
sigma = np.sqrt(sigma_2)
w = np.random.normal(0, sigma, len(y_bar))
y = y_bar + w

# %% setting parameters
m_steps = 10  # m in the paper, range of the maximum order searching
n_trials = 3  # number of trials to average over
alpha = 4
beta = 4
db_vec = np.arange(0, 20, 0.5)
Zsm_upmat = np.zeros((m_steps, len(db_vec)), dtype=np.csingle)
Zsm_lomat = np.zeros((m_steps, len(db_vec)), dtype=np.csingle)
c = 0
# Zsm_mat[0, :] = np.transpose(db_vec)
for db in db_vec:
    # db_temp = 10
    sigma_2 = ((np.sum(np.power(y_bar, 2)))/len(y_bar))/(np.power(10,
                                                                  db/10))
    sigma = np.sqrt(sigma_2)
    Xsm_vec = np.zeros((m_steps, n_trials), dtype=np.csingle)
    Jsm_vec = np.zeros((m_steps, 1), dtype=np.csingle)
    Zsm_upvec = np.zeros((m_steps, n_trials), dtype=np.csingle)
    Zsm_lovec = np.zeros((m_steps, n_trials), dtype=np.csingle)
    for m in range(n_trials):
        for i in range(m_steps):
            Asm = A[:, 0:i+1]
            theta_m = theta_vec[0:i+1]
            theta_hat = np.linalg.inv(Asm.transpose() @ Asm) @ Asm.transpose()\
                @ y
            theta_hat.resize(len(y_bar))
            # Asm_temp = np.hstack((Asm, np.zeros((data_length,
            #                                     data_length-(i+1)))))
            y_hat = A @ theta_hat
            Xsm = (np.power(LA.norm((y - y_hat), 2), 2))/data_length
            Xsm_vec[i, m] = Xsm
            Jsm = np.power(LA.norm((theta_hat - theta_vec), 2), 2)
            Jsm_vec[i] = Jsm
            mw = (1 - ((i+1)/data_length)) * sigma_2
            Ksm = (2 * alpha * sigma / np.sqrt(data_length)) * np.sqrt(
                np.power((alpha*sigma), 2) + Xsm - (mw/2) + 0.j)
            Usm = Xsm - mw + (2 * np.power((alpha * sigma), 2)/data_length)\
                + Ksm
            Lsm = Xsm - mw + (2 * np.power((alpha * sigma), 2)/data_length)\
                - Ksm
            Zsm_up = Usm + (i+1)/data_length * sigma_2 +\
                beta * np.sqrt(2 * m) * sigma_2 / data_length
            Zsm_lo = Lsm + (i+1)/data_length * sigma_2 -\
                beta * np.sqrt(2 * m) * sigma_2 / data_length
            Zsm_upvec[i, m] = Zsm_up
            Zsm_lovec[i, m] = Zsm_lo
    Xsm_mean = np.mean(Xsm_vec, axis=1)
    Zsm_upmean = np.mean(Zsm_upvec, axis=1)
    Zsm_lomean = np.mean(Zsm_lovec, axis=1)
    Zsm_upmat[:, c] = Zsm_upmean
    Zsm_lomat[:, c] = Zsm_lomean
    c = c+1
