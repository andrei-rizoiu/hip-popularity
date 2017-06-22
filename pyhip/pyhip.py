from __future__ import print_function, division
import sys
import os
import bz2
import json
import cPickle as pickle
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# Hawkes Intensity Process model in Python


def get_C(k, alpha=2.016, beta=0.1):
    """
    Get parameter capital C.
    :param k: scaling factor for video quality
    :param alpha: power-law exponent of user influence distribution
    :param beta: user influence component
    :return: parameter capital C
    """
    return k*(alpha-1)/(alpha-beta-1)


def rand_initialize_weights(n):
    """
    Initialize multiple sets of random weights for theta.
    :param n: number of sets of random weights
    :return: n sets of random vectors, in the order of mu, theta, C, c, gamma, eta
    """
    ret = []
    for _ in xrange(n):
        rand_mu = np.random.uniform(0, 505.90)
        rand_theta = np.random.uniform(2.3, 67.7)
        rand_C = get_C(np.random.uniform(0, 52.9))
        rand_c = np.random.uniform(0, 4)
        rand_gamma = np.random.uniform(0, 9947)
        rand_eta = np.random.uniform(0, 289.2)
        ret.append(np.array([rand_mu, rand_theta, rand_C, rand_c, rand_gamma, rand_eta]))
    return ret


def time_decay(i, c):
    """
    Time decay part for series (tau + c).
    :param i: tau value
    :param c: c value
    :return: abbreviated presentation
    """
    return np.arange(1, i+1)[::-1]+c


def predict(params, x):
    """
    Predict viewcount with sharecount sequence x.
    Comments are for vector operation style
    :param params: model parameters, mu, theta, C, c, gamma, eta
    :param x: observed sharecount sequence from beginning
    :return: predict value
    """
    mu, theta, C, c, gamma, eta = params
    n = len(x)
    x_predict = np.zeros(len(x))
    for i in xrange(n):
        if i == 0:
            x_predict[0] = gamma + mu*x[0]
        else:
            x_predict[i] = eta + mu*x[i] + C*np.sum(x_predict[:i]*(time_decay(i, c)**(-1-theta)))
    return x_predict


def cost_function(params, x, y, num_split=None):
    """
    Non-regularized cost function for HIP model
    :param params: model parameters, mu, theta, C, c, gamma, eta
    :param x: observed sharecount
    :param y: observed viewcount
    :param num_split: number of test set
    :return: cost function value
    """
    view_predict = predict(params, x)
    cost_vector = view_predict - y
    if num_split is not None:
        cost_vector = cost_vector[-num_split:]
    cost = np.sum(cost_vector ** 2) / 2
    return cost/len(cost_vector)


def grad_descent(params, x, y):
    """
    Non-regularized gradient function for HIP model
    :param params: model parameters, mu, theta, C, c, gamma, eta
    :param x: observed sharecount
    :param y: observed viewcount
    :return: cost function value
    """
    mu, theta, C, c, gamma, eta = params
    view_predict = predict(params, x)
    n = len(x)
    # partial derivative for mu
    grad_mu_vector = np.zeros(n)
    grad_mu_vector[0] = x[0]
    for i in xrange(1, n):
        grad_mu_vector[i] = x[i] + C*np.sum(grad_mu_vector[:i] * (time_decay(i, c)**(-1-theta)))
    grad_mu = np.sum((view_predict-y)*grad_mu_vector)
    # partial derivative for theta
    grad_theta_vector = np.zeros(n)
    grad_theta_vector[0] = 0
    for i in xrange(1, n):
        grad_theta_vector[i] = C*np.sum((grad_theta_vector[:i]-view_predict[:i]*np.log(time_decay(i, c))) * (time_decay(i, c)**(-1-theta)))
    grad_theta = np.sum((view_predict-y)*grad_theta_vector)
    # partial derivative for C
    grad_C_vector = np.zeros(n)
    grad_C_vector[0] = 0
    for i in xrange(1, n):
        grad_C_vector[i] = np.sum((C*grad_C_vector[:i]+view_predict[:i]) * (time_decay(i, c)**(-1-theta)))
    grad_C = np.sum((view_predict-y)*grad_C_vector)
    # partial derivative for c
    grad_c_vector = np.zeros(n)
    grad_c_vector[0] = 0
    for i in xrange(1, n):
        grad_c_vector[i] = C*np.sum((grad_c_vector[:i]-(1+theta)*view_predict[:i]/time_decay(i, c)) * (time_decay(i, c)**(-1-theta)))
    grad_c = np.sum((view_predict-y)*grad_c_vector)
    # partial derivative for gamma
    grad_gamma_vector = np.zeros(n)
    grad_gamma_vector[0] = 1
    for i in xrange(1, n):
        grad_gamma_vector[i] = C*np.sum(grad_gamma_vector[:i] * (time_decay(i, c)**(-1-theta)))
    grad_gamma = np.sum((view_predict-y)*grad_gamma_vector)
    # partial derivative for eta
    grad_eta_vector = np.zeros(n)
    grad_eta_vector[0] = 0
    for i in xrange(1, n):
        grad_eta_vector[i] = 1 + C*np.sum(grad_eta_vector[:i] * (time_decay(i, c)**(-1-theta)))
    grad_eta = np.sum((view_predict-y)*grad_eta_vector)
    return np.array([grad_mu, grad_theta, grad_C, grad_c, grad_gamma, grad_eta])/n


def train_process(x_train, y_train, initial_weights_sets):
    """
    Train HIP with BFGS optimization tool
    :param x_train: train sharecount
    :param y_train: train viewcount
    :param initial_weights_sets: sets of random initial weights
    :return: best optimization parameters
    """
    best_params = None
    best_cost = np.inf

    for init_idx, initial_weight in enumerate(initial_weights_sets):
        # perform non-regularized optimization with l-bfgs
        optimizer = optimize.minimize(cost_function, initial_weight, jac=grad_descent, method='L-BFGS-B',
                                      args=(x_train, y_train), bounds=bounds)
        if optimizer.fun < best_cost:
            best_cost = optimizer.fun
            best_params = optimizer.x

    return best_params


def plot_func(params, x, y, title, idx):
    """
    Plot trend from R-HIP, PY-HIP and AUTO-HIP parameters
    :param params: model parameters, mu, theta, C, c, gamma, eta
    :param x: observed sharecount
    :param y: observed viewcount
    :param title: figure title, YoutubeID
    :param idx: subplot index
    :return:
    """
    # visualise sample data
    ax1 = fig.add_subplot(121+idx)
    # ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(np.arange(1, age+1), y, 'k--', label='observed #views')
    ax2.plot(np.arange(1, age+1), x, 'r-', label='#share')
    ax1.plot((num_train, num_train), (ax1.get_ylim()[0], ax1.get_ylim()[1]), 'k--')

    ax1.set_ylim(ymin=max(0, ax1.get_ylim()[0]))
    ax2.set_ylim(ymax=3*max(x))
    ax1.set_xlabel('video age (day)')
    ax1.set_ylabel('Number of views', color='k')
    ax1.tick_params('y', colors='k')
    ax2.set_ylabel('Number of shares', color='r')
    ax2.tick_params('y', colors='r')

    mu, theta, C, c, gamma, eta = params
    ax2.text(0.03, 0.85, '$\mu$={0:.2f}, $\\theta$={1:.2f}\nC={2:.2f}, c={3:.2f}\n$\gamma$={4:.2f}, $\eta$={5:.2f}'
             .format(mu, theta, C, c, gamma, eta), transform=ax1.transAxes)
    ax1.set_title(title)

    predidt_x = predict(params, x)
    ax1.plot(np.arange(1, num_train+1), predidt_x[:num_train], 'b-', label='HIP fit')
    ax1.plot(np.arange(num_train+1, age+1), predidt_x[num_train:age], 'm-', label='HIP forecast')


if __name__ == '__main__':
    # == == == == == == == == Part 1: Load ACTIVE dataset == == == == == == == == #
    # First time it gets loaded from the JSON format and writes essential fields into a pickle binary file.
    # check if the binary exists
    if not os.path.exists('../data/active-dataset.p'):
        print('--> Converting ACTIVE dataset from JSON format to pickle... might take a while!')
        test_cases = {}
        with bz2.BZ2File('../data/active-dataset.json.bz2') as f:
            dataset = json.loads(f.readline())
            for video in dataset:
                test_cases[video['YoutubeID']] = (video['numShare'], video['dailyViewcount'])
        pickle.dump(test_cases, open('../data/active-dataset.p', 'wb'))

    print('--> Loading the ACTIVE dataset from pickle...')
    test_cases = pickle.load(open('../data/active-dataset.p', 'rb'))
    # select 2 videos from paper
    test_vids = ['bUORBT9iFKc', 'cG0nQTYd8ck']
    # or random select 2 videos
    # test_videos = np.array(test_cases.keys())
    # random_index = np.random.randint(0, len(test_videos), 2)
    # test_vids = test_videos[random_index]

    # == == == == == == == == Part 2: Set up experiment parameters == == == == == == == == #
    # setting parameters
    fig = plt.figure(figsize=(14, 5))
    age = 120
    num_train = 90
    num_test = 30
    k = 5
    bounds = [(0, None), (0, None), (0, None), (None, None), (0, None), (0, None)]

    for tc_idx, vid in enumerate(test_vids):
        print('fitting and forecasting for video: {0}'.format(vid))
        dailyshare, dailyview = test_cases[vid]
        dailyshare = dailyshare[:age]
        dailyview = dailyview[:age]

        x_train = dailyshare[: num_train]
        y_train = dailyview[: num_train]

        # initialize weights
        # k sets of random params
        initial_weights_sets = rand_initialize_weights(k)

        # == == == == == == == == Part 3: Train with closed form gradient == == == == == == == == #
        best_fitted_params = train_process(x_train, y_train, initial_weights_sets)

        # == == == == == == == == Part 4: Plot fitting and forecast result == == == == == == == == #
        plot_func(best_fitted_params, dailyshare, dailyview, vid, tc_idx)

    plt.tight_layout()
    plt.show()
