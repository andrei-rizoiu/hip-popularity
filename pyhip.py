# -*- coding: utf-8 -*-
"""
==============================================================================
The code for the class pyhip
==============================================================================
This is the main class of python Hawkes Intensity Process model.
It provides fitting and forecasting functions for two time series data.
"""

from __future__ import print_function, division
import numpy as np
from scipy import optimize
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
# pretty font size
mpl.rcParams.update({'axes.titlesize': 16,
                     'font.size': 12})
# disable warning display to avoid flush of scalar overflow messages
warnings.filterwarnings("ignore")


class HIP(object):
    """ The main class of python Hawkes Intensity Process model.
    It provides fitting and forecasting functions for two time series data.
    """
    def __init__(self):
        # data parameters, set by self.initial function
        self.x = None
        self.y = None
        self.num_train = None
        self.num_test = None
        self.num_cv_train = None
        self.num_cv_test = None
        self.num_initialization = None

        # model parameters, learnt by self.fit_with_bfgs function
        self.mu = None
        self.theta = None
        self.C = None
        self.c = None
        self.gamma = None
        self.eta = None
        self.endo = None
        self.viral = None

    # initialize or reset model
    def initial(self, x, y, num_train=90, num_test=30, num_initialization=10):
        """ Initialize or reset model with data.
        :param x: exogenous stimuli signal
        :param y: HIP response signal
        :param num_train: number of training data
        :param num_test: number of test data, immediately after training data
        :param num_initialization: number of initialization sets to avoid local minimal
        """
        self.__init__()
        self.x = x
        self.y = y
        self.num_train = num_train
        self.num_test = num_test
        self.num_cv_train = int(num_train*0.8)
        self.num_cv_test = num_train - self.num_cv_train
        self.num_initialization = num_initialization

    # set parameters to model
    def set_parameters(self, params):
        """ Set array parameters to model.
        :param params: model parameters in an array structure
        """
        self.mu, self.theta, self.C, self.c, self.gamma, self.eta = params
        # set endo and viral
        self.endo = self.get_endo()
        self.viral = self.mu * self.endo

    # get parameters from model
    def get_parameters(self):
        """ Get parameters from model.
        :return: model parameters in the order of mu, theta, C, c, gamma, eta, endo, viral
        """
        return np.array([self.mu, self.theta, self.C, self.c, self.gamma, self.eta, self.endo, self.viral])

    # get parameters from model
    def get_parameters_abbr(self):
        """ Get abbreviate parameters from model.
        :return: model parameters in the order of mu, theta, C, c, gamma, eta
        """
        return np.array([self.mu, self.theta, self.C, self.c, self.gamma, self.eta])

    def get_endo(self):
        """ Get endogenous response from model.
        :return: endogenous response value
        """
        x_predict = np.zeros(10000)
        x_predict[0] = 1
        for i in range(1, 10000):
            x_predict[i] = self.C * np.sum(x_predict[:i] * (self.time_decay_base(i, self.c) ** (-1 - self.theta)))
        return np.sum(x_predict)

    def print_parameters(self):
        """ Print model parameters.
        """
        print('--- mu={0:.2f}, theta={1:.2f}, C={2:.2f}\n'
              '--- c={3:.2f}, gamma={4:.2f}, eta={5:.2f}\n'
              '--- endo={6:.2f}, viral={7:.2f}'
              .format(self.mu, self.theta, self.C, self.c, self.gamma, self.eta, self.endo, self.viral))

    # == == == == == == == == modelling components == == == == == == == == #
    def _rand_initialize_parameters(self):
        """ Initialize random parameters for HIP model.
        The initial bounds are estimated from a large collection of time series data.
        :return: a set of randomized parameters in the order of mu, theta, C, c, gamma, eta
        """
        rand_mu = np.random.uniform(0, 505.90)
        rand_theta = np.random.uniform(2.3, 67.7)
        rand_C = self._get_C(np.random.uniform(0, 52.9))
        rand_c = np.random.uniform(0, 4)
        rand_gamma = np.random.uniform(0, 9947)
        rand_eta = np.random.uniform(0, 289.2)
        return np.array([rand_mu, rand_theta, rand_C, rand_c, rand_gamma, rand_eta])

    @staticmethod
    def _get_C(k, alpha=2.016, beta=0.1):
        """ Get parameter capital C.
        :param k: scaling factor for video quality
        :param alpha: power-law exponent of user influence distribution
        :param beta: user influence component
        :return: parameter capital C
        """
        return k * (alpha - 1) / (alpha - beta - 1)

    @staticmethod
    def time_decay_base(i, c):
        """ Time decay kernel base for series (tau + c).
        :param i: tau value
        :param c: c value
        :return: bounded time lag series from start time to previous time
        """
        return np.arange(i, 0, -1) + c

    def predict(self, params, x):
        """ Predict viewcount given sharecount sequence x.
        :param params: model parameters, mu, theta, C, c, gamma, eta
        :param x: observed sharecount sequence from start time
        :return: predicted viewcount value
        """
        mu, theta, C, c, gamma, eta = params
        n = len(x)
        x_predict = np.zeros(len(x))
        x_predict[0] = gamma + mu * x[0]
        for i in range(1, n):
            x_predict[i] = eta + mu * x[i] + C * np.sum(x_predict[:i] * (self.time_decay_base(i, c) ** (-1 - theta)))
        return x_predict

    def cost_function(self, params, x, y, params0=None, num_split=None):
        """ MSE as cost function for HIP model.
        :param params: model parameters, mu, theta, C, c, gamma, eta
        :param x: observed sharecount
        :param y: observed viewcount
        :param params0: reference values from non-regularized model
        :param num_split: number of test set
        :return: cost function value
        """
        view_predict = self.predict(params, x)
        cost_vector = view_predict - y
        if num_split is not None:
            cost_vector = cost_vector[-num_split:]
        cost = np.sum(cost_vector ** 2) / 2
        if params0 is not None:
            # add one smooth to handle refer parameters equal to zero
            for i in range(4):
                if params0[i] == 0:
                    params0[i] = 1
            mu, theta, C, c, gamma, eta = params
            mu0, C0, gamma0, eta0, w0 = params0
            cost += w0 / 2 * ((mu / mu0) ** 2 + (C / C0) ** 2 + (gamma / gamma0) ** 2 + (eta / eta0) ** 2)
        return cost / len(cost_vector)

    def _compute_fitting_error(self):
        """ Get fitting RMSE of training data.
        :return: fitting RMSE of training data
        """
        return np.sqrt(self.cost_function(self.get_parameters_abbr(), self.x[:self.num_train], self.y[:self.num_train]))

    def _compute_forecast_error(self):
        """ Get forecast RMSE of test data.
        :return: forecast RMSE of test data
        """
        return np.sqrt(self.cost_function(self.get_parameters_abbr(), self.x[:self.num_train + self.num_test],
                                          self.y[:self.num_train + self.num_test], num_split=self.num_test))

    def grad_descent(self, params, x, y, params0=None):
        """ Gradient function for HIP model.
        :param params: model parameters, mu, theta, C, c, gamma, eta
        :param x: observed sharecount
        :param y: observed viewcount
        :param params0: reference values from non-regularized model
        :return: gradient descent function value
        """
        mu, theta, C, c, gamma, eta = params
        if params0 is not None:
            # add one smooth to handle refer parameters equal to zero
            for i in range(4):
                if params0[i] == 0:
                    params0[i] = 1
            mu0, C0, gamma0, eta0, w0 = params0
        else:
            mu0, C0, gamma0, eta0, w0 = 1, 1, 1, 1, 0
        view_predict = self.predict(params, x)
        n = len(x)
        # partial derivative for mu
        grad_mu_vector = np.zeros(n)
        grad_mu_vector[0] = x[0]
        for i in range(1, n):
            grad_mu_vector[i] = x[i] + C * np.sum(grad_mu_vector[:i] * (self.time_decay_base(i, c) ** (-1 - theta)))
        grad_mu = np.sum((view_predict - y) * grad_mu_vector) + w0 * mu / mu0 / mu0
        # partial derivative for theta
        grad_theta_vector = np.zeros(n)
        grad_theta_vector[0] = 0
        for i in range(1, n):
            grad_theta_vector[i] = C * np.sum((grad_theta_vector[:i] - view_predict[:i] * np.log(self.time_decay_base(i, c)))
                                              * (self.time_decay_base(i, c) ** (-1 - theta)))
        grad_theta = np.sum((view_predict - y) * grad_theta_vector)
        # partial derivative for C
        grad_C_vector = np.zeros(n)
        grad_C_vector[0] = 0
        for i in range(1, n):
            grad_C_vector[i] = np.sum((C * grad_C_vector[:i] + view_predict[:i])
                                      * (self.time_decay_base(i, c) ** (-1 - theta)))
        grad_C = np.sum((view_predict - y) * grad_C_vector) + w0 * C / C0 / C0
        # partial derivative for c
        grad_c_vector = np.zeros(n)
        grad_c_vector[0] = 0
        for i in range(1, n):
            grad_c_vector[i] = C * np.sum((grad_c_vector[:i] - (1 + theta) * view_predict[:i] / self.time_decay_base(i, c))
                                          * (self.time_decay_base(i, c) ** (-1 - theta)))
        grad_c = np.sum((view_predict - y) * grad_c_vector)
        # partial derivative for gamma
        grad_gamma_vector = np.zeros(n)
        grad_gamma_vector[0] = 1
        for i in range(1, n):
            grad_gamma_vector[i] = C * np.sum(grad_gamma_vector[:i] * (self.time_decay_base(i, c) ** (-1 - theta)))
        grad_gamma = np.sum((view_predict - y) * grad_gamma_vector) + w0 * gamma / gamma0 / gamma0
        # partial derivative for eta
        grad_eta_vector = np.zeros(n)
        grad_eta_vector[0] = 0
        for i in range(1, n):
            grad_eta_vector[i] = 1 + C * np.sum(grad_eta_vector[:i] * (self.time_decay_base(i, c) ** (-1 - theta)))
        grad_eta = np.sum((view_predict - y) * grad_eta_vector) + w0 * eta / eta0 / eta0
        return np.array([grad_mu, grad_theta, grad_C, grad_c, grad_gamma, grad_eta]) / n

    def fit_with_bfgs(self):
        """ Fit HIP with BFGS optimization tool.
        """
        # bound all parameter to be positive
        bounds = [(0, None), (0, None), (0, None), (0, None), (0, None), (0, None)]
        best_cost = np.inf
        best_initial_weight = None
        best_params0 = None

        # use cross validation to find best initialization and best regularization reference params0
        for t in range(self.num_initialization):
            # perform non-regularized optimization with l-bfgs
            initial_weight = self._rand_initialize_parameters()
            init_optimizer = optimize.minimize(self.cost_function, initial_weight, jac=self.grad_descent,
                                               method='L-BFGS-B',
                                               args=(self.x[:self.num_cv_train], self.y[:self.num_cv_train]),
                                               bounds=bounds)

            mu0, theta0, C0, c0, gamma0, eta0 = init_optimizer.x
            J0 = init_optimizer.fun
            # line search in logspace (10e-4*J0, 10*J0)
            for w0 in np.arange(np.log(10 ** -4 * J0), np.log(10 * J0), 1):
                w0 = np.exp(w0)
                reg_params0 = np.array([mu0, C0, gamma0, eta0, w0])
                reg_optimizer = optimize.minimize(self.cost_function, init_optimizer.x, jac=self.grad_descent,
                                                  method='L-BFGS-B',
                                                  args=(self.x[:self.num_cv_train], self.y[:self.num_cv_train], reg_params0),
                                                  bounds=bounds)
                # model selection by using cv dataset, here we consider fitting error and forecast error
                cv_cost = self.cost_function(reg_optimizer.x, self.x[:self.num_train], self.y[:self.num_train],
                                             params0=reg_params0)
                if cv_cost < best_cost:
                    best_cost = cv_cost
                    best_initial_weight = reg_optimizer.x
                    best_params0 = reg_params0

            # display process bar
            print('--- Finish initialization set {0}...'.format(t+1))

        # re-train on the first self.num_train days
        best_optimizer = optimize.minimize(self.cost_function, best_initial_weight, jac=self.grad_descent,
                                           method='L-BFGS-B',
                                           args=(self.x[:self.num_train], self.y[:self.num_train], best_params0),
                                           bounds=bounds)
        self.set_parameters(best_optimizer.x)
        print('--- Model fitting RMSE: {0:.2f}'.format(self._compute_fitting_error()))
        print('--- Model forecast RMSE: {0:.2f}'.format(self._compute_forecast_error()))

    # plot function for fitting and forecasting process
    def plot_func(self, title):
        """ Plot fitting and forecasting of HIP model
        :param title: figure title, YoutubeID
        """
        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()
        age = self.num_train + self.num_test
        ax1.plot(np.arange(age), self.y[:age], 'k--', label='observed #views')
        ax2.plot(np.arange(age), self.x[:age], 'r-', label='#share')
        ax1.plot((self.num_train, self.num_train), (ax1.get_ylim()[0], ax1.get_ylim()[1]), 'k--')

        ax1.set_xlim(xmin=0)
        ax1.set_xlim(xmax=age)
        ax1.set_ylim(ymin=max(0, ax1.get_ylim()[0]))
        ax2.set_ylim(ymin=max(0, ax2.get_ylim()[0]))
        ax2.set_ylim(ymax=3 * max(self.x))
        ax1.set_xlabel('video age (day)')
        ax1.set_ylabel('Number of views', color='k')
        ax1.tick_params('y', colors='k')
        ax2.set_ylabel('Number of shares', color='r')
        ax2.tick_params('y', colors='r')

        ax2.text(0.03, 0.85, '$\mu$={0:.2f}, $\\theta$={1:.2f}, C={2:.2f}\n'
                             'c={3:.2f}, $\gamma$={4:.2f}, $\eta$={5:.2f}\n'
                             'endo={6:.2f}, viral={7:.2f}'.format(*self.get_parameters()), transform=ax1.transAxes)
        ax1.set_title(title)

        predicted_x = self.predict(self.get_parameters_abbr(), self.x)
        ax1.plot(np.arange(self.num_train), predicted_x[:self.num_train], 'b-', label='HIP fit')
        ax1.plot(np.arange(self.num_train, age), predicted_x[self.num_train:age], 'g-', label='HIP forecast')

        plt.legend([plt.Line2D((0, 1), (0, 0), color='k', linestyle='--'),
                    plt.Line2D((0, 1), (0, 0), color='b'),
                    plt.Line2D((0, 1), (0, 0), color='g'),
                    plt.Line2D((0, 1), (0, 0), color='r')],
                   ['Observed view', 'Fitted view', 'Forecast view', 'Observed share'],
                   frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.125),
                   fancybox=True, shadow=True, ncol=4)

        plt.tight_layout(rect=[0, 0.04, 1, 1])
        plt.show()
