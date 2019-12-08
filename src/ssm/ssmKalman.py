"""
This module implement the ssm from the textbook.
Time Series Analysis by State Space Methods : Time Series Analysis by State Space Methods.
"""
import numpy as np
import copy
import math
import sys

class ssmKalman(object):
    def __init__(self, seed = 0, P = 5, T = 100, M=3, lambda_1=0.1, lambda_2=0.1, lambda_3=0.1, positive_state=False, positive_em=False, verbose=False):
        """
        :param P: #assay 
        :param T: #bp
        :param M: #factor

        """
        self.P = P
        self.T = T
        self.M = M
        
        # initialize the matrix
        np.random.seed(seed)
        self.y_m = np.asmatrix(np.random.rand(self.P, self.T))
        self.alpha_m = np.abs(np.asmatrix(np.random.rand(self.M, self.T)))
        self.z_m = np.abs(np.asmatrix(np.random.rand(self.P, self.M)))
        self.t_m = np.asmatrix(np.eye(self.M))
        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2
        self.lambda_3 = lambda_3
        self.positive_state = positive_state
        self.positive_em = positive_em
        self.precision = sys.float_info.epsilon
        self.verbose = verbose
        self.error_m = []
        
    def initialize(self):
        self.alpha_m[:, 0] = np.abs(np.random.rand(self.M, 1))
        t = 0
        while t < self.T - 1:
            self.y_m[:, t] = (
                np.matrix(np.random.multivariate_normal((self.z_m * self.alpha_m[:, t]).flatten().tolist()[0], np.eye(self.P)))).T
            self.alpha_m[:, t + 1] = np.abs(
                np.matrix(np.random.multivariate_normal((self.t_m * self.alpha_m[:, t]).flatten().tolist()[0], np.eye(self.M)))).T
            t += 1
        self.y_m[:, t] = (
            np.matrix(np.random.multivariate_normal((self.z_m * self.alpha_m[:, t]).flatten().tolist()[0], np.eye(self.P)))).T
    
    def set_y(self, y_m):
        self.y_m = copy.deepcopy(np.asmatrix(y_m))

    def set_z(self, z_m):
        self.z_m = copy.deepcopy(np.asmatrix(z_m))

    def set_t(self, t_m):
        self.t_m = copy.deepcopy(np.asmatrix(t_m))

    def set_alpha(self, alpha_m):
        self.alpha_m = copy.deepcopy(np.asmatrix(alpha_m))

    def filter(self):
        for t in range(self.T):
            v_m = self.y_m[:, t] - self.z_m * self.param_dic["a_m"][t]
            f_m = self.z_m * self.param_dic["p_m"][t] * self.z_m.T + np.asmatrix(np.eye(self.P))
            k_m = self.t_m * self.param_dic["p_m"][t] * self.z_m.T * f_m.I
            a_m_next = self.t_m * self.param_dic["a_m"][t] + k_m * v_m
            p_m_next = self.t_m * self.param_dic["p_m"][t] * (self.t_m - k_m * self.z_m).T + self.lambda_1 * np.asmatrix(np.eye(self.M))    
            self.param_dic["v_m"].append(copy.deepcopy(v_m))
            self.param_dic["f_m"].append(copy.deepcopy(f_m))
            self.param_dic["k_m"].append(copy.deepcopy(k_m))
            self.param_dic["a_m"].append(copy.deepcopy(a_m_next))
            self.param_dic["p_m"].append(copy.deepcopy(p_m_next))

    def smooth(self):
        for t in reversed(range(self.T)):
            l_m = self.t_m - self.param_dic["k_m"][t] * self.z_m
            # from book
            r_m_last = self.z_m.T * self.param_dic["f_m"][t].I * self.param_dic["v_m"][t] + l_m.T * self.param_dic["r_m"][0]

            # # from file:///Users/bowenchen/Downloads/fulton_statsmodels_2017_v1.pdf, actually the old version of text book
            # r_m_last = self.z_m.T * (self.param_dic["f_m"][t].I * self.param_dic["v_m"][t] - self.param_dic["k_m"][t].T * self.param_dic["r_m"][0]) + self.t_m.T * self.param_dic["r_m"][0]
            new_alpha_m = self.param_dic["a_m"][t] + self.param_dic["p_m"][t] * r_m_last
            if self.positive_state:
                for _ in range(self.M):
                    # add the lagrange multiplyers for elements equal to 0
                    A = np.asmatrix(np.zeros((self.M,self.M)))
                    b = np.asmatrix(np.zeros((self.M,1)))
                    active = False
                    for i in range(self.M):
                        if new_alpha_m[i,0] < 0:
                            A[i,i] = 1
                            active = True
                    if active:
                        # print("before lag\n", new_alpha_m)
                        new_alpha_m = new_alpha_m - A.T*np.linalg.pinv(A*A.T)*(A*new_alpha_m - b)
                        # print("after lag", new_alpha_m)
                    else:
                        new_alpha_m = new_alpha_m
                    
                    for i in range(self.M):
                        if -self.precision < new_alpha_m[i, 0] < self.precision:
                            new_alpha_m[i, 0] = 0

                    v = - self.alpha_m[:,t] + new_alpha_m
                    min_k = 1 # alpha_start = alpha_m + kv
                    min_i = 0
                    for i in range(self.M):
                        if v[i,0] < 0 and self.alpha_m[i,t] != 0:
                            k = self.alpha_m[i,t] / (-v[i,0])
                            if k < min_k:
                                min_k = k
                                min_i = i
                    # print('k', min_k, 'i', min_i)
                    # print('origin alpha', self.alpha_m[:,t])
                    # print('v', v)
                    self.alpha_m[:,t] += v*min_k

                    for i in range(self.M):
                        if -self.precision < self.alpha_m[i,t] < self.precision:
                            self.alpha_m[i,t] = 0

                    # print('updated alpha', self.alpha_m[:,t])
                    if not self.check_positive(self.alpha_m[:,t]):
                        print(self.alpha_m[:,t])
                        raise Exception('Non positive state')
            else:
                self.alpha_m[:,t] = new_alpha_m
            self.param_dic["r_m"].insert(0, copy.deepcopy(r_m_last))

    def update_alpha(self):
        # initailize param dic
        self.param_dic = {"v_m":[], "f_m":[], "k_m":[], "a_m":[], "p_m":[],  "r_m":[]}
        # initailize a0 and P0 before filtering 
        # self.param_dic["a_m"].append(np.asmatrix(np.random.rand(self.M, 1)))
        self.param_dic["a_m"].append(copy.deepcopy(np.asmatrix(self.alpha_m[:,0])))
        self.param_dic["p_m"].append(np.asmatrix(np.eye(self.M)))
        # forward filtering
        self.filter()
        # initialize r_{n-1} before somoothing
        self.param_dic["r_m"].insert(0, np.asmatrix(np.zeros((self.M, 1))))
        # backwards smoothing
        self.smooth()
        
    def check_positive(self, a):
        m,n = a.shape
        for i in range(m):
            for j in range(n):
                if a[i,j] < 0:
                    return False
        return True
    
    def update_z(self, lhs_m, rhs_m):
        """
        using the samve active method, but
        - go update first.
        - then calculate the active set
        """
        t = 0
        new_z_m = (lhs_m* rhs_m).T
        z_m_noLag = new_z_m.copy()
        if self.positive_em:
            # print('rhs_m', rhs_m)
            new_z_m = np.asmatrix(np.zeros((self.M, self.P))).T
            for p in range(self.P):
                old_z = self.z_m[p, :].T
                lag_lambda_m = np.asmatrix(np.zeros((self.P, self.M)))
                for _ in range(self.M):
                    # print('[z]:before update\n', old_z)
                    # print('[z]:update without lag\n', z_m_noLag[p, :])
                    # print('[z]:lag term\n', lag_lambda_m[p,:].T)
                    new_z = lhs_m * (rhs_m[:,p] + 1/2*lag_lambda_m[p,:].T)
                    # print('[z]:after lag\n', new_z)

                    # print('after round up\n', new_z)
                    v = new_z - old_z
                    min_k = 1
                    min_m = 0
                    for m in range(self.M):
                        if v[m, 0] < 0:
                            k = self.z_m[p, m].T / (-v[m, 0])
                            if k < min_k:
                                min_k = k
                                min_m = m
                    old_z += min_k*v
                    for m in range(self.M):
                        if -self.precision <= old_z[m,0] <= self.precision:
                            old_z[m,0] = 0
                    # print('[z]:updated z', old_z[:, 0])
                    
                    # get active set
                    active = []
                    for m in range(self.M):
                        if old_z[m, 0] < 0:
                            active.append(m)
                    # print('[z]:active set', active)
                    # set lag term to 0
                    lag_lambda_m = np.asmatrix(np.zeros((self.P, self.M)))
                    
                    if len(active) > 0:
                        active_len = len(active)
                        lag_lambda = np.asmatrix(np.zeros((self.P, active_len)))
                        lag_z = np.asmatrix(np.zeros((self.M, self.P)))
                        lag_lhs = np.asmatrix(np.zeros((active_len, active_len)))
                        lag_lhs_bar = np.asmatrix(np.zeros((active_len, self.M - active_len)))
                        lag_rhs = np.asmatrix(np.zeros((active_len, 1)))
                        lag_rhs_bar = np.asmatrix(np.zeros((self.M - active_len, 1)))
                        for r, val in enumerate(active):
                            counter = 0
                            bar_counter = 0
                            for c in range(self.M):
                                if c in active:
                                    lag_lhs[r, counter] = lhs_m[val, c]
                                    counter += 1
                                else:
                                    lag_lhs_bar[r, bar_counter] = lhs_m[val, c]
                                    bar_counter += 1
                        counter = 0
                        bar_counter = 0
                        for c in range(self.M):
                            if c in active:
                                lag_rhs[counter, 0] = rhs_m[c, p]
                                counter += 1
                            else:
                                lag_rhs_bar[bar_counter, 0] = rhs_m[c, p]
                                bar_counter += 1
                        lag_lambda = (-2 * np.linalg.pinv(lag_lhs) * lag_lhs_bar * lag_rhs_bar - 2 * lag_rhs).T
                        # print('lag_lambda\n', lag_lambda)
                        # print('lag_lhs\n', lag_lhs)
                        # print('lag_lhs_bar\n', lag_lhs_bar)
                        # print('lag_rhs\n', lag_rhs)
                        # print('lag_rhs_bar\n', lag_rhs_bar)
                        expected_theta = lag_lhs * (lag_rhs + 1/2*lag_lambda.T) + lag_lhs_bar * lag_rhs_bar
                        # print('expected theta\n', expected_theta) 
                        
                        for c, val in enumerate(active):
                            lag_lambda_m[:, val] = lag_lambda[:, c]
                if not self.check_positive(old_z):
                    print(old_z)
                    raise Exception('[z]:Non positive state')
                new_z_m[p, :] = old_z[:, 0].T
        self.z_m = copy.deepcopy(new_z_m)

    def get_z_msg(self):
        t = 0
        lhs_m = np.asmatrix(np.zeros((self.M, self.M)))
        rhs_m = np.asmatrix(np.zeros((self.M, self.P)))
        while t < self.T:
            lhs_m += self.alpha_m[:, t] * self.alpha_m[:, t].T
            rhs_m += self.alpha_m[:, t] * self.y_m[:, t].T
            t += 1
        lhs_m += self.lambda_2 * np.eye(self.M)
        lhs_m = np.linalg.pinv(lhs_m)
        return lhs_m, rhs_m
    
    def get_t_msg(self):
        t = 0
        lhs_m = np.asmatrix(np.zeros((self.M, self.M)))
        rhs_m = np.asmatrix(np.zeros((self.M, self.M)))
        while t < self.T - 1:
            lhs_m += self.alpha_m[:,t] * self.alpha_m[:,t].T
            rhs_m += self.alpha_m[:,t] * self.alpha_m[:,t+1].T
            t += 1
        return lhs_m, rhs_m

    def update_t(self):
        lhs_m, rhs_m = self.get_t_msg()
        self.t_m = np.asmatrix(np.linalg.lstsq(lhs_m + self.lambda_3 * np.eye(self.M), rhs_m, rcond=1)[0]).T

    def total_error(self):
        error = 0
        for t in range(self.T):
            error += np.sum((self.y_m[:,t] - self.z_m * self.alpha_m[:, t]).T * (self.y_m[:,t] - self.z_m * self.alpha_m[:, t]))
            error += np.sum(self.lambda_1 * self.alpha_m[:,t].T * self.alpha_m[:,t])
            if t != 0:
                error += np.sum((self.alpha_m[:, t] - self.t_m * self.alpha_m[:, t-1]).T * (self.alpha_m[:, t] - self.t_m * self.alpha_m[:, t-1]))
        error += np.sum(self.lambda_2*np.multiply(self.z_m, self.z_m))
        error += np.sum(self.lambda_3*np.multiply(self.t_m, self.t_m))
        return error

    def optimization(self, iteration=15):
        self.error_m.append(self.total_error())
        for i in range(iteration): 
            self.update_alpha()
        lhs_m, rhs_m = self.get_z_msg()
        self.update_z(lhs_m, rhs_m)
        self.error_m.append(self.total_error())

if __name__ == "__main__":
    seed = 0
    M = 3
    T = 100
    P = 5
    lambda_3 = 0.1
    lambda_2 = 0.1
    lambda_1 = 0.1
    state_sum = 1
    # experiment
    model_pos = ssmKalman(seed=seed, P=P, T=T, M=M, lambda_1=lambda_1, 
        lambda_2=lambda_2, lambda_3=lambda_3, positive_state=True, positive_em=True)

    test_iteration = 40
    model_pos.optimization(test_iteration)

    # plot
    import time
    import numpy as np
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    # plot, predict vs real trans value
    fig, axarr = plt.subplots()
    axarr.plot(model_pos.error_m, label='pos')

    axarr.legend(loc='upper right')
    plt.xlabel('iteration')
    plt.ylabel('error')
    plt.show()
