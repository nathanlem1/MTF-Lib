# This script runs KF example
#
# Version 1.0, January, 2019
#
# This function was written by Nathanael L. Baisa
# email: nathanaellmss@gmail.com
#
#
# State vector is position x,y and velocity vx, vy i.e. [x, vx, y,vy] and observation vector is position x, y i.e.
# [x, y]


from __future__ import division, print_function, unicode_literals  # To support both python 2 and python 3
import numpy as np
import matplotlib.pyplot as plt
from KalmanFilter import *


class KF_demo(object):

    def __init__(self):
        # Dynamical model parameters(Constant Velocity (CV) model)
        self.model_T = {"T": 1}                                          # sampling period
        self.model_A0 = {"A0": np.array([[1, self.model_T["T"]], [0, 1]])}
        self.model_F = {"F": np.concatenate((np.concatenate((self.model_A0["A0"], np.zeros((2,2))),axis=1), np.concatenate((np.zeros((2,2)), self.model_A0["A0"]),axis=1)),axis=0)}  # transition matrix
        self.model_B0 = {"B0": np.array([[float((self.model_T["T"]**2))/2], [self.model_T["T"]]])}
        self.model_B = {"B": np.concatenate((np.concatenate((self.model_B0["B0"], np.zeros((2,1))),axis=0),np.concatenate((np.zeros((2,1)), self.model_B0["B0"]),axis=0)),axis=1)}  # process noise standard deviation
        self.model_sigma_v = {"sigma_v": 2.5}
        self.model_Q = {"Q": self.model_sigma_v["sigma_v"]**2*self.model_B["B"].dot(np.transpose(self.model_B["B"]))}   # process noise covariance

        # Observation model parameters(noisy x / y only)
        self.model_H = {"H": np.array([[1, 0, 0, 0], [0, 0, 1, 0]])}    # observation matrix
        self.model_D = {"D": np.diag(np.array([10, 10]))}              # observation noise standard deviation
        self.model_R = {"R": self.model_D["D"].dot(np.transpose(self.model_D["D"]))}    # observation noise covariance

        # Combine the models together into one model: 'model'
        self.model = self.model_T.copy()
        self.model.update(self.model_A0), self.model.update(self.model_F), self.model.update(self.model_B0),
        self.model.update(self.model_B), self.model.update(self.model_sigma_v), self.model.update(self.model_Q),
        self.model.update(self.model_H), self.model.update(self.model_D), self.model.update(self.model_R)

        print(self.model["F"])
        print(self.model["Q"])
        print(self.model["H"])
        print(self.model["R"])

        # Birth parameters
        self.m_birth = np.array([[10],  [0],  [-10],  [0]])  # mean of Gaussian birth term
        self.B_birth = np.diag([10, 10, 10, 10])  # std of Gaussian birth term
        self.P_birth = self.B_birth.dot(np.transpose(self.B_birth))  # cov of Gaussian birth term

        self.xstart = self.m_birth
        self.targetstate = self.xstart

        self.N_duration = 100 # length of data / number of scans
        self.truth_X = [] # ground truth for state of target
        self.meas_Z = []  # generated measurement of a target

        # Initialization
        self.m_update = self.m_birth
        self.P_update = self.P_birth

        # Storing estimated states
        self.estimated_X = []

    def KF_demo_run(self):

        model = self.model
        truth_X = self.truth_X
        meas_Z = self.meas_Z
        estimated_X = self.estimated_X
        m_update = self.m_update
        P_update = self.P_update
        targetstate = self.targetstate

        for k in range(self.N_duration):

            # simulate dynamic target state
            # targetstate = model["F"].dot(targetstate) + np.sqrt(model["Q"]).dot(np.random.randn(targetstate.shape[0],targetstate.shape[1]))  # generate ground truth for states of a target
            W = model["sigma_v"] * model["B"].dot(np.random.randn(model["B"].shape[1], targetstate.shape[1])) #  equivalent to 'np.sqrt(model["Q"]).dot(np.random.randn(targetstate.shape[0],targetstate.shape[1])) '
            targetstate = model["F"].dot(targetstate) + W  # generate ground truth for states of a target
            truth_X.append(targetstate)

            # simulate target measurement
            # np.sqrt(model["R"]).dot(np.random.randn(model["R"].shape[1],1)) meas_zk = model["H"].dot(targetstate) + np.sqrt(model["R"]).dot(np.random.randn(model["R"].shape[1],1))  # generate measurement of a target
            V = model["D"].dot(np.random.randn(model["D"].shape[1], targetstate.shape[1])) # equivalent to 'np.sqrt(model["R"]).dot(np.random.randn(model["R"].shape[1],1)) '
            meas_zk = model["H"].dot(targetstate) + V # generate measurement of a target
            meas_Z.append(meas_zk)

            # Call KalmanFilter function
            m_update, P_update = KalmanFilter(meas_zk, model, m_update, P_update)
            estimated_X.append(m_update)


        return truth_X, meas_Z, estimated_X



    def KF_demo_plot(self, truth_X, meas_Z, estimated_X):
        truth_X = np.array(truth_X)
        meas_Z = np.array(meas_Z)
        estimated_X = np.array(estimated_X)

        # Plot the ground truth
        plt.plot(truth_X[:, 0], truth_X[:, 2], '.r', label='ground truth')
        # Plot the measurement
        plt.plot(meas_Z[:, 0], meas_Z[:, 1], '.b', label='measurement')
        # Plot the estimated state
        plt.plot(estimated_X[:, 0], estimated_X[:, 2], '.g', label='estimated state')
        plt.xlabel('X coordinate')
        plt.ylabel('Y coordinate')
        plt.legend()

        plt.show()


def main():
    Filter = KF_demo()
    truth_X, meas_Z, estimated_X = Filter.KF_demo_run()
    Filter.KF_demo_plot(truth_X, meas_Z, estimated_X)

if __name__ == "__main__":
    main()
