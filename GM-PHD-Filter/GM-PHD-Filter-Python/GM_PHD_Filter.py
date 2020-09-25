"""
This is the GM-PHD-Filter implementation of the algorithm in [1] with the assumption of measurement-driven birth of
targets which has been used in work [2]. In work [2] i.e. in visual tracking applications, extracting states from the
pruned intensity gives better result than extracting them from the updated intensity.

References:
[1]. B.-N. Vo, W.-K. Ma, "The Gaussian Mixture Probability Hypothesis Density Filter", IEEE Transactions on Signal
Processing, Vol 54, No. 11, November 2006, pp4091-4104
[2]. Nathanael L. Baisa, "Online multi-object visual tracking using a GM-PHD filter with deep appearance learning",
22nd International Conference on Information Fusion (FUSION), July, 2019.

# -------------------------------------------------------------------------
 Nathanael L. Baisa: nathanaellmss@gmail.com
 Original:
 Modified: July, 2019
# -------------------------------------------------------------------------
"""

import copy
import math

import numpy as np


def set_model():

    model = dict()

    # Constant velocity (CV) model - Dynamic model parameters
    model['F_k'] = np.eye(4, dtype=np.float64)  # state transition model
    T = 1.0
    I = T * np.eye(2, dtype=np.float64)
    model['F_k'][0:2, 2:4] = I

    sigma_v = 5
    Q1 = np.array([[T ** 4 / 4, T ** 3 / 2], [T ** 3 / 2, T ** 2]], dtype=np.float64)
    Q = np.zeros((4, 4), dtype=np.float64)
    Q[np.ix_([0, 2], [0, 2])] = Q1
    Q[np.ix_([1, 3], [1, 3])] = Q1
    # cov of process noise
    model['Q_k'] = sigma_v ** 2 * Q

    # Constant velocity (CV) model - Observation model parameters
    model['H_k'] = np.array([[1., 0, 0, 0], [0, 1., 0, 0]], dtype=np.float64)  # observation model
    sigma_r = 6
    model['R_k'] = sigma_r ** 2 * np.eye(2, dtype=np.float64)  # the covariance of observation noise

    # Initial state covariance
    P_k = np.diag([100., 100., 25., 25.])
    model['P_k'] = np.array(P_k, dtype=np.float64)

    # Other important parameters
    model['w_birth_sum'] = 0.0001  # 0.02 # The total weight of birth targets. It is chosen depending on handling false
    # positives.

    model['p_D'] = 0.98  # Probability of target detection,
    model['p_S'] = 0.99  # Probability of target survival (prob_death = 1 - prob_survival)
    model['T'] = 10**-5  # Pruning weight threshold.
    model['U'] = 4  # Merge distance threshold.
    model['w_thresh'] = 0.5  # State extraction weight threshold

    # Compute clutter intensity
    lambda_t = np.random.poisson(10)  # Poisson average rate of uniform clutter (per scan); lambda_t = lambda_c*A
    x_range = [-1000, 1000]  # X range of measurements
    y_range = [-1000, 1000]  # Y range of measurements
    A = (x_range[1] - x_range[0])*(y_range[1]-y_range[0])
    pdf_c = 1.0/A
    clutter_intensity = lambda_t*pdf_c  # Generate clutter intensity.
    model['clutterIntensity'] = clutter_intensity
    model['lambda_t'] = lambda_t
    model['x_range'] = x_range
    model['y_range'] = y_range

    return model


# The probability density function (pdf) of the d-dimensional multivariate normal distribution
def mvn_pdf(x, mean, covariance):

    # x = np.array(x, dtype=np.float64)
    # mean = np.array(mean, dtype=np.float64)
    # covariance = np.array(covariance, dtype=np.float64)

    d = mean.shape[0]
    delta_m = x - mean
    pdf_res = 1.0/(np.sqrt((2*np.pi)**d * np.linalg.det(covariance))) * np.exp(-0.5*np.transpose(delta_m).dot(
        np.linalg.inv(covariance)).dot(delta_m))[0][0]
    # pdf_res = 1.0 / (np.sqrt((2 * np.pi) ** d * np.linalg.det(covariance))) * \
    #           math.exp(-0.5 * np.transpose(delta_m).dot(np.linalg.inv(covariance)).dot(delta_m))

    return pdf_res


# Class implementing a GM-PHD Filter
class GM_PHD_Filter:
    """
    This class implements the Gaussian mixture probability hypothesis density (GM-PHD) filter.
    """
    def __init__(self, model):
        self.model = model  # set_model

    # Step 1 and 2: Birth new targets and predict existing targets (Gaussian components)
    # (According to the original paper)
    def predict(self, Z_k, pruned_intensity):

        # An intensity (Probability Hypothesis Density - PHD) is described using weight, mean and covariance
        w = []  # weight of a Gaussian component
        m = []  # mean of a Gaussian component
        P = []  # Covariance of a Gaussian component
        v_init = [0.0, 0.0]  # initial velocity
        w_birth_sum = self.model['w_birth_sum']

        # Birth new targets
        for i in range(len(Z_k)):
            z_k = copy.deepcopy(Z_k[i])
            w.append(w_birth_sum/len(Z_k))
            m.append(np.array([z_k[0][0], z_k[1][0], v_init[0], v_init[1]]).reshape(-1, 1).astype('float64'))  # Targets
            # are born here with [x, y, vx, vy] state format

            P.append(self.model['P_k'].astype('float64'))

        # Predict existing targets
        num_targets_Jk_minus_1 = len(pruned_intensity['w'])  # Number of Gaussian components after the pruning and
        # merging step

        for i in range(num_targets_Jk_minus_1):
            w.append(self.model['p_S']*pruned_intensity['w'][i])
            m.append(self.model['F_k'].dot(pruned_intensity['m'][i]).astype('float64'))
            P.append(self.model['Q_k'] +
                     self.model['F_k'].dot(pruned_intensity['P'][i]).dot(np.transpose(self.model['F_k'])))

        predicted_intensity = dict()
        predicted_intensity['w'] = w
        predicted_intensity['m'] = m
        predicted_intensity['P'] = P

        return predicted_intensity

    # Step 3 and 4: Construct PHD update components and doing observation update (According to the original paper)
    def update(self, Z_k, predicted_intensity):

        # Construct PHD update components
        eta = []
        S = []
        K = []
        P = []
        for i in range(len(predicted_intensity['w'])):
            eta.append(self.model['H_k'].dot(predicted_intensity['m'][i]).astype('float64'))
            hph_transpose = self.model['H_k'].dot(predicted_intensity['P'][i]).dot(np.transpose(self.model['H_k']))
            S.append(self.model['R_k'] + hph_transpose.astype('float64'))
            Si = copy.deepcopy(S[i])
            invSi = np.linalg.inv(np.array(Si, dtype=np.float64))  # Using normal inverse function
            # Vs = np.linalg.cholesky(np.array(Si, dtype=np.float64)); inv_sqrt_S = np.linalg.inv(Vs);
            # invSi = inv_sqrt_S.dot(np.transpose(inv_sqrt_S))  # Using Cholesky method
            K.append(predicted_intensity['P'][i].dot(np.transpose(self.model['H_k'])).dot(invSi).astype('float64'))
            P.append(predicted_intensity['P'][i] -
                     K[i].dot(self.model['H_k']).dot(predicted_intensity['P'][i]).astype('float64'))

        construct_update_intensity = dict()
        construct_update_intensity['eta'] = eta
        construct_update_intensity['S'] = S
        construct_update_intensity['K'] = K
        construct_update_intensity['P'] = P

        # Miss-detection part of GM-PHD update
        # We scale all weights by probability of missed detection (1 - p_D)
        w = []
        m = []
        P = []
        for i in range(len(predicted_intensity['w'])):
            w.append((1.0 - self.model['p_D'])*predicted_intensity['w'][i])
            m.append(predicted_intensity['m'][i])
            P.append(predicted_intensity['P'][i])

        # Detection part of GM-PHD update
        # Every observation updates every target
        num_targets_Jk_k_minus_1 = len(predicted_intensity['w'])  # Number of Gaussian components after the prediction
        # step

        l = 0
        for z in range(len(Z_k)):
            l = l + 1
            for j in range(num_targets_Jk_k_minus_1):
                z_k = copy.deepcopy(Z_k[z])
                # w.append(self.model['p_D'] * predicted_intensity['w'][j] *
                #          mvn_pdf(z_k[0:2],
                #                  construct_update_intensity['eta'][j][0:2],
                #                  construct_update_intensity['S'][j][0:2, 0:2]))  # Hoping multivariate_normal.pdf is
                # # the right one to use; this is for video [x, y, w, h]

                w.append(self.model['p_D'] * predicted_intensity['w'][j] *
                         mvn_pdf(z_k, construct_update_intensity['eta'][j], construct_update_intensity['S'][j]))  #
                # Hoping multivariate_normal.pdf is the right one to use; this is for simulation [x, y]

                m.append(predicted_intensity['m'][j] +
                         construct_update_intensity['K'][j].dot(z_k -
                                                                construct_update_intensity['eta'][j]).astype('float64'))
                P.append(construct_update_intensity['P'][j])

            total_w_d = 0.0
            for j in range(num_targets_Jk_k_minus_1):
                total_w_d += w[l*num_targets_Jk_k_minus_1 + j]

            for j in range(num_targets_Jk_k_minus_1):
                k_k = self.model['clutterIntensity']
                w[l*num_targets_Jk_k_minus_1 + j] = w[l*num_targets_Jk_k_minus_1 + j] / (k_k + total_w_d)  # Updated
                # weight

        # Combine both miss-detection and detection part of the GM-PHD update
        updated_intensity = dict()
        updated_intensity['w'] = w
        updated_intensity['m'] = m
        updated_intensity['P'] = P

        return updated_intensity

    # Step 5: Pruning and Merging
    def prune_and_merge(self, updated_intensity):
        w = []
        m = []
        P = []

        # Prune out the low-weighted components
        I = [index for index, value in enumerate(updated_intensity['w']) if value > self.model['T']]  # Indices of large
        # enough weights

        # I = np.where(np.array(updated_intensity['w']) >= self.model['T'])[0]

        # Merge the close-together components
        while len(I) > 0:
            high_weights = np.array(updated_intensity['w'])[I]
            j = np.argmax(high_weights)
            j = I[j]
            # Find all points with Mahalanobis distance less than U from point updated_intensity['m'][j]
            L = []  # A vector of indices of merged Gaussian components.
            for iterI in range(len(I)):
                this_I = copy.deepcopy(I[iterI])
                delta_m = updated_intensity['m'][this_I] - updated_intensity['m'][j]
                mahal_dist = np.transpose(delta_m).dot(np.linalg.inv(
                    np.array(updated_intensity['P'][this_I], dtype=np.float64))).dot(delta_m)
                if mahal_dist <= self.model['U']:
                    L.append(this_I)  # Indices of merged Gaussian components

            # The new weight of the resulted merged Gaussian is the summation of the weights of the Gaussian components.
            w_bar = sum(np.array(updated_intensity['w'])[L])
            w.append(w_bar)

            # The new mean of the merged Gaussian is the weighted average of the merged means of Gaussian components.
            m_val = []
            for i in range(len(L)):
                this_I = copy.deepcopy(L[i])
                m_val.append(updated_intensity['w'][this_I]*updated_intensity['m'][this_I])
            m_bar = sum(m_val)/w_bar
            m.append(m_bar.astype('float64'))

            # Calculating covariance P_bar is a bit trickier
            P_val = []
            for i in range(len(L)):
                this_I = copy.deepcopy(L[i])
                delta_m = m_bar - updated_intensity['m'][this_I]
                P_val.append(updated_intensity['w'][this_I]*(updated_intensity['P'][this_I] +
                                                             delta_m.dot(np.transpose(delta_m))))
            P_bar = sum(P_val)/w_bar
            P.append(P_bar.astype('float64'))

            # Now delete the elements in L from I
            for i in L:
                I.remove(i)

        pruned_merged_intensity = dict()
        pruned_merged_intensity['w'] = w
        pruned_merged_intensity['m'] = m
        pruned_merged_intensity['P'] = P

        return pruned_merged_intensity

    # Step 6: extracting estimated states
    def extract_states(self, pruned_and_merged):
        w = []
        m = []
        P = []
        # pruned_and_merged = self.prune_and_merge()
        for i in range(len(pruned_and_merged['w'])):
            if pruned_and_merged['w'][i] > self.model['w_thresh']:
                for j in range(int(round(pruned_and_merged['w'][i]))):  # If a target has a rounded weight greater
                    # than 1, output it multiple times.

                    w.append(pruned_and_merged['w'][i])
                    m.append(pruned_and_merged['m'][i])
                    P.append(pruned_and_merged['P'][i])

        extracted_states = dict()
        extracted_states['w'] = w
        extracted_states['m'] = m
        extracted_states['P'] = P

        return extracted_states
