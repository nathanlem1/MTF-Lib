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
import numpy as np
import copy
from scipy.stats import multivariate_normal
import math


def set_model():
    model = {}

    # Constant velocity (CV) model

    # Dynamic model parameters
    model['F_k'] = np.eye(4, dtype=np.float64)  # state transition model
    T = 1.0
    I = T*np.eye(2, dtype=np.float64)
    model['F_k'][0:2, 2:4] = I

    sigma_v = 5
    Q1 = np.array([[T ** 4 / 4, T ** 3 / 2], [T ** 3 / 2, T ** 2]], dtype=np.float64)
    Q = np.zeros((4, 4), dtype=np.float64)
    Q[np.ix_([0, 2], [0, 2])] = Q1
    Q[np.ix_([1, 3], [1, 3])] = Q1
    # cov of process noise
    model['Q_k'] = sigma_v ** 2 * Q

    # Observation model parameters
    model['H_k'] = np.array([[1., 0, 0, 0], [0, 1., 0, 0]], dtype=np.float64)  # observation model
    sigma_r = 6
    model['R_k'] = sigma_r ** 2 * np.eye(2,dtype=np.float64)  # the covariance of observation noise (change with the size of detection?)

    # Initial state covariance
    P_k = np.diag([100., 100., 25., 25.])
    model['P_k'] = np.array(P_k, dtype=np.float64)

    # Other important parameters
    model['w_birthsum'] = 0.0001 #0.02 # The total weight of birth targets. It is chosen depending on handling false positives.
    model['p_D'] = 0.98  # Probability of target detection,
    model['p_S'] = 0.99 # Probability of target survival (prob_death = 1 - prob_survival)
    model['T'] = 10**-5  # Pruning weight threshold.
    model['U'] = 4  # Merge distance threshold.
    model['w_thresh'] = 0.5 # State extraction weight threshold

    # Compute clutter intensity
    lambda_t = np.random.poisson(10) # Poisson average rate of uniform clutter (per scan); lambda_t = lambda_c*A
    x_range = [-1000, 1000]  # X range of measurements
    y_range = [-1000, 1000]  # Y range of measurements
    A = (x_range[1] - x_range[0])*(y_range[1]-y_range[0])
    pdf_c = 1.0/A
    clutterIntensity = lambda_t*pdf_c  # Generate clutter intensity.
    model['clutterIntensity'] = clutterIntensity
    model['lambda_t'] = lambda_t
    model['xrange'] = x_range
    model['yrange'] = y_range

    return model

# The probability density function (pdf) of the d-dimensional multivariate normal distribution
def mvnpdf(x, mean, covariance):
    # x = np.array(x, dtype=np.float64)
    # mean = np.array(mean, dtype=np.float64)
    # covariance = np.array(covariance, dtype=np.float64)

    d = mean.shape[0]
    delta_m = x - mean
    pdf_res = 1.0/(np.sqrt((2*np.pi)**d *np.linalg.det(covariance))) * np.exp(-0.5*np.transpose(delta_m).dot(np.linalg.inv(covariance)).dot(delta_m))[0][0]
    # pdf_res = 1.0 / (np.sqrt((2 * np.pi) ** d * np.linalg.det(covariance))) * math.exp(-0.5 * np.transpose(delta_m).dot(np.linalg.inv(covariance)).dot(delta_m))

    return pdf_res


class GM_PHD_Filter:
    def __init__(self, model):
        self.model = model #set_model

    # Step 1 and 2: Birth new targets and predict existing targets (Gaussian components) (According to the orignal paper)
    def Predict(self, Z_k, prunedIntensity):
        # An intensity (Probability Hypothesis Density - PHD) is described using weight, mean and covariance
        w = []  # weight of a Gaussian component
        m = []  # mean of a Gaussian component
        P = []  # Covariance of a gausssian component
        v_init = [0.0, 0.0] # initial velocity
        w_birthsum = self.model['w_birthsum']

        # Birth new targets
        for i in range(len(Z_k)):
            z_k = copy.deepcopy(Z_k[i])
            w.append(w_birthsum/len(Z_k))
            m.append(np.array([z_k[0], z_k[1], v_init[0], v_init[1]]).reshape(-1,1).astype('float64'))  # Targets are born here with [x, y, vx, vy] state format
            P.append(self.model['P_k'].astype('float64'))

        # Predict existing targets
        numTargets_Jk_minus_1 = len(prunedIntensity['w']) # Number of Gaussian components after the pruning and merging step
        for i in range(numTargets_Jk_minus_1):
            w.append(self.model['p_S']*prunedIntensity['w'][i])
            m.append(self.model['F_k'].dot(prunedIntensity['m'][i]).astype('float64'))
            P.append(self.model['Q_k'] + self.model['F_k'].dot(prunedIntensity['P'][i]).dot(np.transpose(self.model['F_k'])))

        predictedIntensity = {}
        predictedIntensity['w'] = w
        predictedIntensity['m'] = m
        predictedIntensity['P'] = P

        return predictedIntensity

    # Step 3 and 4: Construct PHD update components and doing observation update (According to the orignal paper)
    def Update(self, Z_k, predictedIntensity):

        # Construct PHD update components
        eta = []
        S = []
        K = []
        P = []
        for i in range(len(predictedIntensity['w'])):
            eta.append(self.model['H_k'].dot(predictedIntensity['m'][i]).astype('float64'))
            S.append(self.model['R_k'] + self.model['H_k'].dot(predictedIntensity['P'][i]).dot(np.transpose(self.model['H_k'])).astype('float64'))
            Si = copy.deepcopy(S[i])
            invSi = np.linalg.inv(np.array(Si, dtype=np.float64))  # Using normal inverse function
            # Vs = np.linalg.cholesky(np.array(Si, dtype=np.float64)); inv_sqrt_S = np.linalg.inv(Vs);
            # invSi = inv_sqrt_S.dot(np.transpose(inv_sqrt_S))  # Using Cholesky method
            K.append(predictedIntensity['P'][i].dot(np.transpose(self.model['H_k'])).dot(invSi).astype('float64'))
            P.append(predictedIntensity['P'][i] - K[i].dot(self.model['H_k']).dot(predictedIntensity['P'][i]).astype('float64'))

        constructUpdateIntensity = {}
        constructUpdateIntensity['eta'] = eta
        constructUpdateIntensity['S'] = S
        constructUpdateIntensity['K'] = K
        constructUpdateIntensity['P'] = P

        # Miss-detection part of GM-PHD update
        # We scale all weights by probability of missed detection (1 - p_D)
        w = []
        m = []
        P = []
        for i in range(len(predictedIntensity['w'])):
            w.append((1.0 - self.model['p_D'])*predictedIntensity['w'][i])
            m.append(predictedIntensity['m'][i])
            P.append(predictedIntensity['P'][i])

        # Detection part of GM-PHD update
        # Every observation updates every target
        numTargets_Jk_k_minus_1 = len(predictedIntensity['w']) # Number of Gaussian components after the prediction step
        l = 0
        for z in range(len(Z_k)):
            l = l + 1
            for j in range(len(predictedIntensity['w'])):
                z_k = copy.deepcopy(Z_k[z])
                # w.append(self.model['p_D'] * predictedIntensity['w'][j] * mvnpdf(z_k[0:2], constructUpdateIntensity['eta'][j][0:2], constructUpdateIntensity['S'][j][0:2, 0:2]))  # Hoping multivariate_normal.pdf is the right one to use; this is for video [x, y, w, h]
                w.append(self.model['p_D'] * predictedIntensity['w'][j] * mvnpdf(z_k, constructUpdateIntensity['eta'][j], constructUpdateIntensity['S'][j]))  # Hoping multivariate_normal.pdf is the right one to use; this is for simulation [x, y]
                m.append(predictedIntensity['m'][j] + constructUpdateIntensity['K'][j].dot(z_k - constructUpdateIntensity['eta'][j]).astype('float64'))
                P.append(constructUpdateIntensity['P'][j])

            total_w_d = 0.0
            for j in range(len(predictedIntensity['w'])):
                total_w_d += w[l*numTargets_Jk_k_minus_1 + j]

            for j in range(len(predictedIntensity['w'])):
                k_k = self.model['clutterIntensity']
                w[l*numTargets_Jk_k_minus_1 + j] = w[l*numTargets_Jk_k_minus_1 + j] / (k_k + total_w_d) # Updated weight

        # Combine both miss-detection and detection part of the GM-PHD update
        updatedIntensity = {}
        updatedIntensity['w'] = w
        updatedIntensity['m'] = m
        updatedIntensity['P'] = P

        return updatedIntensity

    # Step 5: Pruning and Merging
    def PruneAndMerge(self, updatedIntensity):
        w = []
        m = []
        P = []

        # Prune out the low-weighted components
        I = [index for index,value in enumerate(updatedIntensity['w']) if value > self.model['T']]  # Indices of large enough weights
        # I = np.where(np.array(updatedIntensity['w']) >= self.model['T'])[0]

        # Merge the close-together components
        while len(I) > 0:
            highWeights = np.array(updatedIntensity['w'])[I]
            j = np.argmax(highWeights)
            j = I[j]
            # Find all points with Mahalanobis distance less than U from point updatedIntensity['m'][j]
            L = []  # A vector of indices of merged Gaussians.
            for iterI in range(len(I)):
                thisI = copy.deepcopy(I[iterI])
                delta_m = updatedIntensity['m'][thisI] - updatedIntensity['m'][j]
                mahal_dist = np.transpose(delta_m).dot(np.linalg.inv(np.array(updatedIntensity['P'][thisI],dtype=np.float64))).dot(delta_m)
                if mahal_dist <= self.model['U']:
                    L.append(thisI)  # Indices of merged Gaussians

            # The new weight of the resulted merged Guassian is the summation of the weights of the Gaussian components.
            w_bar = sum(np.array(updatedIntensity['w'])[L])
            w.append(w_bar)

            # The new mean of the merged Gaussian is the weighted average of the merged means of Gaussian components.
            m_val = []
            for i in range(len(L)):
                thisI = copy.deepcopy(L[i])
                m_val.append(updatedIntensity['w'][thisI]*updatedIntensity['m'][thisI])
            m_bar = sum(m_val)/w_bar
            m.append(m_bar.astype('float64'))

            # Calculating covariance P_bar is a bit trickier
            P_val = []
            for i in range(len(L)):
                thisI = copy.deepcopy(L[i])
                delta_m = m_bar - updatedIntensity['m'][thisI]
                P_val.append(updatedIntensity['w'][thisI]*(updatedIntensity['P'][thisI] + delta_m.dot(np.transpose(delta_m))))
            P_bar = sum(P_val)/w_bar
            P.append(P_bar.astype('float64'))

            # Now delete the elements in L from I
            for i in L:
                I.remove(i)

        prunedMergedIntensity = {}
        prunedMergedIntensity['w'] = w
        prunedMergedIntensity['m'] = m
        prunedMergedIntensity['P'] = P

        return prunedMergedIntensity


    # Step 6: extracting estimated states
    def ExtractStates(self, PrunedAndMerged):
        w = []
        m = []
        P = []
        # PrunedAndMerged = self.PruneAndMerge()
        for i in range(len(PrunedAndMerged['w'])):
            if PrunedAndMerged['w'][i] > self.model['w_thresh']:
                for j in range(int(round(PrunedAndMerged['w'][i]))):  # If a target has a rounded weight greater than 1, output it multiple times.
                   w.append(PrunedAndMerged['w'][i])
                   m.append(PrunedAndMerged['m'][i])
                   P.append(PrunedAndMerged['P'][i])

        extractedStates = {}
        extractedStates['w'] = w
        extractedStates['m'] = m
        extractedStates['P'] = P

        return extractedStates





