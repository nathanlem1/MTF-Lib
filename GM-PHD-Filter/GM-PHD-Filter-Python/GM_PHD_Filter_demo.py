from GM_PHD_Filter import set_model, GM_PHD_Filter
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
plt.ion()
fig = plt.figure()

model = set_model()

# Starting target states (positions , velocities)
m_start1 = np.array([250, 250, 0, 0]).reshape(4,1)
m_start2 = np.array([10,-10, 5, -5]).reshape(4,1)
m_start3 = np.array([-250, -250, 0, 0]).reshape(4,1)


def gen_newstates(model, targetStates):
    truthStates = []
    for i in range(len(targetStates)):
        W = np.sqrt(model['Q_k']).dot(np.random.randn(model['Q_k'].shape[0], targetStates[i].shape[1]))#.reshape(-1,1) # Process noise
        truthState = model['F_k'].dot(targetStates[i]) + W   # New target true state
        truthStates.append(truthState)

    return truthStates

def gen_observations(model, truthStates):
    obserations = []
    # We are not guaranteed to detect the target - there is only a probability
    for i in range(len(truthStates)):
        detect_i = np.random.rand() # Uniformly distribution.
        if detect_i <= model['p_D']:
            V = np.sqrt(model['R_k']).dot(np.random.randn(model['R_k'].shape[0], truthStates[i].shape[1]))#.reshape(-1,1) # Observation noise
            obs = model['H_k'].dot(truthStates[i]) + V
            obserations.append(obs)

    return obserations

def gen_clutter(model):
    lambda_t = model['lambda_t']
    x_range = model['xrange']
    y_range = model['yrange']
    clutter = []
    for n in range(lambda_t):
        clutterX =  np.random.rand() * (x_range[1] - x_range[0]) + x_range[0] # Random number between x_range[0] and x_range[1], uniformly distributed.
        clutterY = np.random.rand() * (y_range[1] - y_range[0]) + y_range[0]
        clutter.append(np.array([clutterX, clutterY]).reshape(-1,1))
    return clutter

def GMPHD_Filter_demo_plot(truthStates, observations, estimatedStates, clutter):

    # plt.clf()  # Comment this line for cummulative show or uncomment it for frame-wise show

    # Plot the ground truth
    for i in range(len(truthStates)):
        truthState = truthStates[i]
        # plt.plot(truthState[0], truthState[1], '.b', markersize = 10.0, label='ground truth')
        plt.plot(truthState[0], truthState[1], '.b', markersize = 10.0)

    # Plot the measurements
    for i in range(len(observations)):
        observation = observations[i]
        # plt.plot(observation[0], observation[1], '.r', markersize = 10.0, label='measurement')
        if len(observation) > 0:
            plt.plot(observation[0], observation[1], '.r', markersize = 10.0)

    # Plot the clutter
    for i in range(len(clutter)):
        clut = clutter[i]
        # plt.plot(observation[0], observation[1], 'xk', markersize = 5.0, label='clutter')
        plt.plot(clut[0], clut[1], 'xk', markersize = 5.0)

    # Plot the estimated state
    for i in range(len(estimatedStates)):
        estimatedState = np.array(estimatedStates[i], dtype=np.float64)
        # plt.plot(estimatedState[0], estimatedState[1], '.g', markersize = 10.0, label='estimated state')
        plt.plot(estimatedState[0], estimatedState[1], '.g', markersize = 10.0)

    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.title('Ground truth (blue), true observations (red), estimated states(green) and clutter (black x)', fontsize=8)
    plt.xlim((model['xrange'][0], model['xrange'][1]))
    plt.ylim((model['yrange'][0], model['yrange'][1]))
    plt.legend()

    fig.canvas.draw()

# Initialize the initial pruned intensity
prunedIntensity = {}
prunedIntensity['w'] = []
prunedIntensity['m'] = []
prunedIntensity['P'] = []

targetStates = [m_start1, m_start2, m_start3]

n_scan = 120 # Duration of iterations (simulation)

# For analysis
prunedIntensity_list = []
estimates_list = []
targetStates_list  = []
observations_list = []

for i in range(n_scan):
    print(i)
    # Generate target states, actual observations, and clutter
    targetStates = gen_newstates(model,targetStates)
    observations = gen_observations(model, targetStates)
    clutter = gen_clutter(model)
    Z_k = observations + clutter   # Add clutter to the observations to mimic cluttered environment

    # Apply GM-PHD-Filter
    Filter = GM_PHD_Filter(model)
    predictedIntensity = Filter.Predict(Z_k, prunedIntensity)
    updatedIntensity = Filter.Update(Z_k, predictedIntensity)
    prunedIntensity = Filter.PruneAndMerge(updatedIntensity)
    estimates = Filter.ExtractStates(prunedIntensity)  # extracting estimates from the pruned intensity this gives better result than extracting them from the updated intensity!
    estimatedStates = estimates['m']

    # Plot ground-truth states, true observations, estiamted states and clutter
    GMPHD_Filter_demo_plot(targetStates, observations, estimatedStates, clutter)

    # Store output intensity and estimates for analysis
    prunedIntensity_list.append(prunedIntensity)
    estimates_list.append(estimates)
    targetStates_list.append(targetStates)
    observations_list.append(Z_k) # observations




