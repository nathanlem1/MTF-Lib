"""
This a demo of GM-PHD filter simulation for filtering multiple targets.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from GM_PHD_Filter import set_model, GM_PHD_Filter

matplotlib.use('tkagg')
plt.ion()
fig = plt.figure()

model = set_model()

# Starting target states (positions , velocities)
m_start1 = np.array([250, 250, 0, 0]).reshape(4, 1)
m_start2 = np.array([10, -10, 5, -5]).reshape(4, 1)
m_start3 = np.array([-250, -250, 0, 0]).reshape(4, 1)


def gen_new_states(model, target_states):
    truth_states = []
    for i in range(len(target_states)):
        # Process noise
        W = np.sqrt(model['Q_k']).dot(np.random.randn(model['Q_k'].shape[0], target_states[i].shape[1]))
        truth_state = model['F_k'].dot(target_states[i]) + W   # New target true state
        truth_states.append(truth_state)

    return truth_states


def gen_observations(model, truth_states):
    obserations = []

    # We are not guaranteed to detect the target - there is only a probability
    for i in range(len(truth_states)):
        detect_i = np.random.rand()  # Uniformly distribution.
        if detect_i <= model['p_D']:
            # Observation noise
            V = np.sqrt(model['R_k']).dot(np.random.randn(model['R_k'].shape[0], truth_states[i].shape[1]))
            obs = model['H_k'].dot(truth_states[i]) + V
            obserations.append(obs)

    return obserations


def gen_clutter(model):
    lambda_t = model['lambda_t']
    x_range = model['x_range']
    y_range = model['y_range']
    clutter = []
    for n in range(lambda_t):
        clutter_X = np.random.rand() * (x_range[1] - x_range[0]) + x_range[0]  # Random number between x_range[0] and
        # x_range[1], uniformly distributed.

        clutter_Y = np.random.rand() * (y_range[1] - y_range[0]) + y_range[0]
        clutter.append(np.array([clutter_X, clutter_Y]).reshape(-1, 1))

    return clutter


def demo_plot(truth_states, observations, estimated_states, clutter):

    # plt.clf()  # Comment this line for cumulative show or uncomment it for frame-wise show

    # Plot the ground truth
    for i in range(len(truth_states)):
        truth_state = truth_states[i]
        # plt.plot(truth_state[0], truth_state[1], '.b', markersize=10.0, label='ground truth')
        plt.plot(truth_state[0], truth_state[1], '.b', markersize=10.0)

    # Plot the measurements
    for i in range(len(observations)):
        observation = observations[i]
        # plt.plot(observation[0], observation[1], '.r', markersize=10.0, label='measurement')
        if len(observation) > 0:
            plt.plot(observation[0], observation[1], '.r', markersize=10.0)

    # Plot the clutter
    for i in range(len(clutter)):
        clutter_i = clutter[i]
        # plt.plot(observation[0], observation[1], 'xk', markersize=5.0, label='clutter')
        plt.plot(clutter_i[0], clutter_i[1], 'xk', markersize=5.0)

    # Plot the estimated state
    for i in range(len(estimated_states)):
        estimated_state = np.array(estimated_states[i], dtype=np.float64)
        # plt.plot(estimatedState[0], estimatedState[1], '.g', markersize=10.0, label='estimated state')
        plt.plot(estimated_state[0], estimated_state[1], '.g', markersize=10.0)

    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.title('Ground truth (blue), true observations (red), estimated states(green) and clutter (black x)', fontsize=8)
    plt.xlim((model['x_range'][0], model['x_range'][1]))
    plt.ylim((model['y_range'][0], model['y_range'][1]))
    # plt.legend()

    fig.canvas.draw()


# Initialize the initial pruned intensity
pruned_intensity = dict()
pruned_intensity['w'] = []
pruned_intensity['m'] = []
pruned_intensity['P'] = []

target_states = [m_start1, m_start2, m_start3]

n_scan = 120  # Duration of iterations (simulation)

# For analysis
pruned_intensity_list = []
estimates_list = []
target_states_list = []
observations_list = []

# Apply the GM-PHD filter for filtering of multiple targets (simulation)
for i in range(n_scan):
    print(i)

    # Generate target states, actual observations, and clutter
    target_states = gen_new_states(model, target_states)
    observations = gen_observations(model, target_states)
    clutter = gen_clutter(model)
    Z_k = observations + clutter   # Add clutter to the observations to mimic cluttered environment

    # Apply GM-PHD Filter
    Filter = GM_PHD_Filter(model)
    predicted_intensity = Filter.predict(Z_k, pruned_intensity)
    updated_intensity = Filter.update(Z_k, predicted_intensity)
    pruned_intensity = Filter.prune_and_merge(updated_intensity)
    estimates = Filter.extract_states(pruned_intensity)  # extracting estimates from the pruned intensity this gives
    # better result than extracting them from the updated intensity!

    estimated_states = estimates['m']

    # Plot ground-truth states, true observations, estimated states and clutter
    demo_plot(target_states, observations, estimated_states, clutter)

    # Store output intensity and estimates for analysis
    pruned_intensity_list.append(pruned_intensity)
    estimates_list.append(estimates)
    target_states_list.append(target_states)
    observations_list.append(Z_k)  # observations
