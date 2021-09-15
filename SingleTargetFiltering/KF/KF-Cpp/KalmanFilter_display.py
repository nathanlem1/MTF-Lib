import numpy as np
import matplotlib.pyplot as plt
import pprint

groundTruthPositions = []
estimatedPositions = []

with open("output.txt") as f:
    contents = f.readlines()
    pprint.pprint(contents)

    for i in range(len(contents)):
        row = contents[i].split()
        rowf = [float(s) for s in row]
        groundTruthPositions.append(rowf[:2])
        estimatedPositions.append(rowf[2:])

groundTruthPositions = np.array(groundTruthPositions)
estimatedPositions = np.array(estimatedPositions)

# Plot the ground truth
plt.plot(groundTruthPositions[:, 0], groundTruthPositions[:, 1], '-r', label='ground truth')
# Plot the estimated state
plt.plot(estimatedPositions[:, 0], estimatedPositions[:, 1], '.g', label='estimated state')
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.legend()

plt.show()


print('ok')
