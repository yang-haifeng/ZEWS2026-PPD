import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("image.out", skiprows = 6)
data.shape = 100, 100

plt.imshow(data)
plt.show()
