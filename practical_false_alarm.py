import matplotlib.pyplot as plt
import numpy as np
import periodicity2



for sigma in [.001, .01, .1, 1., 2., 3.]:

    model_data = np.random.normal(0, sigma, 60000)

