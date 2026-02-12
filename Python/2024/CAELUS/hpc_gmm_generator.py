# Initializing GMM Class and Imports
import numpy as np
from matplotlib import pyplot as plt # Utilizing NVIDIA CUDA through Numba just-in-time compiler
from numba import jit, cuda, njit, prange
from scipy.stats import multivariate_normal as mvn
from sklearn.datasets import make_spd_matrix
import math
import time

class GMM:

    def __init__(self, ws, mus, covs):
      """
      Create GMM object.
      """
      self.ws = ws
      self.mus = mus
      self.covs = covs 
      self.pdfObjs = self.__getObjs()

    def __getObjs(self):
        objs = []
        for ii in range(len(self.ws)):
            objs.append(mvn(mean= self.mus[:, ii], cov = self.covs[:, :, ii]))
        return objs


    def __sample(self):
        acc_ws = [np.sum(self.ws[:i]) for i in range(1, len(self.ws)+1)]
        assert np.isclose(acc_ws[-1], 1)
        r = np.random.uniform(0, 1)
        k = 0
        for i, threshold in enumerate(acc_ws):
            if r < threshold:
                k = i
                break
        return self.pdfObjs[k].rvs()

    def getSamples(self, N):
      """
      Returns random sample of the pdf.
      """
      return np.array([self.__sample() for _ in range(N)]).T

# Create GMM Object
nComp  = 10
nPixels = 100
maxPixelValue = 10

ws_unorm = np.random.uniform(0.0, 1.0, size = nComp)
ws = ws_unorm/np.sum(ws_unorm)

means_x = np.random.uniform(-maxPixelValue, maxPixelValue, size = nComp)
means_y = np.random.uniform(-maxPixelValue, maxPixelValue, size = nComp)
means = np.vstack((means_x, means_y))

covs = np.zeros((2, 2, nComp))
for ii in range(nComp):
    spd = make_spd_matrix(n_dim=2)
    covs[:, :, ii] = spd

# Generate the GMM object
gmm = GMM(ws, means, covs)

# Setting Objects
objs = gmm.pdfObjs
ws = gmm.ws
mus = gmm.mus
covs = gmm.covs

# Initialize grid
x_values = np.linspace(-maxPixelValue, maxPixelValue, nPixels)
y_values = np.linspace(-maxPixelValue, maxPixelValue, nPixels)
X, Y = np.meshgrid(x_values, y_values)

xFlat = X.flatten()
yFlat = Y.flatten()

xs = np.vstack((xFlat, yFlat))

start = time.time()

# Creating PDF Creation Function
@njit
def create_pdf(weights, means, covar, grid):

  """
  Input
    wights: matrix of weights as initialized above
    means: matrix of means as initialized above
    covar: covariance matrix as initialized above
    grid: x-y values vertically stacked as above

  Output: Outputs the full pdf grid in "xs" format
  """

  pdfarray = np.zeros(grid.shape[1])
  for ii in range(grid.shape[1]):
    # Private PDF Function
    pdf = 0.0
    for jj in range(len(weights)):
      w = weights[jj]

      """
      obj = mvn(mus[:, jj], covs[:, :, jj])
      y = obj.pdf(grid[:, ii])
      """

      # Implementing mvn function
      k = len(means)
      xx = grid[:, ii]
      uu = means[:, jj]
      cov = covar[:, :, jj]

      t1 = (2 * np.pi)**k          
      t2 = np.linalg.det(cov)      
      t3 = 1.0 / np.sqrt(t1 * t2)  

      t4 = (xx - uu).T          
      t5 = np.linalg.inv(cov)   
      t6 = (xx - uu)            
      t7 = -0.5 * (np.dot(t4,t5).dot(t6))
      result = t3 * np.exp(t7) 
      y = result
      
      pdf = pdf + w*y
    # Public PDF Function
    pdfarray[ii] = pdf
    
  # Output the correct format of the array in 2D form
  return pdfarray

# Calling the parallel optimization function
pdf_vals_flat = create_pdf(ws, mus, covs, xs)
pdf_vals_grid = np.reshape(pdf_vals_flat, X.shape)

end = time.time()
print(f"The runtime was {end-start}")

# Save Data
np.savetxt('data.csv', pdf_vals_grid, delimiter=',')
