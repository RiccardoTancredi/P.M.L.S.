import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import PIL.Image as Image 
import matplotlib.image as img
from tqdm import tqdm


# read an image and convert it to a binary pattern os size [Lx, Ly]: 
# arguments = figurename, final size=[Lx, Ly]
def readPatterns(fname, size):
	this_img = Image.open(fname) # open colour image
	this_img = this_img.convert('1') # convert image to black and white
	this_img = this_img.resize(size=[size[0],size[1]]) # resize it to have the dimensions [Lx, Ly]
	this_img.save("%s_converted.png"%fname) # save converted image

	# plot the original and the converted file
	fig = plt.figure()

	# subplot on the left - original figure
	fig.add_subplot(1,2,1)		
	original = img.imread(fname)
	plt.imshow(original)
	plt.title("original")
	plt.axis('off')

	# subplot on the right - converted figure
	matrix = img.imread("%s_converted.png"%fname) # re-read from the file and convert it to a matrix of [0.,1.]
	matrix = 2*matrix-1 # convert from [0.,1.] -> [-1,1]
	fig.add_subplot(1,2,2)		
	plt.matshow(matrix, cmap=plt.cm.gray,fignum=0)
	plt.axis('off')
	plt.title("simplified")
	plt.show()
	
	# Return the pattern as a 1D vector (not as a 2D a matrix)
	return matrix.flatten()
	
#####################################################################################



class HopfieldNet:
	# initialize: 
	# arguments = number of neurons, list of patterns (vector of M components, each element of the pattern has to be an array of -1,+1 of size N)
	def __init__(self, N, patterns):
		self.N = N
		self.time_elapsed = 0.
		
		self.w = np.zeros([N,N]) # weights
		self.h = np.zeros(N) # threshold functions
		
		self.s = -np.ones(N) # default configuration = s[i]=-1
		
		# HEBBIAN RULE (h_i = 0., w_{ij} = sum_{k=1,...,M} s_i^k*s_j^k / M)
		print("The network is learning...")
		self.M = len(patterns)
		for k in range(self.M):
			print("pattern ", k)
			self.w += np.outer(patterns[k],patterns[k])/(1.*self.M)


		print("Done!")
	
	
#		# COMPUTE THE ENERGY - As before, I avoid loops and use efficient functions
		self.E = -0.5*np.sum(self.w) - np.sum(self.h) # energy for s_i = -1
		
		return

	# given and input s=[s_1,s_2,...,s_N], set the state of the network and recompute the energy
	def set_state(self, sinput):
		self.energy = []

		self.s = np.copy(sinput)

		# COMPUTE THE ENERGY - I use efficient functions rather than loops
		s2 = np.outer(self.s, self.s) # this returns a matrix s2[i,j]=s[i]*s[j]
		self.E = -0.5*np.sum(self.w*s2) + np.sum(self.h*self.s)
		self.energy.append(self.E)
	
		return

	# evolve the state of the networks doing a number "steps" of Monte Carlo steps
	def evolve(self, steps):
		for _ in tqdm(range(steps)):
			i = np.random.randint(self.N) # choose one node randomly
			
			sum_wijsj = np.sum(self.w[i,:]*self.s) # compute the argument of the activation function			
			if sum_wijsj < self.h[i]: # below the threshold
				self.s[i] = -1
			else: # above the threshold
				self.s[i] = 1

			s2 = np.outer(self.s, self.s) # this returns a matrix s2[i,j]=s[i]*s[j]
			self.energy.append(-0.5*np.sum(self.w*s2) + np.sum(self.h*self.s))
				
		return
#######################################################################
