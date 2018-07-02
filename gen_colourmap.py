#!/usr/bin/python3
# Imports
import numpy as np

# number of discrete colours
N = 5

# generate R array
R = np.zeros(N).astype(np.uint8)
R[:int(N/4)] = 255
R[int(N/4):int(2*N/4)+1] = np.linspace(255,0,num=(N/4)+1,endpoint=True)

# generate G array
G = 255*np.ones(N).astype(np.uint8)
G[0:int(N/4)+1] = np.linspace(0,255,num=(N/4)+1,endpoint=True)
G[int(3*N/4):] = np.linspace(255,0,num=(N/4)+1,endpoint=True)

# generate B array
B = np.zeros(N).astype(np.uint8)
B[int(2*N/4):int(3*N/4)+1] = np.linspace(0,255,num=(N/4)+1,endpoint=True)
B[int(3*N/4)+1:] = 255

# stack arrays
RGB = np.dstack((R,G,B,255*np.ones(N).astype(np.uint8)))[0]
print(RGB)
