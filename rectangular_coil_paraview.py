#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:41:40 2024

@author: akaptanoglu
"""

import numpy as np
from pyevtk.hl import polyLinesToVTK
filename = 'paraview_coil'
N = 1000
I = np.ones(N)

def cot(y):
  if not np.all(np.isclose(np.tan(y), 0.0)):
    return 1/np.tan(y)
  else:
      return 0

def floor(x): 
    return (-0.5 + x + np.arctan(cot(np.pi * x))/np.pi)

def l(phi):
    return (((np.cos(np.pi/4))/(np.cos(phi- (np.pi/2) * floor((4*phi + np.pi)/(2*np.pi))))) * np.array([np.cos(phi), np.sin(phi), np.zeros(len(phi))]))

phi = np.linspace(0, 2 * np.pi, N)
X, Y, Z = l(phi)
data = np.concatenate([I[i] * np.ones((N, )) for i in range(N)])
polyLinesToVTK(filename, X, Y, Z, pointsPerLine=np.asarray([N]), pointData={'idx': data})
