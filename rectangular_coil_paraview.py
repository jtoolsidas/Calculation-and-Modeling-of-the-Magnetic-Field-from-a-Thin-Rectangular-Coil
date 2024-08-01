#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:41:40 2024

@author: akaptanoglu
"""

import numpy as np
from pyevtk.hl import polyLinesToVTK, pointsToVTK

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

a = 0.7071067811865477
mu0 = np.pi * 4 * 1e-7
x = np.linspace(-1, 1, 5)
xv, yv, zv = np.meshgrid(x, x, x)
x = np.ravel(xv)
y = np.ravel(yv)
z = np.ravel(zv)
def mu(a):
    return (1/np.sqrt(2*a**(-2)))
def b(y):
    return (mu0 * np.sqrt(2*a**(-2)))/np.pi
def r1(x,y,z,a): 
    return np.sqrt((x+a)**2 + (y+a)**2 + z**2)
def r2(x,y,z,a): 
    return np.sqrt((x-a)**2 + (y+a)**2 + z**2)
def r3(x,y,z,a): 
    return np.sqrt((x+a)**2 + (y-a)**2 + z**2)
def r4(x,y,z,a): 
    return np.sqrt((x-a)**2 + (y-a)**2 + z**2)

bx = ((-mu(a) * b(a) * z)/4) * (
    (1/(r1(x,y,z,a) * (r1(x,y,z,a) - y - a))) - (
    1/(r2(x,y,z,a) * (r2(x,y,z,a) - y - a))) - (
    1/(r3(x,y,z,a) * (r3(x,y,z,a) - y + a))) + (
    1/(r4(x,y,z,a) * (r4(x,y,z,a) - y + a))))
by = ((-mu(a) * b(a) * z)/4) * (
    (1/(r1(x,y,z,a) * (r1(x,y,z,a) - x - a))) - (
    1/(r2(x,y,z,a) * (r2(x,y,z,a) - x + a))) - (
    1/(r3(x,y,z,a) * (r3(x,y,z,a) - x - a))) + (
    1/(r4(x,y,z,a) * (r4(x,y,z,a) - x + a))))
bz = ((mu(a) * b(a))/4) * (
    ((x + a)/(r1(x,y,z,a) * (r1(x,y,z,a) - y - a)) + \
    (y + a)/(r1(x,y,z,a) * (r1(x,y,z,a) - x - a)) \
    -(x - a)/(r2(x,y,z,a) * (r2(x,y,z,a) - y - a)) \
    -(y + a)/(r2(x,y,z,a) * (r2(x,y,z,a) - x + a)) \
    -(x + a)/(r3(x,y,z,a) * (r3(x,y,z,a) - y + a)) \
    -(y - a)/(r3(x,y,z,a) * (r3(x,y,z,a) - x - a)) \
    +(x - a)/(r4(x,y,z,a) * (r4(x,y,z,a) - y + a)) \
    +(y - a)/(r4(x,y,z,a) * (r4(x,y,z,a) - x + a))))
print(bx, by, bz)

data = {"B": (bx, by, bz)}
pointsToVTK(filename + '_B', X, Y, Z, data=data)
