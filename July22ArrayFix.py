#!/usr/bin/env python3
"""
Test file for 3D FE class
"""
import numpy as np
phi = np.linspace(0, 2*np.pi, 1000)
print(phi.size)

#print(phi)

B = np.array([[np.cos(phi), np.sin(phi), np.zeros(1000)]])


A = np.zeros((3,1000))
A[0,:] = np.cos(phi)
A[1,:] = np.sin(phi)
A[2,:] = np.zeros((1000))

#print(np.transpose(A))

B = np.transpose(A)
C = (((np.cos(np.pi/4))/(np.cos(phi- (np.pi/2) * ((4*phi + np.pi)//(2*np.pi))))))


l = np.dot((((np.cos(np.pi/4))/(np.cos(phi- (np.pi/2) * ((4*phi + np.pi)//(2*np.pi)))))),B)
l = np.dot(C,B)
print(l)


# B = ((np.cos(np.pi/4))/(np.cos(phi- (np.pi/2) * ((4*phi + np.pi)//(2*np.pi))))).shape

# print(B)


# #print(np.array([[np.cos(phi), np.sin(phi),0
# A = np.array([[np.cos(phi), np.sin(phi), np.zeros(1000)]])
# print(np.array([[np.cos(phi), np.sin(phi)]]))
# print(np.shape(A))


# #A = np.reshape()


# #print(np.cos(phi))