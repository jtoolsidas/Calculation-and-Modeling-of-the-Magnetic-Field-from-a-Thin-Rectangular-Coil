#!/usr/bin/env python3
"""

"""
import numpy as np

def cot(t):
  if not np.isclose(np.tan(t), 0.0):
    return 1/np.tan(t)


  else:
      return 0

#Other option to work with vector input t


def cot(t):
  if not np.all(np.isclose(np.tan(t), 0.0)):
    return 1/np.tan(t)

  else:
      return 0
