import numpy as np
import pandas as pd
import random

def distance(nUsers,dOuterRadius,dInnerRadius):
    theta = 2*np.pi * (np.random.rand(nUsers,1));
    r = (dOuterRadius-dInnerRadius) * np.sqrt(np.random.rand(nUsers, 1)) + dInnerRadius;
    x = r * np.cos(theta);
    y = r * np.sin(theta);
    return np.sqrt(x**2 + y**2);