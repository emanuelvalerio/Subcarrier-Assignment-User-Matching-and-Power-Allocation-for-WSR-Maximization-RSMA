import numpy as np
hn1 = np.array([1, 2])
hn2 = np.array([4, 5])
y1 = 2
y2 = 3

H1 = hn1 / y1
H2 = hn2 / y2

alphaMatrix = np.dot(np.array([H1, H2]).T, np.array([H1, H2]))
print(alphaMatrix)
