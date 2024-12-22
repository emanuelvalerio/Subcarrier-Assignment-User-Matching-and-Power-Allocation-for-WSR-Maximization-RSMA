import numpy as np

def calculationFc(hn1,hn2,Pn1,Pn2,rho):
    y1 = np.sqrt(1+(np.linalg.norm(hn1)**2)*rho*Pn1);
    y2 = np.sqrt(1+(np.linalg.norm(hn2)**2)*rho*Pn2);
    H1 = hn1/y1; # vector of antennas with dimension 4x1
    H2 = hn2/y2; # vector of antennas with dimension 4x1
    H1 = H1.reshape((-1,1));
    H2 = H2.reshape((-1,1));
    alphaMatrix = np.dot(np.vstack((np.conj(H1).T, np.conj(H2).T)),np.hstack((H1,H2)));
    a = (alphaMatrix[1,1] - abs(alphaMatrix[0,1])) ;
    b = (alphaMatrix[0,0] - abs(alphaMatrix[0,1]));
    muMatrix = (1 / (alphaMatrix[0, 0] + alphaMatrix[1, 1] - 2 * np.abs(alphaMatrix[0, 1]))) * np.array([a, b]);
    lambda1 = (alphaMatrix[0,0]*alphaMatrix[1,1] - np.abs(alphaMatrix[0,1])**2) / ((alphaMatrix[0,0] + alphaMatrix[1,1] - 2*np.abs(alphaMatrix[0,1])));
    fc = (1 / np.sqrt(lambda1)) * (muMatrix[0] * H1 + muMatrix[1] * H2 * np.exp(-1j * np.angle(alphaMatrix[0, 1])));
    return fc;
