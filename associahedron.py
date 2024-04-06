import numpy as np
import matplotlib.pyplot as plt
from codewords import *

## Step 1: Compute all codewords
n = 4
an = Associahedron(n)

## Step 2: Compute triangle area vectors for each codeword
X = np.zeros((len(an.codewords), n+2))
t = np.linspace(0, 2*np.pi, n+3)[0:n+2]
Y = np.zeros((n+2, 2))
Y[:, 0] = np.cos(t)
Y[:, 1] = np.sin(t)
for i, c in enumerate(an.codewords):
    tris = c["c"].get_tris()
    for [a, b, c] in tris:
        Yt = Y[[a, b, c], :]
        u = Yt[1, :] - Yt[0, :]
        v = Yt[2, :] - Yt[0, :]
        area = np.cross(u, v)
        X[i, [a, b, c]] += area

## Step 3: Compute principal axes
X -= np.mean(X, axis=0, keepdims=True) # Center at the origin before projection
D = (X.T).dot(X)
lam, U = np.linalg.eigh(D)
XU = X.dot(U)[:, -3:] # Take out the 3 dimensions with the highest variance
plt.figure(figsize=(10, 10))
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.scatter(XU[:, 0], XU[:, 1], XU[:, 2], c='k')


for i in range(X.shape[0]):
    w1 = an.codewords[i]["c"].w
    for j in range(X.shape[0]):
        w2 = an.codewords[j]["c"].w
        if are_rotation_neighbors(w1, w2):
            Xij = XU[[i, j], :]
            plt.plot(Xij[:, 0], Xij[:, 1], Xij[:, 2], c='k', linewidth=1)

for i in range(XU.shape[0]-1):
    c = [1, i/XU.shape[0], 1]
    plt.plot(XU[[i, i+1], 0], XU[[i, i+1], 1], XU[[i, i+1], 2], c=c, linewidth=2)
    
    
for i, c in enumerate(an.codewords):
    w = c["c"].w
    #ax.text(XU[i, 0], XU[i, 1], XU[i, 2], "{}".format(w))
    
plt.figure()
plt.stem(lam)
plt.show()
