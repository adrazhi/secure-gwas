import sys
import numpy as np

f1 = sys.argv[1]
f2 = sys.argv[2]
x1 = np.loadtxt(f1)
x2 = np.loadtxt(f2)
print("Correlation: ", np.corrcoef(x1, x2))