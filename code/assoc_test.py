import sys
import numpy as np

f1 = sys.argv[1]
f2 = sys.argv[2]
x1 = np.loadtxt(f1)
x2 = np.loadtxt(f2)
print("Average Absolute Difference: ", np.average(np.abs(x1 - x2)))
print("Correlation: ", np.corrcoef(x1, x2)[0, 1])
print("Distance: ", np.linalg.norm(x1 - x2))