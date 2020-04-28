import sys
import numpy as np

pref1 = sys.argv[1]
pref2 = sys.argv[2]
x1 = np.loadtxt(pref1 + "_assoc.txt")
x2 = np.loadtxt(pref2 + "_assoc.txt")

def get_gkeep(pref):
	gkeep = np.loadtxt(pref + "_gkeep1.txt")
	gkeep2 = np.loadtxt(pref + "_gkeep2.txt")
	ind = 0;
    for j in range(gkeep.shape[0]):
    	if gkeep[j] == 1:
    		gkeep[j] = gkeep2[ind] == 1
    		ind += 1
    return gkeep

if (x1.shape != x2.shape):
	gkeep1 = get_gkeep(pref1)
	gkeep2 = get_gkeep(pref2)
	ind1 = 0
	ind2 = 0
	for i in range(gkeep1.shape[0]):
		if gkeep1[i] == 1 and gkeep2[i] == 1:
			ind1 += 1
			ind2 += 1
		elif gkeep1[i] == 1:
			x1 = np.delete(x1, ind1)
		elif gkeep2[i] == 1:
			x2 = np.delete(x2, ind2)

print("len of x1: ", x1.shape[0])
print("len of x2: ", x2.shape[0])

print("Top 10 Positions for X1: ", x1.argsort()[-10:][::-1])
print("Top 10 Positions for X2: ", x2.argsort()[-10:][::-1])
print("Average Absolute Difference: ", np.average(np.abs(x1 - x2)))
print("Correlation: ", np.corrcoef(x1, x2)[0, 1])
print("Distance: ", np.linalg.norm(x1 - x2))