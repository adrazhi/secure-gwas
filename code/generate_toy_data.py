import sys
import numpy as np

dir_name = sys.argv[1]
n = int(sys.argv[2])
m = int(sys.argv[3])

def save(name, vals):
	np.savetxt(name + ".txt", vals, fmt="%i " * (vals.shape[1] - 1) + "%i")

cov = np.random.choice(2, (n, 10))
pheno = np.random.choice(2, (n, 1))
geno = np.random.choice(3, (n, m))

chrom = 22
upper = 14000000
lower = 19000000
pos = np.hstack((np.full((n, 1), chrom), np.sort(upper + np.random.choice(lower - upper, (n, 1), replace=False), axis=0)))

save('../' + dir_name + '/cov', cov)
save('../' + dir_name + '/pheno', pheno)
save('../' + dir_name + '/geno', geno)
save('../' + dir_name + '/pos', pos)
