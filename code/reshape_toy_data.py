import sys
import numpy as np
import os

dir_name = sys.argv[1]
n = int(sys.argv[2])
m = int(sys.argv[3])

def save(name, vals):
	np.savetxt(name + ".txt", vals, fmt="%i " * (vals.shape[1] - 1) + "%i")

cov = np.loadtxt('../' + dir_name + '/cov.txt')
pheno = np.reshape(np.loadtxt('../' + dir_name + '/pheno.txt'), (-1, 1))
geno = np.loadtxt('../' + dir_name + '/geno.txt')
pos = np.loadtxt('../' + dir_name + '/pos.txt')

new_cov = cov[0:n,:]
new_pheno = pheno[0:n]
new_geno = geno[0:n,0:m]
new_pos = pos[0:m]

new_dir = '../' + dir_name + '2'
if not os.path.exists(new_dir):
	os.mkdir(new_dir)

save(new_dir + '/cov', new_cov)
save(new_dir + '/pheno', new_pheno)
save(new_dir + '/geno', new_geno)
save(new_dir + '/pos', new_pos)
