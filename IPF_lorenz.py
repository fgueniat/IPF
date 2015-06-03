import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import dyn as d
import plot_tools
import statsmodels.distributions.empirical_distribution as stats
from scipy.stats.kde import gaussian_kde
import time
point = np.array([ -4.4090617 ,   0.94099541,  31.65011104])


ref_param = d.particle_parameters(True)
n_particle = 15
ps_param = [d.particle_parameters(False) for i in range(n_particle)]
ref_param.set_verbose(False)
ndim = ref_param.get_dim()


n_t = 100

rec_traj = np.zeros((ndim,n_t))
traj = np.zeros((ndim,n_t))

isverbose = False

lpath = [np.zeros((ndim,n_t)) for i in range(n_particle)]

density = d.density_particle(n_particle, ref_param,ps_param)
for i_p in range(0,n_particle):
	lpath[i_p][:,0] = density.get_p(i_p).get_current_position()

	rec_traj[:,0] = rec_traj[:,0] + lpath[i_p][:,0]
rec_traj[:,0] = rec_traj[:,0]/n_particle
traj[:,0] = density.get_current_position()


for i_t in range(1,n_t):
	print(i_t)
	density.next_step() # physical part

# for verifications and plot
	rec_traj[:,i_t] = density.get_estimate_position()
	traj[:,i_t] = density.get_current_position()
#	X_ips = density.get_current_positions()
	for i_p in range(0,n_particle):
		lpath[i_p][:,i_t] = density.get_current_positions()[:,i_p]


x = (traj[0,:],rec_traj[0,:])
y = (traj[1,:],rec_traj[1,:])
z = (traj[2,:],rec_traj[2,:])
plot_tools.multiplot3(x,y,z,['-r','--k'])

for i_p in range(n_particle):
	x = x + (lpath[i_p][0,:],)
	y = y + (lpath[i_p][1,:],)
	z = z + (lpath[i_p][2,:],)
plot_tools.multiplot3(x,y,z,['-r','--k'])

plot_tools.multiplot1(((traj[0,:],rec_traj[0,:])),['-r','--k'])

plot_tools.multiplot1(x,['-r','-k'])
issave = False
if issave is True:
	mydate = time.strftime("%Y_%m_%d_%Hh%M")
	s = '../' + mydate
	np.savez(s, lp=lp,rec=rec,ref_p=ref_p,ref_param=ref_param)

isload = False
if isload is True:
	mydate = '2015_06_01_19h37' + 'npz'
	s = '../'
	data = nl.load(s)
	lp = data['lp']
	rec = data['rec']
	ref_p = data['ref_p']
	ref_param = data['ref_param']
	n_particle = len(lp)
	n_t = ref_p.path[:,1].size - 1
	
