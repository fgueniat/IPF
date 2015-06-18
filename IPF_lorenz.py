import numpy as np
import matplotlib.pyplot as plt
import dyn as d
import plot_tools

import time
point = np.array([ -4.4090617 ,   0.94099541,  31.65011104])


ref_param = d.particle_parameters(True)
n_particle = 10
ps_param = [d.particle_parameters(False) for i in range(n_particle)]
ref_param.set_verbose(False)
ndim = ref_param.get_dim()


n_t = 1000

rec_traj = np.zeros((ndim,n_t))
traj = np.zeros((ndim,n_t))
obs = np.zeros((n_t))

isverbose = False

lpath = [np.zeros((ndim,n_t)) for i in range(n_particle)]

density = d.density_particle(n_particle, ref_param,ps_param)
for i_p in range(0,n_particle):
	lpath[i_p][:,0] = density.get_p(i_p).get_current_position()

	rec_traj[:,0] = rec_traj[:,0] + lpath[i_p][:,0]
rec_traj[:,0] = rec_traj[:,0]/n_particle
traj[:,0] = density.get_current_position()

obs_strat = np.zeros(n_t)
for i_t in range(1,n_t):
	if np.mod(i_t,5)==0:
		obs_strat[i_t] = 1
compteur_obs = 0
for i_t in range(1,n_t):
	print(i_t)
	
	if obs_strat[i_t]==1:
		density.new_data_provided(True)
	else:
		density.new_data_provided(False)
		compteur_obs = compteur_obs+1

	density.next_step() # physical part

# for verifications and plot
	traj[:,i_t] = density.get_current_position()
	obs[i_t] = density.p_ref.get_obs()
	if obs_strat[i_t]==1:
		rec_traj[:,i_t] = density.get_estimate_position()

		intermediate_steps = density.get_estimate_intermediate_position()

		for i_obs in range(0,compteur_obs):

			rec_traj[:,i_t - (compteur_obs - i_obs)] = intermediate_steps[i_obs]
			
#	X_ips = density.get_current_positions()
		for i_p in range(0,n_particle):
			lpath[i_p][:,i_t] = density.get_current_positions()[:,i_p]

			intermediate_steps = density.get_estimate_intermediate_positions(i_p)
			for i_obs in range(0,compteur_obs):
				lpath[i_p][:,i_t - (compteur_obs - i_obs)] = intermediate_steps[i_obs]
#				rec_traj[:,i_t - (compteur_obs - i_obs)] = rec_traj[:,i_t - (compteur_obs - i_obs)] + density.get_p(i_p).get_weight() * lpath[i_p][:,i_t - (compteur_obs - i_obs)]

				

# Initialize compteur_obs
	if obs_strat[i_t]==1:
		compteur_obs = 0


x = (traj[0,:-10],rec_traj[0,:-10])
y = (traj[1,:-10],rec_traj[1,:-10])
z = (traj[2,:-10],rec_traj[2,:-10])
plot_tools.multiplot3(x,y,z,['-r','--k'])

for i_p in range(n_particle):
	x = x + (lpath[i_p][0,:-10],)
	y = y + (lpath[i_p][1,:-10],)
	z = z + (lpath[i_p][2,:-10],)

plot_tools.multiplot3(x,y,z,['-r','--k'])


plot_tools.multiplot1(((traj[0,:-10],rec_traj[0,:-10])),['-r','--k'])

plot_tools.multiplot1(x,['-r','-k'])

ind = np.where(obs_strat[:-10]==1)
y = (traj[0,ind[0]],traj[0,:-10])
x = (ind[0],np.array(range(0,y[1].size)))

plot_tools.multiplot2(x,y,['or','.k'])

y = (obs[ind[0]],obs[:-10])
x = (np.array(ind[0]),range(0,y[1].size))

plot_tools.multiplot2(x,y,['or','.k'])


issave = False
if issave is True:
	mydate = time.strftime("%Y_%m_%d_%Hh%M")
	s = '../' + mydate
	np.savez(s, obs = obs, traj=traj,rec_traj=rec_traj,density=density,lpath=lpath,ref_param=ref_param)

isload = False
if isload is True:
	mydate = '2015_06_07_12h53' + '.npz'
	s = '../' + mydate
	data = np.load(s)
	traj = data['traj']
	rec_traj = data['rec_traj']
	rec_traj = data['rec_traj']
	density = data['density']
	lpath = data['lpath']
	ref_param = data['ref_param']
	
