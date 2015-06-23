import numpy as np
import matplotlib.pyplot as plt
import dyn2 as d
import plot_tools
from plot_tools import closeall
import time


scenario = 'burgers'

if scenario == 'lorenz':
	ref_param = d.particle_parameters(True)
	n_particle = 10
	ps_param = [d.particle_parameters(False) for i in range(n_particle)]
	ref_param.set_verbose(False)
if scenario == 'roessler':
	ref_param = d.particle_parameters(True)
	n_particle = 10
	ps_param = [d.particle_parameters(False) for i in range(n_particle)]
	ref_param.set_verbose(False)
if scenario == 'burgers':
	x = np.linspace(-1,1,30)
#	X_init = np.exp(- x**2)
##	X_init = np.exp(- x**2) + 5*np.sinc(np.pi*x)
#	X_init = X_init - X_init.min()
#	X_init = X_init/X_init.max()
#	X_init = (2*X_init-1)/2
	X_init = np.sinc(np.pi * x)
	X_init = X_init/np.max(np.abs(X_init))

	isverbose = False
	dt = 0.01
	t0 = 0
	s_obs = np.sqrt(0.01)
	g_int = np.sqrt(0.01)
	objective = 'filter'

	ref_param = d.particle_parameters(True,X_init,isverbose, d.f_burgers,d.h_burgers,dt,t0,s_obs,g_int,objective)

	n_particle = 10
	ps_param = [d.particle_parameters(False,X_init,isverbose, d.f_burgers,d.h_burgers,dt,t0,s_obs,g_int,objective) for i in range(n_particle)]
	ref_param.set_verbose(False)


ndim = ref_param.get_dim()


n_t = 15

rec_traj = np.zeros((ndim,n_t))
traj = np.zeros((ndim,n_t))
obs = np.zeros((ref_param.get_dim_obs(),n_t))

isverbose = False

lpath = [np.zeros((ndim,n_t)) for i in range(n_particle)]

density = d.density_particle(n_particle, ref_param, ps_param)
for i_p in range(0,n_particle):
	lpath[i_p][:,0] = density.get_p(i_p).get_current_position()

	rec_traj[:,0] = rec_traj[:,0] + lpath[i_p][:,0]
rec_traj[:,0] = rec_traj[:,0]/n_particle
traj[:,0] = density.get_current_position()

obs_strat = np.zeros(n_t)
for i_t in range(1,n_t):
	if np.mod(i_t,1)==0:
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
	obs[:,i_t] = density.p_ref.get_obs()
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



isplot = True


if isplot is True:
	if scenario == 'lorenz':
		x = (traj[0,:-10],rec_traj[0,:-10])
		y = (traj[1,:-10],rec_traj[1,:-10])
		z = (traj[2,:-10],rec_traj[2,:-10])
		plot_tools.multiplot3(False,x,y,z,['-r','--k'])
		for i_p in range(n_particle):
			x = x + (lpath[i_p][0,:-10],)
			y = y + (lpath[i_p][1,:-10],)
			z = z + (lpath[i_p][2,:-10],)
		plot_tools.multiplot3(False,x,y,z,['-r','--k'])
		plot_tools.multiplot1(False,((traj[0,:-10],rec_traj[0,:-10])),['-r','--k'])
		plot_tools.multiplot1(False,x,['-r','-k'])
		ind = np.where(obs_strat[:-10]==1)
		y = (traj[0,ind[0]],traj[0,:-10])
		x = (ind[0],np.array(range(0,y[1].size)))
		plot_tools.multiplot2(False,x,y,['or','.k'])
		y = (obs[0,ind[0]],obs[0,:-10])
		x = (np.array(ind[0]),range(0,y[1].size))
		plot_tools.multiplot2(False,x,y,['or','.k'])
	if scenario == 'burgers':
		for it in range(i_t+1):
			x = (traj[:,it],rec_traj[:,it],obs[:,it])
			plot_tools.multiplot1(False,x,['b','k','r'])
			plot_tools.save('/media/DATA_MIXTE/code/python/IPF/fig/',it)
			closeall()


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
	
