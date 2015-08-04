# biblio plot3
import pylab
from matplotlib import pyplot
import pylab as pylab
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mtick

def plot3c(figure,a,b,c,z,a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):
	p = np.argsort(z)
	macol = pyplot.cm.seismic(np.arange(z.size))
	pylab.ion()
	if figure is False:
		fig = pylab.figure()
		ax = fig.add_subplot(111, projection='3d')
	else:
		fig = figure[0]
		ax = figure[1]
	ax.scatter(a[p], b[p], c[p], c=macol,edgecolor='')
	ax.set_xlabel(a_label)
	ax.set_ylabel(b_label)
	ax.set_zlabel(c_label)

	return (fig,ax)


def plot3(figure,a,b,c,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):

	pylab.ion()

	if figure is False:
		fig = pylab.figure()
		ax = Axes3D(fig)
	else:
		fig = figure[0]
		ax = figure[1]
	ax.plot(a, b, c,mark_col)
	ax.set_xlabel(a_label)
	ax.set_ylabel(b_label)
	ax.set_zlabel(c_label)
	return (fig,ax)

#def plot2(figure,a,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):

#	pylab.ion()
#	if figure is False:
#		fig, ax = pyplot.subplots()
#	else:
#		fig = figure[0]
#		ax = figure[1]
##	fig = pylab.figure()
#	pylab.plot(a, b,mark_col)
#	return (fig,ax)

#def plot1(figure,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):

#	a = range(0,np.size(b))
#	if figure is False:
#		fig, ax = pyplot.subplots()
#	else:
#		fig = figure[0]
#		ax = figure[1]
#	figure = (fig,ax)
#	plot2(figure,a,b,mark_col,a_label,b_label)
#	return  (fig,ax)


def multiplot1(figure,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):
	a = (np.array(range(0,b[0].size)),)
	for i in range(1,len(b)):
		a = a + (np.array(range(0,b[i].size)),)
	if figure is False:
		fig, ax = pyplot.subplots()
	else:
		fig = figure[0]
		ax = figure[1]
	figure = (fig,ax)
	figure = multiplot2(figure,a,b,mark_col,a_label,b_label)
	return figure

def multiplot2(figure,a,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):


	pylab.ion()
	if figure is False:
		fig, ax = pyplot.subplots()
	else:
		fig = figure[0]
		ax = figure[1]

	for i in range(0,len(b)):
		j = len(b)-i-1
		try:
			mci = mark_col[j]
		except:
			mci = '-k'
		ax.plot(a[j],b[j],mci)
	ax.set_xlabel(a_label)
	ax.set_ylabel(b_label)

	return (fig,ax)


def multiplot3(figure,a,b,c,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):

	pylab.ion()
	if figure is False:
		fig = pylab.figure()
		ax = Axes3D(fig)
	else:
		fig = figure[0]
		ax = figure[1]
	for i in range(0,len(a)):
		j = len(a)-i-1
		try:
			mci = mark_col[j]
		except:
			mci = '-k'
		ax.plot(a[j], b[j], c[j],mci)
	ax.set_xlabel(a_label)
	ax.set_ylabel(b_label)
	ax.set_zlabel(c_label)
	return (fig,ax)

def closeall(n=50):
	for i in range(n):
		pyplot.close()
		

def flim1(figure,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):

	if figure is False:
		fig = pyplot.figure( 1 )
		ax = fig.add_subplot( 111 )
		im = ax.imshow( np.zeros( ( 256, 256, 3 ) ) ) # Blank starting image
		figure = (fig,ax,im)
		fig.show()
		im.axes.figure.canvas.draw()
	else:
		fig = figure[0]
		ax = figure[1]
		im = figure[2]
	im.set_data( b )
	im.axes.figure.canvas.draw()
	return (fig,ax,im)

def save(path,i=-1):
	if i == -1:
		s = path
	else:
		s = path + 'temp_0000'
		s2 = str(i)
		s = s[:-len(s2)] + s2
		s = s + '.png'
	pyplot.savefig(s)


class Paradraw():
	def __init__(self,marks="-k",x_label = 'axis 1',y_label = 'axis 2',z_label = 'axis 3',c_label = 'colors', colmap = 'hot_r',xlim = False,ylim = False,zlim = False, x_scale = 'linear',y_scale = 'linear',z_scale = 'linear',title='',iscolorbar = True,fontsize = 20,ticksize = 16, x_tick_label = False,y_tick_label = False,cbar_tick_label = False,axformat = '%.1e',cbformat = '%.1e'):

#plot style
		self.marks = marks
		self.colmap = colmap

#axis labels
		self.x_label = x_label
		self.x_tick_label = x_tick_label
		self.y_tick_label = y_tick_label
		self.cbar_tick_label = y_tick_label
		self.y_label = y_label
		self.z_label = z_label
		self.c_label = c_label
		self.title = title

#axis style
		self.x_scale = x_scale
		self.y_scale = y_scale
		self.z_scale = z_scale

#axis limits
		self.xlim = xlim
		self.ylim = ylim
		self.zlim = zlim

		self.iscolorbar = iscolorbar

#fonts
		self.fontsize = fontsize
		self.ticksize = ticksize
		self.cbformat = cbformat
		self.axformat = axformat
		
def pcolor(data,param):

	x = data[0]
	y = data[1]
	z = data[2]
	if param.zlim is False:
		zmin = z.min()
		zmax = z.max()
	else:
		zmin = param.zlim[0]
		zmax = param.zlim[1]

	pylab.ion()
	fig = pylab.figure()
	ax = fig.add_subplot(111) 
	cax = ax.pcolor(x,y,z, cmap=param.colmap, vmin=zmin, vmax=zmax)
	# set the limits of the plot to the limits of the data

	cbar = fig.colorbar(cax,format=param.cbformat)
	if param.iscolorbar is False:
		fig.delaxes(fig.axes[1]) #cbar IS axes[1] while ax is axes[0]
		Options(ax,(x,y),param)
	else:
		Options(ax,(x,y,z),param,cbar)
	
	pylab.show()
	return fig,ax,cbar
	

def plot1(data,param):
	x = np.arange(len(data))
	y = data
	data2 = (x,y)
	fig,ax = plot2(data2,param)
	return fig,ax

def plot2(data,param):
	x = data[0]
	y = data[1]

	pylab.ion()
	fig = pylab.figure()
	ax = fig.add_subplot(111) 
	ax.plot( x, y, param.marks)
	# set the limits of the plot to the limits of the data
	Options(ax,(x,y),param)
	pylab.show()
	return fig,ax

def Options(ax,X,param,cbar = False):
	x=X[0]
	y=X[1]
	
	if param.xlim is False:
		ax.set_xlim(  [ x.min(),x.max() ]  )
	else:
		ax.set_xlim(  [ param.xlim[0],param.xlim[1] ]  )
	if param.ylim is False:
		ax.set_ylim(  [ y.min(),y.max() ]  )
	else:
		ax.set_ylim(  [ param.ylim[0],param.ylim[1] ]  )
	ax.set_xscale(param.x_scale)
	ax.set_yscale(param.y_scale)
	ax.set_xlabel(param.x_label,None,None,fontsize=param.fontsize)
	ax.set_ylabel(param.y_label,None,None,fontsize=param.fontsize)
	if param.x_tick_label is False:
		if param.x_scale == 'log':
			xlabel = np.exp(np.linspace(np.log(x.min()),np.log(x.max()),3))
		else:
			xlabel = np.linspace(x.min(),x.max(),3)
	else:
		xlabel = param.x_tick_label
	if param.y_tick_label is False:
		if param.y_scale == 'log':
			ylabel = np.exp(np.linspace(np.log(y.min()),np.log(y.max()),3))
		else:
			ylabel = np.linspace(y.min(),y.max(),3)
	else:
		ylabel = param.y_tick_label
	ax.set_xticks(xlabel)
	ax.set_yticks(ylabel)
	ax.set_xticklabels(ax.get_xticks(),None,None,fontsize=param.ticksize)
	ax.set_yticklabels(ax.get_yticks(),None,None,fontsize=param.ticksize)

	ax.xaxis.set_major_formatter(mtick.FormatStrFormatter(param.axformat))
	ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

	ax.set_title(param.title)
	if cbar is not False:
		z=X[2]
		cbar.set_label(param.c_label,fontsize=param.fontsize)
		if param.cbar_tick_label is False:
			if param.zlim is False:
				clabel = np.linspace(z.min(),z.max(),3)
			else:
				clabel = np.linspace(param.zlim[0],param.zlim[1],3)
				
		else:
			clabel = param.cbar_tick_label
		print(clabel)


#		cbar.set_ticklabels(clabel)
		for t in cbar.ax.get_yticklabels():
			t.set_fontsize(param.ticksize)
		cbar.set_ticks(clabel)
#		cbar.set_ticklabels(clabel)
#		cbar.ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

	pyplot.tight_layout()
#	ax.tick_params(axis='both', which='minor', labelsize=param.ticksize)


def change_options():
	ax = pyplot.gca()
#	ax.tick_params(axis='both', which='minor', labelsize=param.ticksize)
	ax.xaxis.set_minor_locator(pyplot.FixedLocator([2,4]))

