# biblio plot3
import pylab
from matplotlib import pyplot
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def plot3c(a,b,c,z,a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):
	p = np.argsort(z)
	macol = pyplot.cm.seismic(np.arange(z.size))
	ylab.ion()
	fig = pylab.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(a[p], b[p], c[p], c=macol,edgecolor='')
	ax.set_xlabel(a_label)
	ax.set_ylabel(b_label)
	ax.set_zlabel(c_label)

	return fig


def plot3(a,b,c,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):

	pylab.ion()
	fig = pylab.figure()
	ax = Axes3D(fig)
	ax.plot(a, b, c,mark_col)
	ax.set_xlabel(a_label)
	ax.set_ylabel(b_label)
	ax.set_zlabel(c_label)
	return fig

def plot2(a,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):

	pylab.ion()
	fig = pylab.figure()
	pylab.plot(a, b,mark_col)
	return fig

def plot1(b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):

	a = range(0,np.size(b))
	fig = plot2(a,b,mark_col,a_label,b_label)
	return fig

def multiplot2(a,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):


	pylab.ion()
	fig, ax = pyplot.subplots()
	for i in range(0,len(b)):
		j = len(b)-i-1
		try:
			mci = mark_col[j]
		except:
			mci = '-k'
		ax.plot(a,b[j],mci)
	ax.set_xlabel(a_label)
	ax.set_ylabel(b_label)

def multiplot1(b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):
	a = np.array(range(0,b[0].size))
	fig = multiplot2(a,b,mark_col,a_label,b_label)
	return fig

def multiplot3(a,b,c,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):

	pylab.ion()
	fig = pylab.figure()
	ax = Axes3D(fig)
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
	return fig
