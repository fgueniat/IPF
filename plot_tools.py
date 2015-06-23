# biblio plot3
import pylab
from matplotlib import pyplot
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

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

def plot2(figure,a,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):

	pylab.ion()
	if figure is False:
		fig, ax = pyplot.subplots()
	else:
		fig = figure[0]
		ax = figure[1]
#	fig = pylab.figure()
	pylab.plot(a, b,mark_col)
	return (fig,ax)

def plot1(figure,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):

	a = range(0,np.size(b))
	figure = plot2(figure,a,b,mark_col,a_label,b_label)
	return figure


def multiplot1(figure,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):
	a = (np.array(range(0,b[0].size)),)
	for i in range(len(b)-1):
		a = a + (np.array(range(0,b[0].size)),)
	fig = multiplot2(figure,a,b,mark_col,a_label,b_label)
	return fig

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


