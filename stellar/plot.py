#!/usr/bin/python

# import stuff
import os.path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Change matplotlib fonts
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('figure', autolayout=True)
plt.rc('font', size=18)

# Supernova prob.
sn = np.loadtxt('supernova.dat')
plt.plot(sn[:,0],sn[:,1])
plt.xlabel('Cluster size')
plt.ylabel('Prob. of supernova')
plt.savefig('supernova.pdf',transparent=True)
plt.close()

# quartiles
quarts = np.loadtxt('quartiles.dat')
print 'LOWER QUARTILE: '+str(np.percentile(quarts,25))
print 'MEDIAN: '+str(np.percentile(quarts,50))
print 'UPPER QUARTILE: '+str(np.percentile(quarts,75))
print 'MEAN: '+str(np.mean(quarts))

# big cluster
def plotit(datafile,ystr):
	data = np.loadtxt(datafile+'.log')
	mean = []
	upper = []
	lower = []
	times = range(0,16)
	for t in times:
		sub = data[:,t]
		mean.append(np.mean(sub))
		if t == 12:
			print datafile
			print np.mean(sub)
		lower.append(np.percentile(sub,25))
		upper.append(np.percentile(sub,75))
	plt.close()
	p1, = plt.plot(times,mean,label='mean')
	p2, = plt.plot(times,upper,label='Upper quart.')
	p3, = plt.plot(times,lower,label='Lower quart.')
	plt.xlabel('Time, Gyr')
	plt.ylabel(ystr)
	plt.legend(handles=[p1, p2, p3])
	plt.savefig(datafile+'.pdf',transparent=True)

plotit('ms','No. of stars in MS')
plotit('ns','No. of NS')
plotit('bh','No. of BH')
plotit('lm','Lum.-mass ratio, Gyr')
plotit('sm','Cluster Mass, $M_\\mathrm{sun}$')
plotit('wd','No. of WD')
