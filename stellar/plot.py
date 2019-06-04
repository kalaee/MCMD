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

## Supernova prob.
#sn = np.loadtxt('supernova.dat')
#plt.plot(sn[:,0],sn[:,1])
#plt.xlabel('Cluster size')
#plt.ylabel('Prob. of supernova')
#plt.savefig('supernova.pdf',transparent=True)
#plt.close()
#
## quartiles
#quarts = np.loadtxt('quartiles.dat')
#print 'LOWER QUARTILE: '+str(np.percentile(quarts,25))
#print 'MEDIAN: '+str(np.percentile(quarts,50))
#print 'UPPER QUARTILE: '+str(np.percentile(quarts,75))
#print 'MEAN: '+str(np.mean(quarts))

# big cluster
def plotit(datafile,ystr):
	print datafile
	data = np.loadtxt(datafile+'.log')
	mean = []
	upper = []
	lower = []
	if datafile == 'sm':
		rem = []
	times = np.loadtxt('time.log')
	for t in range(len(times)):
		if datafile == 'sm':
			sub = data[:,t]
			remn = data[:,0] - data[:,t]
		else:
			sub = data[:,t]
		if datafile == 'sm':
			rem.append(np.mean(remn)/1e5)
			mean.append(np.mean(sub)/1e5)
		elif datafile == 'lm':
			mean.append(1e5*np.mean(sub))
		else:
			mean.append(np.mean(sub))
		lower.append(np.percentile(sub,25))
		upper.append(np.percentile(sub,75))
	plt.close()
	p1, = plt.plot(times,mean,label='Cluster')
	if datafile == 'sm':
		p2, = plt.plot(times,rem,linestyle=':',label='Loss')
#		p2, = plt.plot(times,upper,color='tab:orange',linestyle=':',label='Quartiles')
#		p2, = plt.plot(times,lower,color='tab:orange',linestyle=':',label='Quartiles')
	plt.xscale('log')
	plt.xlabel('Time, years')
	plt.ylabel(ystr)
	if datafile == 'sm':
		plt.legend(handles=[p1, p2])
	if datafile == 'lm':
		plt.yscale('log')
	plt.savefig(datafile+'.pdf',transparent=True)
	return times, mean, upper, lower

mst, msm, msu, msl = plotit('ms','No. of stars in MS')
nst, nsm, nsu, nsl = plotit('ns','No. of NS')
bht, bhm, bhu, bhl = plotit('bh','No. of BH')
wdt, wdm, wdu, wdl = plotit('wd','No. of WD')
plotit('lm','mass/luminosity, arb.u.')
plotit('sm','Mass, $10^{5}M_\\mathrm{sun}$')

plt.close()

pms, = plt.plot(mst,msm,color='tab:blue',linestyle='-',label='MS')
pns, = plt.plot(nst,nsm,color='tab:orange',linestyle='--',label='NS')
pbh, = plt.plot(bht,bhm,color='tab:red',linestyle=':',label='BH')
pwd, = plt.plot(wdt,wdm,color='tab:green',linestyle='-.',label='WD')

plt.xscale('log')
plt.yscale('log')
plt.legend(handles=[pms,pns,pbh,pwd])
plt.xlabel('Time, yr.')
plt.ylabel('Population')
plt.savefig('remn.pdf',transparent=True)

# Distribution
t1 = np.linspace(0.08,0.5,1000)
t2 = np.linspace(0.5,1,1000)
t3 = np.linspace(1,120,1000)
f1 = np.power(t1,-1.3)*np.power(0.5,-2.2)/np.power(0.5,-1.3)
f2 = np.power(t2,-2.2)
f3 = np.power(t3,-2.7)
plt.close()
p1, = plt.plot(t1,f1,color='tab:blue',linestyle='-',label='$f(m)$')
plt.plot(t2,f2,color='tab:blue',linestyle='-',label='$f(m)$')
plt.plot(t3,f3,color='tab:blue',linestyle='-',label='$f(m)$')
t = np.linspace(0.08,120,3000)
g = np.power(t,-2.2)
p2, = plt.plot(t,g,color='tab:orange',linestyle=':',label='$g(m)$')
plt.xlabel('Stellar mass, $M_\\mathrm{sun}$')
plt.ylabel('Prob. density, arb.u.')
plt.xscale('log')
plt.yscale('log')
plt.legend(handles=[p1, p2])
plt.savefig('dist.pdf',transparent=True)

files = ['ms', 'wd', 'ns', 'bh']
for f in files:
	data = np.loadtxt(f+'12,log')
	sub = data[:,1]
	print f
	print np.mean(sub)
