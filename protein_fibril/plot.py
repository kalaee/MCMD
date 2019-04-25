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

# Question 0, convergence
def q0(proj):
	plt.close()
	bot = 20
	handles = []
	for i in range(bot,25):
		name = 'proj'+proj+'/lng_'+str(i)
		name_old = 'proj'+proj+'/lng_'+str(i-1)
		data = np.loadtxt(name)
		data_old = np.loadtxt(name_old)
		energy = data[:,0]
		dlng = data[:,1] - data_old[:,1]
		p, = plt.plot(energy,dlng,color=plt.cm.summer(1-(i-bot)/(25.-bot)),label='$f = $'+str(i))
		if (i == bot or i == 24):
			handles.append(p)
	plt.xlabel('Energy')
	plt.ylabel('ln($g(E)$)$_{f=i}$ - ln($g(E)$)$_{f=i-1}$')
	plt.legend(handles=handles)
	plt.savefig('proj'+proj+'_0.pdf')
	return

#q0('A')
#q0('B')

# Question 1, heat capacity
def pdf(energy,lng,t,zpart):
	return np.exp(lng-energy/t)/zpart

def avg(func,energy,lng,t,zpart):
	s = 0
	p = pdf(energy,lng,t,zpart)
	for i in range(len(energy)):
		s += func(energy[i])*p[i]
	return s

def unit(x):
	return 1

def first(x):
	return x

def second(x):
	return x**2

def cvt(temp,energy,lng):
	c = []
	tnew = []
	for t in temp:
		zpart = avg(unit,energy,lng,t,1)
		p = pdf(energy,lng,t,zpart)
		if p[-1] > 1e-4 or p[-2] > 1e-4:
			continue
		tnew.append(t)
		x1t = avg(first,energy,lng,t,zpart)
		x2t = avg(second,energy,lng,t,zpart)
		c.append((x2t-x1t**2)/t**2)
	return tnew, c

def q1(proj,bot,top):
	name = 'proj'+proj+'/lng_'+str(24)
	data = np.loadtxt(name)
	energy = data[:,0]
	lng = data[:,1]
	tpre = np.linspace(bot,top,1000)
	temp, c = cvt(tpre,energy,lng)
	# find index of largest element
	cmax = max(c)
	tmax = temp[c.index(cmax)]
	plt.close()
	plt.plot([tmax, tmax],[0,cmax],linestyle=':',color="tab:orange")
	plt.plot(temp,c,color='tab:blue')
	plt.xlabel('Temperature')
	plt.ylabel('$C_V(T)$')
	s = '{:0.3f}'.format(tmax)
	print 'Question 1, proj'+proj+': T_max = '+s
	plt.savefig('proj'+proj+'_1.pdf')
	return tmax

maxA = q1('A',0.5,0.8)
maxB = q1('B',0.5,0.7)

# Question 2, PDF for Tmax
def q2(proj,tmax,width):
	name = 'proj'+proj+'/lng_'+str(24)
	data = np.loadtxt(name)
	energy = data[:,0]
	lng = data[:,1]
	zpart = avg(unit,energy,lng,tmax,1)
	p = pdf(energy,lng,tmax,zpart)
	enew = []
	pnew = []
	for i in range(len(energy)/width):
		enew.append((energy[i*width]+energy[(i+1)*width-1])/2)
		ps = 0
		for j in range(width):
			ps += p[i*width+j]
		pnew.append(ps)
	test = 0
	for pi in pnew:
		test += pi
	plt.close()
	plt.plot(enew,pnew)
	plt.xlabel('Energy')
	plt.ylabel('$P_{T_\mathrm{max}}(E)$')
	plt.savefig('proj'+proj+'_2.pdf')
	return

q2('A',maxA,4)
q2('B',maxB,4)
