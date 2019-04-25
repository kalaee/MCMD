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
def avg(func,energy,lng,temp,zpart):
	s = 0
	for i in range(len(energy)):
		s += func(energy[i])*np.exp(lng[i]-energy[i]/temp)
	return s / zpart

def unit(x):
	return 1

def first(x):
	return x

def second(x):
	return x**2

def cvt(temp,energy,lng):
	c = []
	for t in temp:
		zpart = avg(unit,energy,lng,t,1)
		x1t = avg(first,energy,lng,t,zpart)
		x2t = avg(second,energy,lng,t,zpart)
		c.append((x2t-x1t**2)/t**2)
	return c

def q1(proj):
	name = 'proj'+proj+'/lng_'+str(24)
	data = np.loadtxt(name)
	energy = data[:,0]
	lng = data[:,1]
	temp = np.linspace(0.5,1,1000)
	c = cvt(temp,energy,lng)
	# find index of largest element
	cmax = max(c)
	tmax = temp[c.index(cmax)]
	plt.close()
	plt.plot([tmax, tmax],[0,cmax],linestyle=':',color="tab:orange")
	plt.plot(temp,c,color='tab:blue')
	plt.xlabel('Temperature')
	plt.ylabel('$C_V(T)$')
	s = '{:0.2f}'.format(tmax)
	print 'Question 1, proj'+proj+': T_max = '+s
	plt.savefig('proj'+proj+'_1.pdf')
	return

q1('A')
q1('B')

# Question 2, PDF for Tmax
def pdf(energy,lng,tmax):

