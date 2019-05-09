#!/usr/bin/python

# import stuff
import os.path
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Change matplotlib fonts
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('figure', autolayout=True)
plt.rc('font', size=18)

# First
string = ['water.out', 'al.out', 'iodine.out', 'lead.out']
col = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
label = ['Water', 'Aluminium', 'Iodine', 'Lead']
plt.close()
hand = []
for i in range(len(string)):
	data = np.loadtxt(string[i])
	energy = data[:,0]
	mean = data[:,1]
	med = data[:,2]
	low = data[:,3]
	high = data[:,4]
	plt.fill_between(energy,low,high,color=col[i],alpha=0.15)
	p, = plt.plot(energy,mean,color=col[i],linestyle="-",label=label[i])
	plt.plot(energy,med,color=col[i],linestyle=":")
	hand.append(p)
plt.xlabel('Energy, MeV')
plt.ylabel('Typical path length, cm')
plt.legend(handles=hand)
plt.savefig('first.pdf',transparent=True)
plt.close()

# Second, I
#-- Generate Data -----------------------------------------
# Using linspace so that the endpoint of 360 is included...
data = np.loadtxt('secondi.dat')
r = data[:,0]
n = len(data[0,:])-1
displ = 0.5 * 180 / (n+1)
theta = np.radians(np.linspace(displ, 180-displ, n))
values = data[:,range(1,len(data[0,:]))]

peaka = []
for row in range(len(r)):
	index = np.argmax(values[row,:])
	peaka.append(theta[index])

#-- Plot... ------------------------------------------------
#fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
fig = plt.figure()
ax = fig.add_subplot(111,polar=True)
ax.contourf(theta, r, values)
p, = plt.polar(peaka,r,color='black',linestyle='--',label='Peak')
plt.legend(handles=[p],loc='lower center')
plt.savefig('secondi.pdf')
plt.close()

# Second, II
#-- Generate Data -----------------------------------------
data = np.loadtxt('secondii.dat')
ine = data[:,0]
n = len(data[0,:])-1
oute = np.linspace(0,1,n)
values = data[:,range(1,len(data[0,:]))]

peake = []
for row in range(len(ine)):
	index = np.argmax(values[row,:])
	peake.append(oute[index])

#-- Plot... ------------------------------------------------
fig, ax = plt.subplots()
ax.contourf(ine, oute, values.T)
plt.xlabel('$h\\nu$, MeV')
plt.ylabel('$1-h\\nu\'/h\\nu$')
p, = plt.plot(ine,peake,linestyle='--', color='black',label='Peak')
plt.legend(handles=[p])
plt.savefig('secondii.pdf')
plt.close()

# Second II part 2
data = np.loadtxt('pairs.dat')
ine = data[:,0]
max_angle = data[:,1]
plt.plot(ine,max_angle)
plt.xlabel('$h\\nu$, MeV')
plt.ylabel('$\\theta$ for most likely energy, deg')
plt.savefig('scatter.pdf')
