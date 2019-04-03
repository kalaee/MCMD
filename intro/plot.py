#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize as optimization

# Change matplotlib fonts
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('figure', autolayout=True)
plt.rc('font', size=18)

# Plotting for question 1
## data
q1 = np.loadtxt('q1.out')
zoom = []
for q in q1:
	if q < 50:
		zoom.append(q)
dense = len(zoom) / (1.0 * len(q1)) # weight of zoom

## prob. function
t = np.linspace(0.0001,50,1000)
f = np.sin(t)*np.sin(t)/t/t / 1.571
for i in range(len(f)):
	if f[i] < 1e-6:
		f[i] = 1e-6
	f[i] /= dense # compensate for zoom

plt.hist(zoom,200,log=True,normed=True)
plt.plot(t,f,color='tab:orange', label='Function')
plt.xlabel('$x$')
plt.ylabel('Probability, $p(x)\\propto \\sin^2 x / x^2$')
plt.savefig('q1.pdf')

# Plotting for question 2

