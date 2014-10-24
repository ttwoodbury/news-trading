
from __future__ import division
import random
import matplotlib.pyplot as pyplot
import numpy


SECS = 23400
times = range(SECS)
times = [x/SECS for x in times]

#------MODEL CONSTANTS ---------
sigma_mu = 1
sigma_u =1
sigma_v = 1
sigma_e = 1
Sigma_0 = 1

random.seed(10)
#SIMULATING RANDOMNESS
d_mu = [random.normalvariate(0, sigma_mu) for x in range(SECS)]

d_e = [random.normalvariate(0, sigma_e) for x in range(SECS)]

d_v = [random.normalvariate(0, sigma_v) for x in range(SECS)]

#Find the root of a cubic used in setting the CONSTANTS
A = sigma_v**2/(sigma_v**2+Sigma_0)
p = [A*sigma_e**2/sigma_v**2-sigma_e**4/sigma_v**4, \
 A*(2*sigma_e**2/sigma_v**2+1)-2*(sigma_e**2/sigma_v**2), A*sigma_e**2/sigma_v**2-(2+sigma_e**2/sigma_v**2)**2,A]
roots = numpy.roots(p)
for x in roots:
	if 0<= x <= 1:
		g = x
		break

C_BF = sigma_u/(Sigma_0+sigma_v**2)**(1/2)*1/(1+sigma_e**2/sigma_v**2*g)**(1/2)* \
(1+(1-g)*sigma_v**2/Sigma_0*(1+sigma_e**2/sigma_v**2+sigma_e**2/sigma_v**2*g)/(2+sigma_e**2/sigma_v**2+sigma_e**2/sigma_v**2*g))
B_F = [1/(1-x)*C_BF for x in times]

gamma_f = sigma_u/sigma_v*g**(1/2)
lamba_f = (Sigma_0+sigma_v**2)**(1/2)/sigma_u*(1/((1+sigma_e**2/sigma_v**2*g)**(1/2)*(1+g)))
mu_f = (1+g)/(2+sigma_e**2/sigma_v**2+sigma_e**2/sigma_v**2*g)
roe_f = sigma_v/sigma_e*g**(1/2)/(1+g)

v = []
v.append(random.normalvariate(0,Sigma_0))
for i in range(SECS-1):
	v.append(v[i]+d_v[i])

q = []
q.append(0)
for i in range(SECS-1):
	q.append(q[i]+mu_f*(d_v[i]+d_e[i]-roe_f*(d_mu[i]+B_F[i]*(v[i]-q[i])*(1/SECS)+gamma_f*d_v[i])) \
		+lamba_f*d_mu[i]+lamba_f*(B_F[i]*(v[i]-q[i])*(1/SECS)+gamma_f*d_v[i]))

x = []
x.append(0)
for i in range(SECS-1):
	x.append(x[i]+B_F[i]*(v[i]-q[i])*(1/SECS)+gamma_f*d_v[i])

p = []
p.append(0)
for i in range(SECS-1):
	p.append(q[i]+lamba_f*(B_F[i]*(v[i]-q[i])*1/SECS +gamma_f*d_v[i] +d_mu[i]))

orders = [x[i+1]-x[i] for i in range(SECS-1)]
orders.append(orders[SECS-2])

errors = [v[i]-q[i] for i in range(SECS)]

pyplot.plot(times, errors)
pyplot.show()
