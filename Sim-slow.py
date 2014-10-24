
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

#SIMULATING RANDOMNESS
d_mu = [random.normalvariate(0, sigma_mu) for x in range(SECS)]

d_e = [random.normalvariate(0, sigma_e) for x in range(SECS)]

d_v = [random.normalvariate(0, sigma_v) for x in range(SECS)]

#------ SLOW TRADING ---------
C_BS = sigma_mu/(Sigma_0)**(1/2)*(1+sigma_v**2*sigma_e**2/(Sigma_0*(sigma_v**2+sigma_e**2)))**(1/2)
B_S = [1/(1-x)*C_BS for x in times]

lamba_s = Sigma_0**(1/2)/sigma_u*(1+sigma_v**2*sigma_e**2/(Sigma_0*(sigma_v**2+sigma_e**2)))**(1/2)

mu_s = sigma_v**2/(sigma_v**2+sigma_e**2)


v = []
v.append(random.normalvariate(0,Sigma_0))
for i in range(SECS-1):
	v.append(v[i]+d_v[i])

q = []
q.append(0)
for i in range(SECS-1):
	q.append(q[i]+mu_s*(d_v[i]+d_e[i])+lamba_s*d_mu[i]+lamba_s*B_S[i]*(v[i]-q[i])*(1/SECS))

x = []
x.append(0)
for i in range(SECS-1):
	x.append(x[i]+B_S[i]*(v[i]-q[i])*1/SECS)

p = []
p.append(0)
for i in range(SECS-1):
	p.append(q[i]+mu_s*(d_v[i]+d_e[i])+lamba_s*(B_S[i]*(v[i]-q[i])*1/SECS +d_mu[i]))

pyplot.plot(times, v)
pyplot.plot(times, x)
pyplot.plot(times, p)
pyplot.show()

