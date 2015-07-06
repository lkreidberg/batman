import matplotlib.pyplot as plt
import numpy as np
import math

plt.figure(figsize = (4,4))
theta = np.linspace(0, 2*math.pi, 1000)

r = 0.15
fac = 7.0e-2 
rs = []

while r < 1.0: 
	rs.append(r)
	x = r*np.cos(theta)
	y = r*np.sin(theta)

	dr = fac*math.acos(r)
	r = r + dr
	plt.plot(x,y, color='k', linestyle='dotted')

r = 0.15
d = 0.64
x = r*np.cos(theta) + d
y = r*np.sin(theta)
plt.plot(x,y, color='k')

r1 = rs[5]
r2 = rs[6]

xb = x[(x>=r1)&(x<=r2)&(y<0.)] 
yb = y[(x>=r1)&(x<=r2)&(y<0.)] 
xt = x[(x>=r1)&(x<=r2)&(y>=0.)] 
yt = y[(x>=r1)&(x<=r2)&(y>=0.)] 
xl = r1*np.cos(theta)
yl = r1*np.sin(theta)
ind = (yl<=yt[0])&(yl>=yb[0])&(xl>=0.)
xl = xl[ind]
yl = yl[ind]
ind0 = np.argsort(yl)
yl = yl[ind0]
xl = xl[ind0]
ind4 = np.argsort(yl)
ind4 = ind4[::-1]
yl = yl[ind4]
xl = xl[ind4]

xr = r2*np.cos(theta)
yr = r2*np.sin(theta)
ind1 = (yr<yt[-1])&(yr>yb[-1])&(xr>0.)
xr = xr[ind1]
yr = yr[ind1]
ind2 = np.argsort(yr)
yr = yr[ind2]
xr= xr[ind2]

xfill = np.append(xb,xr)
yfill = np.append(yb,yr)
xfill = np.append(xfill, xt)
yfill = np.append(yfill, yt)
xfill = np.append(xfill, xl)
yfill = np.append(yfill, yl)

for i in range(len(xfill)): print(xfill[i], yfill[i])
#for i in range(len(xr)): print(xr[i], yr[i])

plt.fill(xfill,yfill, color='orange')


plt.xlim((0., 1.0))
plt.ylim((-0.5, 0.5))
plt.axhline(0., color='0.5', zorder=-5, linestyle = 'dashed')

plt.show()
