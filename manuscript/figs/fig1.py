import matplotlib.pyplot as plt
import numpy as np
import math
from pylab import *

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

plt.figure(figsize = (6,6))
theta = np.linspace(0, 2*math.pi, 1000)

r = 0.15
fac = 7.7e-2 
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

xmin = (d**2 + r1**2 - r**2)/(2.*d)
xmax = (d**2 + r2**2 - r**2)/(2.*d)

xb = x[(x>=xmin)&(x<xmax)&(y<0.)] 
yb = y[(x>=xmin)&(x<xmax)&(y<0.)] 
xt = x[(x>=xmin)&(x<xmax)&(y>0.)] 
yt = y[(x>=xmin)&(x<xmax)&(y>0.)] 
xl = r1*np.cos(theta)
yl = r1*np.sin(theta)
ind = (yl<yt[0])&(yl>yb[0])&(xl>0.)
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

"""plt.clf()
plt.plot(xb, yb, label="b")
plt.plot(xr, yr, label="r")
plt.plot(xt, yt, label="t")
plt.plot(xl, yl, label="l")
plt.legend()
plt.show()"""

xfill = np.append(xb,xr)
yfill = np.append(yb,yr)
xfill = np.append(xfill, xt)
yfill = np.append(yfill, yt)
xfill = np.append(xfill, xl)
yfill = np.append(yfill, yl)

#for i in range(len(xfill)): print(xfill[i], yfill[i])
#for i in range(len(xr)): print(xr[i], yr[i])

plt.fill(xfill,yfill, color='orange')


r = 1.0
x = r*np.cos(theta)
y = r*np.sin(theta)

plt.plot(x,y, color='k')

plt.xlim((-0.06, 1.0))
plt.ylim((-0.53, 0.53))

plt.axhline(0., color='k', zorder=5, linestyle = 'dashed')
plt.axvline(0., color='k', zorder=5, linestyle = 'dashed')

plt.plot([0., d],[0., 0.], color='k', linewidth=1.5, zorder=6)
plt.plot(0., 0., ms= 5, marker='o', color='k')
plt.plot(d, 0., ms= 5, marker='o', color='k')
plt.gca().text(0.3,0.04, "d", fontsize=16, weight='bold')
#plt.gca().annotate("", xy = (0.05, 0.03), xycoords='data', xytext =  (0.6, 0.03), arrowprops = dict(arrowstyle="->", fc="0.4", ec="0.4"))
#plt.gca().annotate("", xy = (0.05, 0.03), xycoords='data', xytext =  (0.6, 0.03), arrowprops = dict(arrowstyle="<-", fc="0.4", ec = "0.4"))

plt.xlabel("x (stellar radii)")
plt.ylabel("y (stellar radii)")

plt.savefig("f1.eps", dpi=300)
plt.show()
