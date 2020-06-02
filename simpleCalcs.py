import matplotlib.pyplot as plt
import numpy as np

# 100 linearly spaced numbers
x = np.linspace(0,10,200)
y0 = x/(1.+1.31e-2*x)
y1 = x/(1.+1.29e-2*x+9.59e-6*x*x)

# setting the axes at the centre
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# plot the function
plt.plot(x,x)
plt.plot(x,y0, 'r')
plt.plot(x,y1, 'b')

# show the plot
plt.show()

geomSet=[0,1]
bgTrgRate=[50.0,70.0] # kHz
Len=150.0
nSeg=[48,64]
sW=[55.0,42.6]
cW=[55.0,46.8]
a2=np.zeros(2)
area=np.zeros(2)
scale=1.0
for i in range(2):
	area[i]=nSeg[i]*(sW[i]+cW[i])*Len*4.0
	if(i==0):scale=area[i]
	a2[i]=area[i]/scale
	print(area[i])
	if(i>0):
		xx=sW[i]*nSeg[i]/(sW[0]*nSeg[0])
		yy=cW[i]*nSeg[i]/(cW[0]*nSeg[0])
		print(xx," x ",yy," = ",xx*yy," and ^2=",xx*xx*yy*yy)
		a2[i]=xx*yy
plt.plot(geomSet,area)
plt.xlabel('geom Id')
plt.ylabel('Area [a.u.]')
plt.show()

normRate=bgTrgRate/(a2*a2)

plt.plot(geomSet,bgTrgRate)
plt.plot(geomSet,normRate)
plt.xlabel('geom Id')
plt.ylabel('Normalised Trig.Rate[kHz]')
plt.show()
