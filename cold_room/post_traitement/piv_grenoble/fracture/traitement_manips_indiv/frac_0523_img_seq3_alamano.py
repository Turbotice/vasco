#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#%%
# ici on va juste voir le resultat que ca donne avec image j, en regardant l'amplitude 
# à plusieurs points pour estimer la courbure
frame_frac = 4603

frame1 = 4572
frame2 = 4582

x1_vals = np.array([1152,833,587,1289,1053,724,942,791,884,848,823,784])
x2_vals = np.array([1154,837,586,1290,1053,724,944,789,884,849,823,783])

y1_vals = np.array([283,285,291,289,288,280,295,310,312,283,306,304])
y2_vals = np.array([271,264,277,281,273,261,276,290,292,262,285,284])

x2fit = (x1_vals+x2_vals)/2
y2fit = (y1_vals-y2_vals)/2

indices2fit = np.where((x2fit>500)&(x2fit<1200))

a,b,c = np.polyfit(x2_vals[indices2fit],y2fit[indices2fit],2)


xth=np.linspace(np.min(x2fit),np.max(x2fit))
plt.figure()
plt.plot((x1_vals+x2_vals)/2, (y1_vals-y2_vals)/2,'o',label='(y1_vals-y2_vals)/2')
plt.plot(xth,a*xth**2+b*xth+c)
plt.ylim(0,11)
plt.legend()
plt.show()


dcmperdpx = 0.065
xth=np.linspace(np.min(x2fit)*dcmperdpx,np.max(x2fit)*dcmperdpx) * 1e-2

a,b,c = np.polyfit(x2_vals[indices2fit] * dcmperdpx * 1e-2,y2fit[indices2fit] * dcmperdpx * 1e-2,2)

plt.figure()
plt.plot(1e-2 * dcmperdpx * (x1_vals+x2_vals)/2,1e-2 * dcmperdpx * (y1_vals-y2_vals)/2,'o',label='(y1_vals-y2_vals)/2')
plt.plot(xth,a*xth**2+b*xth+c)
#plt.ylim(0,11*1e-2)
popt,pcov = curve_fit((lambda x,A,phi:A*np.cos((2*np.pi/1.4)*x+phi)),x2_vals * dcmperdpx * 1e-2,y2fit * dcmperdpx * 1e-2)
plt.plot(xth,popt[0]*np.cos((2*np.pi/1.4)*xth+popt[1]))
plt.legend()
plt.show()



kappa = 2*a
print(kappa)
# %%