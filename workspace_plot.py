import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#sand out: upward outter spiral
t = np.linspace(-200,200,400)
x = 100*np.sin(t/5) + (100 + 224)
y = 100*np.cos(t/5) + (100 + 37)
z = t

#salut out
t = np.linspace(-200,200,400)
x = t + (100 + 224)
y = t + (200 + 37)
z = t

#forward extension
t = np.linspace(100,400,31)
x = t
y = 0*t+350
z = 0*t-150

# curve
t = np.linspace(-np.pi,np.pi,10)
x = 0*t + 100
y = -250 - 250*np.cos(t)
z = 80 + 250*np.sin(t)

# curve
t = np.linspace(-np.pi,np.pi,10)
x = 0*t + 100
y = -250 - 100*np.cos(t)
z = 80 + 100*np.sin(t)

#try:
#x = ln(t), y = sin(t), z = tcos(t) forward-facing cone
#x = sin(t)cos(2.3t), y = sin(t)sin(2.3t), z = cos(t) ball of yarn

ax = plt.axes(projection = '3d')
ax.view_init(25,165)

ax.plot(x,y,z)
plt.show()