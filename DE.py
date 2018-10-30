import numpy as np
from matplotlib import pyplot as plt
from math import sin, cos, exp

	
def f(x,y):       								 #fuction of my equation
	return sin(x)*sin(x)+y*(cos(x)/sin(x))

def euler(x0,y0,xf,n):                           #function for Euler’s Method
	x = np.linspace(x0, xf, n + 1)
	h = x[1] - x[0]
	y = [y0]
	for i in range(len(x)-1):
		y.append(y[i]+h*f(x[i], y[i]))
	return x,y     

def imp_euler(x0,y0,xf,n):                       #function for Improved Euler’s Method
	x = np.linspace(x0, xf, n + 1)
	h = x[1] - x[0]
	y = [y0]
	for i in range(len(x)-1):
		m1=f(x[i], y[i])
		m2=f(x[i+1], y[i]+h*m1)
		y.append(y[i] + h*(m1+m2)/2)
	return x,y

def runge_kutta(x0,y0,xf,n):                     #function for Runge-Kutta Method
	x = np.linspace(x0, xf, n + 1)
	h = x[1] - x[0]
	y = [y0]
	for i in range(len(x)-1):
		k1 = h*f(x[i], y[i])
		k2 = h*f(x[i] + h/2, y[i] + k1/2)
		k3 = h*f(x[i] + h/2, y[i] + k2/2)
		k4 = h*f(x[i] + h, y[i] + k3)
		y.append(y[i] + (k1 + 2*k2 + 2*k3 + k4)/6)
	return x,y

def exact(x0,y0,xf,n):                           #function where i caluclate  with help of my exact solution
	x = np.linspace(x0, xf, n + 1)
	c = (y0+sin(x[0])*cos(x[0]))/sin(x[0])
	y = []
	for i in range(len(x)):
		y.append(c*sin(x[i])-sin(x[i])*cos(x[i]))
	return x,y

def calc_errors(x0,y0,xf,n1,n2):                 #function where i calculate global truncation error
	errors1 = [0]*(n2-n1)
	errors2 = [0]*(n2-n1)
	errors3 = [0]*(n2-n1)
	steps = [0]*(n2-n1)

	for i in range(0,n2-n1):
		steps[i]=n1+i
		x,y = exact(x0,y0,xf,steps[i])
		ey = y[-1]
		x,y = euler(x0,y0,xf,steps[i])
		errors1[i]=abs(ey-y[-1])
		x,y = imp_euler(x0,y0,xf,steps[i])
		errors2[i]=abs(ey-y[-1])
		x,y = runge_kutta(x0,y0,xf,steps[i])
		errors3[i]=abs(ey-y[-1])
	return steps,errors1,errors2,errors3



	
x0=float(input("x0 = "))
y0=float(input("y0 = "))
xmax=float(input("X = "))
steps=int(input("Steps = "))
n1=int(input("For Errors\nMin number of steps = "))
n2=int(input("Max number of steps = "))

x1,y1=euler(x0,y0,xmax,steps)
x2,y2=imp_euler(x0,y0,xmax,steps)
x3,y3=runge_kutta(x0,y0,xmax,steps)
x4,y4=exact(x0,y0,xmax,steps)

ax = plt.subplot(211)
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.3,
                 box.width, box.height*0.9])

plt.plot(x1,y1, 'r-',label="Euler’s Method")
plt.plot(x2,y2, 'b-',label="Improved Euler’s Method")
plt.plot(x3,y3, 'g-',label="Runge-Kutta Method")
plt.plot(x4,y4, 'k',label="Analysis Method")

plt.xlabel('Value of x')
plt.ylabel('Value of y')
plt.title('Approximate Solution')

n,e1,e2,e3=calc_errors(1,1,3,n1,n2+1)
bx = plt.subplot(212)
plt.plot(n,e1, 'r-',label="Euler’s Method")
plt.plot(n,e2, 'b-',label="Improved Euler’s Method")
plt.plot(n,e3, 'g-',label="Runge-Kutta Method")


box = bx.get_position()
bx.set_position([box.x0, box.y0 + box.height * 0.3,
                 box.width, box.height*0.9])

plt.xlabel('Number of steps')
plt.ylabel('Error')

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -1.6),
          fancybox=True, shadow=True, ncol=2)



plt .show()
