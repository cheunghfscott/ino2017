"""
===========================
The slidercrank problem
===========================

"""

from numpy import sin, cos, arctan
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i * dt))
    return line, time_text


if __name__ == "__main__":


    L1 = 15  # length of pendulum 1 in m
    L2 = 21.45  # length of pendulum 2 in m


    # create a time array from 0..100 sampled at 0.05 second steps
    dt = 0.05
    td=30.0
    t = np.arange(0.0, td, dt)

    # th1 and th2 are the initial angles (degrees)
    # w10 and w20 are the initial angular velocities (degrees per second)
    th1 = 45.0
    w1 = 0.0
    th2 = -30.0
    w2 = 0.0
    # r2,r3 are linkage length
    r2 = L1
    r3 = L2

    # initial state
    # state = np.radians([th1, w1, th2, w2])
    Num = int(td/dt)
    # integrate your ODE using scipy.integrate.
    # y = integrate.odeint(derivs, state, t)

    # define variables
    time = np.linspace(0.0, td, Num)
    x1 = np.zeros(Num)
    y1 = np.zeros(Num)
    x2 = np.zeros(Num)
    y2 = np.zeros(Num)
    th2t = np.zeros(Num)
    #th2dot = np.zeros(Num)
    th3 = np.zeros(Num)
    th3dot = np.zeros(Num)

    # input as r1 as slider length
    #r1 = np.linspace(1.7, 1.8, Num)
    r1 = np.zeros(Num+1)
    r1[0]=10
    #constant angular speed
    th2dot=-0.05


    for i in range(0, Num):
        # position equations
        A = -2 * r1[i] * r2 * cos(0)
        # B=-2*r1*r2*sin(0)
        C = r1[i] ** 2 + r2 ** 2 - r3 ** 2
        tt = ((-C ** 2 + A ** 2) ** (1.0 / 2.0)) / (C - A)
        th2t[i] = 2 * arctan(tt)

        th3[i] = arctan((-r2 * sin(th2t[i]) / (r1[i] - r2 * cos(th2t[i]))))
        #velocity equations
        th3dot[i] = (-r2*cos(th2t[i])*th2dot)/(r3*cos(th3[i]))
        xedot = th2dot*r2*(-sin(th2t[i])+sin(th3[i])*cos(th2t[i])/cos(th3[i]))

        #ODE solve with first order

        xe= r1[i] + xedot*dt
        r1[i+1]=xe


        # print(th2d)
        # print (len(y))
        # print(y.shape)


        # joint locations

        x1[i] = L1 * cos(th2t[i])
        y1[i] = L1 * sin(th2t[i])

        x2[i] = L2 * cos(th3[i]) + x1[i]
        y2[i] = L2 * sin(th3[i]) + y1[i]
        #print(x2[i])
    th2d = th2t * (180.0 / np.pi)
    y = th2d
   # print(y2)
    # print(x2)
    # print(y2)
    fig = plt.figure(1)
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-10,30), ylim=(-1,20))
    ax.set_aspect('equal')
    ax.grid()

    line, = ax.plot([], [], 'o-', lw=2)
    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

  #  ani.save('slidercrank', fps=15)
    ani = animation.FuncAnimation(fig, animate, np.arange(1, Num),
                                  interval=10, blit=True, init_func=init)
    print(time)
    plt.figure(2)
    plt.plot(time, y)
    plt.suptitle("output angle plot ")
    plt.xlabel("time")
    plt.ylabel("theta in degree" )
    

    plt.figure(3)
    plt.plot(time,x2)
    plt.suptitle("input distance")
    plt.xlabel("time")
    plt.ylabel("Input distance" )

    plt.show()
