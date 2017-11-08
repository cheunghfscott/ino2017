"""
===========================
The double pendulum problem
===========================

This animation illustrates the double pendulum problem.

Double pendulum formula translated from the C code at
http://www.physics.usyd.edu.au/~wheat/dpend_html/solve_dpend.c
"""

from numpy import sin, cos, arctan
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation


def derivs(state, t):
    dydx = np.zeros_like(state)
    dydx[0] = state[1]

    del_ = state[2] - state[0]
    den1 = (M1 + M2) * L1 - M2 * L1 * cos(del_) * cos(del_)
    dydx[1] = (M2 * L1 * state[1] * state[1] * sin(del_) * cos(del_) +
               M2 * G * sin(state[2]) * cos(del_) +
               M2 * L2 * state[3] * state[3] * sin(del_) -
               (M1 + M2) * G * sin(state[0])) / den1

    dydx[2] = state[3]

    den2 = (L2 / L1) * den1
    dydx[3] = (-M2 * L2 * state[3] * state[3] * sin(del_) * cos(del_) +
               (M1 + M2) * G * sin(state[0]) * cos(del_) -
               (M1 + M2) * L1 * state[1] * state[1] * sin(del_) -
               (M1 + M2) * G * sin(state[2])) / den2

    return dydx


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

    G = 9.8  # acceleration due to gravity, in m/s^2
    L1 = 0.3  # length of pendulum 1 in m
    L2 = 1.5  # length of pendulum 2 in m
    M1 = 1.0  # mass of pendulum 1 in kg
    M2 = 1.0  # mass of pendulum 2 in kg

    # create a time array from 0..100 sampled at 0.05 second steps
    dt = 0.05
    td=0.6
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
    time = np.linspace(0.0, 10.0, Num)
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
    r1[0]=1.7
    #constant angular speed
    th2dot=-1


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
    print(y2)
    # print(x2)
    # print(y2)
    fig = plt.figure(1)
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-3, 3), ylim=(-3, 3))
    ax.set_aspect('equal')
    ax.grid()

    line, = ax.plot([], [], 'o-', lw=2)
    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    # ani.save('double_pendulum.mp4', fps=15)
    ani = animation.FuncAnimation(fig, animate, np.arange(1, Num),
                                  interval=100, blit=True, init_func=init)

    plt.figure(2)
    plt.plot(time, x2)

    plt.show()
