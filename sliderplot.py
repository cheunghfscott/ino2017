import matplotlib.animation as animation
import matplotlib.pyplot as plt
import visual as v
import numpy as np
import time

plt.style.use('ggplot')

class My_mechanism(object):
  
    # The init function
    
    def __init__(self,a,b,alpha_dot,alpha0=0):
        self.a = a                            # Rod a
        self.b = b                            # Rod b
        self.alpha_dot = alpha_dot            # Angular speed of rod a in rad/s
        self.alpha0 = alpha0                  # Initial angular position of rod a
        self.k = np.array([0.01])             # Initial time for animation
        self.c_position = []                  # Piston position for animation
        self.conn_rod_angular_speed = []      # Classes for the graphs animation
        self.c_speed = []

    # Angular position of rod a as a function of time
    def alpha(self,t):
        if not all(t):
            raise ValueError('Each time t must be greater than 0')
        alpha = self.alpha0 + self.alpha_dot*t
        return alpha

    # Angular position of rod b as a function of time
    def beta(self,t):
        beta =  np.arcsin((self.a/self.b)*np.sin(self.alpha(t)))
        return beta

    # Piston position of rod a as a function of time
    def piston_position(self,t):
        c_y = np.zeros(t.shape[0])
        c_x = self.a*np.cos(self.alpha(t)) + self.b*np.cos(-np.arcsin((self.a/self.b)*np.sin(self.alpha(t))))
        return c_x,c_y

    # Rod a position in the complex plane
    def rod_a_position(self,t):
        a_y = self.a*np.sin(self.alpha(t))
        a_x = self.a*np.cos(self.alpha(t))
        return a_x,a_y

    # Angular speed of rod b
    def beta_dot(self,t):
        beta_dot = -(self.a*self.alpha_dot*np.cos(self.alpha(t)))/(self.b*np.cos(self.beta(t)))
        return beta_dot

    # Piston speed
    def c_dot(self,t):
        c_dot = self.a*self.alpha_dot*(np.tan(self.beta(t))*np.cos(self.alpha(t))-np.sin(self.alpha(t)))
        return c_dot

    # Animation using vpython
    def animation_vp(self,t0,tn,points):
        if (tn-t0) < 0:
            raise ValueError('tn must be greater than t0')
        if (tn < 0) or (tn < 0) or (points < 0):
            raise ValueError('tn, t0, and points must be greater than 0')

        t = np.linspace(t0,tn,num=points)
        beta = self.beta(t)
        a_x,a_y = self.rod_a_position(t)
        c_x,c_y = self.piston_position(t)

        crank = v.arrow(pos=(0,0,0),axis=(1,0,0),color=(1,1,1),length=m.a,make_trail=False)
        connectingRod = v.arrow(pos=(0,0,0),axis=(-1,0,0),color=(0,0,1),length=m.b,make_trail=False)
        piston = v.cylinder(pos=(c_x[0],c_y[0],0),radius=5,color=(1,0,0),length=10,make_trail=True)
        ball = v.sphere(pos=(0,0,0),radius=3,color=(1,1,1))

        theta0 = beta[1]

        for x_a,y_a,x_c,y_c,b in zip(a_x,a_y,c_x,c_y,beta):

            crank.axis = (x_a,y_a,0)
            connectingRod.pos = (x_a,y_a,0)
            ball.pos = (x_a,y_a,0)

            # Change in angular position
            dtheta = theta0-b
            theta0 = b
            connectingRod.rotate(angle=dtheta,axis=(0,0,1))

            connectingRod.pos = (x_c,y_c,0)
            piston.pos = (x_c,y_c,0)
            
            time.sleep(0.1)

    # Animation using matplotlib
    def animation_m(self):
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)
        
        def animate(i):
            a_x,a_y = self.rod_a_position(self.k)
            c_x,c_y = self.piston_position(self.k)
            ax1.clear()
            # Connecting rod (b)
            plt.xlim(-15,35)
            plt.ylim(-15,15)
            ax1.plot([0,a_x[0]],[0,a_y[0]],linewidth=3,color='blue')
            # Crankshaft rod (a)
            ax1.plot([a_x[0],c_x[0]],[a_y[0],c_y[0]],linewidth=3,color='green')
            # Piston (c)
            ax1.plot([0,c_x[0]],[0,c_y[0]],'o',markersize=20,color='red')
            
            self.k += 0.1
            
        ani = animation.FuncAnimation(fig,animate,interval=50)
        plt.show()

    # Matplotlib animation with graphs
    def animation_m_plus(self):
        fig = plt.figure()
        ax1 = fig.add_subplot(3,1,1)
        ax2 = fig.add_subplot(3,1,2)
        ax3 = fig.add_subplot(3,1,3)

        def animate(i):
            a_x,a_y = self.rod_a_position(self.k)
            c_x,c_y = self.piston_position(self.k)
            ax1.clear()
            ax2.clear()
            ax3.clear()
           
            # Crankshaft rod (a)
            ax1.plot([0,a_x[0]],[0,a_y[0]],linewidth=3,color='blue')
            # Connecting rod (b)
            ax1.plot([c_x[0],a_x[0]],[c_y[0],a_y[0]],linewidth=3,color='green')
            # Piston
            ax1.plot([0,c_x[0]],[0,c_y[0]],'o',markersize=20,color='red')
            ax1.set_xlim(-15,35)
            ax1.set_ylim(-11,11)
            ax1.set_title('Crankshaft, connecting rod and piston mechanism')
            # Piston position
            self.c_position.append(c_x)
            ax2.plot(self.c_position,color='green')
            ax2.set_xlim(0,300)
            ax2.set_ylim(8,32)
            ax2.set_title('Piston position')
            # Piston speed
            self.c_speed.append(self.c_dot(self.k))
            ax3.plot(self.c_speed,color='green')
            ax3.plot([0,300],[0,0],linewidth=1,color='black')
            ax3.set_xlim(0,300)
            ax3.set_ylim(-15,15)
            ax3.set_title('Piston speed')

            self.k += 0.1
            
        ani = animation.FuncAnimation(fig,animate,interval=50)
        plt.show()
