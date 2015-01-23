import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
# from clawpack.visclaw.JSAnimation import IPython_display

x = np.linspace(0,1,1000)  # Spatial grid
t = np.linspace(0,1)       # Temporal grid
a = 1.0                    # Advection speed

def q_0(x):                # Initial condition
    return np.exp(-200.*(x-0.2)**2)
fig = plt.figure(figsize=(8,4))      # Create an empty figure
ax  = plt.axes()
line, = ax.plot([], [],linewidth=2)  # Create an empty line plot
plt.axis((0,1,-0.1,1.1))             # Set the bounds of the plot

def plot_q(t):
    line.set_data(x,q_0(x-a*t))  # Replace the line plot with the solution at time t
    
animation.FuncAnimation(fig, plot_q, frames=t)  # Animate the solution
