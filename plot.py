import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

directory = '/Users/chan/Desktop/data/'

N, J, t = open(directory + 'parameters.txt', 'r').readline().split()
N = int(N)
J = int(J)
t = float(t)

data = np.fromfile(directory + 'output', dtype='float64')
data.shape = (N + 1, J + 1)

fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(-1.5, 1.5))
x = np.linspace(0, 1, J + 1)
line, = ax.plot([], [], lw=2)

def animate(i):
    line.set_data(x, data[i])
    return line,

def init():
    line.set_data([], [])
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(0, N, N / t / 50), init_func=init, interval=1, repeat=False)

plt.show()
