from __future__ import division, print_function
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk
from numpy import *

from matplotlib.figure import Figure
import matplotlib.cm as cm
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas

myfirstwindow = Gtk.Window()
myfirstwindow.connect("delete-event", Gtk.main_quit)
myfirstwindow.set_default_size(400, 400)

#fig = Figure(figsize=(5,5), dpi=100)
fig = Figure()
ax = fig.add_subplot(111, projection='polar')

N = 20
theta = linspace(0.0, 2 * pi, N, endpoint=False)
radii = 10 * random.rand(N)
width = pi / 4 * random.rand(N)

bars = ax.bar(theta, radii, width=width, bottom=0.0)

for r, bar in zip(radii, bars):
    bar.set_facecolor(cm.jet(r / 10.))
    bar.set_alpha(0.5)

ax.plot()

sw = Gtk.ScrolledWindow()
myfirstwindow.add(sw)

canvas = FigureCanvas(fig)
#canvas.set_size_request(400,400)
sw.add_with_viewport(canvas)

myfirstwindow.show_all()
Gtk.main()
