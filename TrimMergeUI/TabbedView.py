from __future__ import division, print_function

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

from TrimMergeUI.TrimPage import TrimPage


class TabbedView(Gtk.Notebook):

    def __init__(self):
        super(TabbedView, self).__init__()
        self.trim_page = TrimPage()
        self.trim_page.set_border_width(10)
        self.append_page(self.trim_page, Gtk.Label('Gentle Trim'))

        self.merge_page = Gtk.Box()
        self.merge_page.set_border_width(10)
        self.merge_page.add(Gtk.Label('Coming soon...'))
        self.append_page(self.merge_page, Gtk.Label('Merge'))
        self.trim_page.show()
        self.merge_page.show()