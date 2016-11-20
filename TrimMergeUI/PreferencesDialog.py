# coding=utf-8
from __future__ import division, print_function
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import GLib, Gio, Gtk

from ScientificProjects.Config import read_config, write_config


class PreferencesDialog(Gtk.Dialog):

    def __init__(self, parent, config_file_name):
        self.config_file_name = config_file_name
        Gtk.Dialog.__init__(self, 'Preferences', parent, 0,
                            (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                             Gtk.STOCK_OK, Gtk.ResponseType.OK))
        self.set_default_size(150, 100)
        self.set_border_width(6)
        label = Gtk.Label("Database settings")
        box = self.get_content_area()
        box.add(label)
        box.set_spacing(15)

        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=15)
        box.add(hbox)

        vbox_left = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
        vbox_right = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
        vbox_left.set_homogeneous(True)
        vbox_right.set_homogeneous(True)
        hbox.pack_start(vbox_left, False, True, 0)
        hbox.pack_start(vbox_right, True, True, 0)

        label = Gtk.Label("name")
        label.set_halign(Gtk.Align.END)
        self.db_name = Gtk.Entry()
        vbox_left.pack_start(label, True, True, 0)
        vbox_right.pack_start(self.db_name, True, True, 0)

        label = Gtk.Label("backend")
        label.set_halign(Gtk.Align.END)
        self.db_backend = Gtk.Entry()
        vbox_left.pack_start(label, True, True, 0)
        vbox_right.pack_start(self.db_backend, True, True, 0)

        label = Gtk.Label("host")
        label.set_halign(Gtk.Align.END)
        self.db_host = Gtk.Entry()
        vbox_left.pack_start(label, True, True, 0)
        vbox_right.pack_start(self.db_host, True, True, 0)

        label = Gtk.Label("port")
        label.set_halign(Gtk.Align.END)
        self.db_port = Gtk.SpinButton()
        self.db_port.set_digits(0)
        self.db_port.set_range(0, 100000)
        self.db_port.set_increments(1, 10);
        vbox_left.pack_start(label, True, True, 0)
        vbox_right.pack_start(self.db_port, True, True, 0)

        label = Gtk.Label("user")
        label.set_halign(Gtk.Align.END)
        self.db_user = Gtk.Entry()
        vbox_left.pack_start(label, True, True, 0)
        vbox_right.pack_start(self.db_user, True, True, 0)

        label = Gtk.Label("password")
        label.set_halign(Gtk.Align.END)
        self.db_password = Gtk.Entry()
        self.db_password.set_visibility(False)
        vbox_left.pack_start(label, True, True, 0)
        vbox_right.pack_start(self.db_password, True, True, 0)

        self.show_all()
        self.config = None
        self.read_in_config()

    def read_in_config(self):
        try:
            self.config = read_config(self.config_file_name)
            self.db_name.set_text(self.config['db_name'])
            self.db_backend.set_text(self.config['backend'])
            self.db_host.set_text(self.config['host'])
            port = self.config['port']
            try:
                port = int(port)
                self.config['port'] = port
            except ValueError:
                if self.config['port'] == '':
                    self.config['port'] = 0
                pass
            try:
                self.db_port.set_value(port)
            except TypeError:
                pass
            self.db_user.set_text(self.config['user'])
            self.db_password.set_text(self.config['password'])
        except (IOError, ValueError):
            self.config = None

    def save_config(self):
        config = {'db_name': self.db_name.get_text(), 'backend': self.db_backend.get_text(),
                  'host': self.db_host.get_text(), 'port': int(self.db_port.get_value()),
                  'user': self.db_user.get_text(), 'password': self.db_password.get_text()}
        if config != self.config:
            print('config changed')
            write_config(config, self.config_file_name)
            return True
        return False
