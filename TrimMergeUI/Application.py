from __future__ import division, print_function

from os import pardir
from os.path import dirname, realpath, join, isfile

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import GLib, Gio, Gtk

from TrimMergeUI.MainWindow import MainWindow
from TrimMergeUI.AboutWindow import AboutWindow, _version
#from TrimMergeUI.PreferencesDialog import PreferencesDialog


class TrimMergeApplication(Gtk.Application):

    def __init__(self, *args, **kwargs):
        super(TrimMergeApplication, self).__init__(*args, application_id="org.projectx_bio.trimmerge",
                                                   **kwargs)
        self.window = None

        dir_path = join(dirname(realpath(__file__)), pardir)
        self.config_file_name = join(dir_path, 'config.ini')

    def do_startup(self):
        dir_path = join(dirname(realpath(__file__)), 'xml')
        menu_ui_file = join(dir_path, 'app_menu.glade')

        Gtk.Application.do_startup(self)

        action = Gio.SimpleAction.new("preferences", None)
        action.connect("activate", self.on_preferences)
        self.add_action(action)

        action = Gio.SimpleAction.new("about", None)
        action.connect("activate", self.on_about)
        self.add_action(action)

        action = Gio.SimpleAction.new("quit", None)
        action.connect("activate", self.on_quit)
        self.add_action(action)

        builder = Gtk.Builder.new_from_file(menu_ui_file)
        self.set_app_menu(builder.get_object("app-menu"))

    def do_activate(self):
        # We only allow a single window and raise any existing ones
        if not self.window:
            # Windows are associated with the application
            # when the last one is closed the application shuts down
            self.window = MainWindow(application=self, title="TrimMerge")
            self.window.connect("delete-event", self.on_quit)
        self.window.present()

    def on_about(self, action, param):
        about_dialog = AboutWindow(transient_for=self.window, modal=True)
        about_dialog.present()

    def on_preferences(self, action, param):
        print('Coming soon')
        return
        '''
        preferences_dialog = PreferencesDialog(self.window, config_file_name=self.config_file_name)
        preferences_dialog.present()
        response = preferences_dialog.run()
        if response == Gtk.ResponseType.OK:
            print("The OK button was clicked")
        elif response == Gtk.ResponseType.CANCEL:
            print("The Cancel button was clicked")
            preferences_dialog.destroy()
        '''

    def on_quit(self, action, param):
        print('Doing cleanup!')
        self.quit()
