# coding=utf-8
from __future__ import division, print_function
from os import pardir
from os.path import dirname, realpath, join
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import GLib, Gio, Gtk


_version = '0.0.1'


class AboutWindow(Gtk.AboutDialog):

    def __init__(self, *args, **kwargs):
        super(AboutWindow, self).__init__(*args, **kwargs)
        self.set_program_name('TrimMerge')
        self.set_version(_version)
        self.set_logo_icon_name('applications-utilities')

        dir_path = join(dirname(realpath(__file__)), pardir)
        licence_file = join(dir_path, 'LICENSE')
        try:
            f = open(licence_file)
            lic = f.read()
            f.close()
            self.set_license(lic)
            self.set_wrap_license(True)
            self.set_copyright('Copyright © 2016 Natalya Bondarenko & Anton Bondarenko')
        except IOError:
            self.set_copyright('Copyright © 2016 Natalya Bondarenko & Anton Bondarenko' +
                               '\nLicensed under the Apache License, Version 2.0')
            pass
        self.set_website('https://github.com/nilannik/TrimMerge')
        self.set_comments('Tool for adapter trimming and sequence assembly using Nextera pair-end libraries')
        self.add_credit_section('Created by', ['Natalya Bondarenko', 'Anton Bondarenko'])
        self.show_all()
