from __future__ import division, print_function

from os import pardir
from os.path import join, dirname, basename, realpath, splitext

try:
    import itertools.izip as zip
except ImportError:
    pass

from Bio import SeqIO
from Bio.Alphabet import IUPAC

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

from TrimMergeUI.utils import find_insert_in_record, clean_records


class TrimPage(Gtk.Box):

    def __init__(self):
        super(TrimPage, self).__init__()

        grid = Gtk.Grid(column_spacing=10, row_spacing=10)
        self.add(grid)

        self.select_fr_button = Gtk.Button("Choose FR File")
        self.select_fr_button.connect("clicked", self.on_file_clicked)
        self.fr_file_name = None
        self.fr_file_label = Gtk.Label('Please choose FR file   ')

        self.select_rf_button = Gtk.Button("Choose RF File")
        self.select_rf_button.connect("clicked", self.on_file_clicked)
        self.rf_file_name = None
        self.rf_file_label = Gtk.Label('Please choose RF file   ')

        file_formats = ["fastq", "fasta", ]
        self.file_format_combo = Gtk.ComboBoxText()
        for file_format in file_formats:
            self.file_format_combo.append_text(file_format)
        self.file_format_combo.set_active(0)

        alphabets = ["Ambiguous DNA", "Unambiguous DNA", "Extended DNA"]
        self.alphabet_combo = Gtk.ComboBoxText()
        for alphabet in alphabets:
            self.alphabet_combo.append_text(alphabet)
        self.alphabet_combo.set_active(0)

        self.select_out_dir_button = Gtk.Button("Choose output directory")
        self.select_out_dir_button.connect("clicked", self.on_folder_clicked)
        self.output_dir_name = None
        self.output_dir_label = Gtk.Label('Please choose output dir')

        self.adapters_dict = None
        adapters = ["Nextera PE",]
        self.adapter_combo = Gtk.ComboBoxText()
        self.adapter_combo.connect("changed", self.on_adapter_combo_changed)
        for adapter in adapters:
            self.adapter_combo.append_text(adapter)
        self.adapter_combo.set_active(0)

        self.adapter_min_length = 20
        adjustment = Gtk.Adjustment(self.adapter_min_length, 3, 100, 1, 10, 0)
        self.adapter_min_length_button = Gtk.SpinButton()
        self.adapter_min_length_button.set_adjustment(adjustment)
        self.adapter_min_length_button.set_numeric(True)
        self.adapter_min_length_button.set_update_policy(Gtk.SpinButtonUpdatePolicy.IF_VALID)
        self.adapter_min_length_button.connect("changed", self.on_adapter_min_length_changed)
        self.adapter_min_length_button.set_value(self.adapter_min_length)

        self.adapter_similarity = 60
        adjustment = Gtk.Adjustment(self.adapter_similarity, 10, 100, 1, 10, 0)
        self.adapter_similarity_button = Gtk.SpinButton()
        self.adapter_similarity_button.set_adjustment(adjustment)
        self.adapter_similarity_button.set_numeric(True)
        self.adapter_similarity_button.set_update_policy(Gtk.SpinButtonUpdatePolicy.IF_VALID)
        self.adapter_similarity_button.connect("changed", self.on_adapter_similarity_changed)
        self.adapter_similarity_button.set_value(self.adapter_similarity)

        self.run_button = Gtk.Button("Run")
        self.run_button.connect("clicked", self.on_run_clicked)

        self.status_label = Gtk.Label('Idle')

        self.progressbar = Gtk.ProgressBar()
        self.progressbar.set_fraction(1.0)

        grid.attach(self.select_fr_button, 0, 0, 1, 1)
        grid.attach(self.fr_file_label, 1, 0, 1, 1)
        grid.attach(self.select_rf_button, 0, 1, 1, 1)
        grid.attach(self.rf_file_label, 1, 1, 1, 1)

        grid.attach(Gtk.Label('Select file format:'), 0, 2, 1, 1)
        grid.attach(self.file_format_combo, 1, 2, 1, 1)

        grid.attach(Gtk.Label('Select alphabet:'), 0, 3, 1, 1)
        grid.attach(self.alphabet_combo, 1, 3, 1, 1)

        grid.attach(self.select_out_dir_button, 0, 4, 1, 1)
        grid.attach(self.output_dir_label, 1, 4, 1, 1)
        grid.attach(Gtk.Label('Select adapter:'), 0, 5, 1, 1)
        grid.attach(self.adapter_combo, 1, 5, 1, 1)
        grid.attach(Gtk.Label('Min adapter length:'), 0, 6, 1, 1)
        grid.attach(self.adapter_min_length_button, 1, 6, 1, 1)
        grid.attach(Gtk.Label('Adapter similarity %:'), 0, 7, 1, 1)
        grid.attach(self.adapter_similarity_button, 1, 7, 1, 1)
        grid.attach(self.run_button, 0, 8, 2, 1)
        grid.attach(self.progressbar, 0, 9, 2, 1)
        grid.attach(self.status_label, 0, 10, 2, 1)

    def on_file_clicked(self, widget):
        dialog = Gtk.FileChooserDialog("Please choose a file", self.get_toplevel(),
                                       Gtk.FileChooserAction.OPEN,
                                       (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                        Gtk.STOCK_OPEN, Gtk.ResponseType.OK))

        self.add_filters(dialog)

        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            print("Open clicked")
            print("File selected: " + dialog.get_filename())
            if widget == self.select_fr_button:
                self.fr_file_name = dialog.get_filename()
                self.file_name_label_update(self.fr_file_label, self.fr_file_name)
            elif widget == self.select_rf_button:
                self.rf_file_name = dialog.get_filename()
                self.file_name_label_update(self.rf_file_label, self.rf_file_name)
        elif response == Gtk.ResponseType.CANCEL:
            print("Cancel clicked")

        dialog.destroy()

    def add_filters(self, dialog):
        filter_text = Gtk.FileFilter()
        filter_text.set_name("Text files")
        filter_text.add_mime_type("text/plain")
        dialog.add_filter(filter_text)

        filter_any = Gtk.FileFilter()
        filter_any.set_name("Any files")
        filter_any.add_pattern("*")
        dialog.add_filter(filter_any)

    def on_folder_clicked(self, widget):
        dialog = Gtk.FileChooserDialog("Please choose a folder", self.get_toplevel(),
                                       Gtk.FileChooserAction.SELECT_FOLDER,
                                       (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                        "Select", Gtk.ResponseType.OK))
        dialog.set_default_size(800, 400)

        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            print("Select clicked")
            print("Folder selected: " + dialog.get_filename())
            self.output_dir_name = dialog.get_filename()
            self.file_name_label_update(self.output_dir_label, self.output_dir_name)
        elif response == Gtk.ResponseType.CANCEL:
            print("Cancel clicked")

        dialog.destroy()

    def on_adapter_combo_changed(self, combo):
        text = combo.get_active_text()
        if text != None:
            print("Selected adapter=%s" % text)
            self.read_in_adapters(text)

    def on_adapter_min_length_changed(self, widget):
        value = int(widget.get_value())
        print("Adapter min length changed=%d" % value)
        self.adapter_min_length = value

    def on_adapter_similarity_changed(self, widget):
        value = int(widget.get_value())
        print("Adapter similarity changed=%d" % value)
        self.adapter_similarity = value

    def file_name_label_update(self, label, file_name):
        if len(file_name) > 24:
            short_name = '...' + file_name[-21:]
        else:
            short_name = file_name
        label.set_text(short_name)

    def read_in_adapters(self, adapter_name):
        dir_path = join(dirname(realpath(__file__)), pardir)
        adapters_dir = join(dir_path, 'adapters')
        print(adapters_dir)
        if adapter_name == 'Nextera PE':
            adapters_file_name = join(adapters_dir, 'NexteraPE-PE_2.txt')
            print(adapters_file_name)
            nextera_records = SeqIO.parse(adapters_file_name, 'fasta', alphabet=IUPAC.ambiguous_dna)
            self.adapters_dict = dict()
            for record in nextera_records:
                if record.id == 'Trans2_rc':
                    self.adapters_dict['FR'] = record
                elif record.id == 'Trans1_rc':
                    self.adapters_dict['RF'] = record
        print(self.adapters_dict)

    def on_run_clicked(self, widget):
        self.status_label.set_text('Starting...')
        file_format = self.file_format_combo.get_active_text()
        alphabet = self.alphabet_combo.get_active_text()
        if alphabet == "Ambiguous DNA":
            alphabet = IUPAC.ambiguous_dna
        elif alphabet ==  "Unambiguous DNA":
            alphabet = IUPAC.unambiguous_dna
        elif alphabet == "Extended DNA":
            alphabet = IUPAC.extended_dna
        self.status_label.set_text('Reading in sequences...')
        records_FR = SeqIO.parse(self.fr_file_name, file_format, alphabet=alphabet)
        max_count = 0
        for record in records_FR:
            max_count += 1
        print(max_count)

        records_FR = SeqIO.parse(self.fr_file_name, file_format, alphabet=alphabet)
        records_RF = SeqIO.parse(self.rf_file_name, file_format, alphabet=alphabet)
        clean_FR, clean_RF, bad_FR, bad_RF = self.clean_PE_reads(records_FR, records_RF, max_count)

        fr_base_name, fr_extension = splitext(basename(self.fr_file_name))
        rf_base_name, rf_extension = splitext(basename(self.rf_file_name))

        file_out_FR_clean = join(self.output_dir_name, fr_base_name + '_clean' + fr_extension)
        file_out_FR_bad = join(self.output_dir_name, fr_base_name + '_bad' + fr_extension)
        file_out_RF_clean = join(self.output_dir_name, rf_base_name + '_clean' + fr_extension)
        file_out_RF_bad = join(self.output_dir_name, rf_base_name + '_bad' + fr_extension)

        SeqIO.write(clean_FR, file_out_FR_clean, file_format)
        SeqIO.write(clean_RF, file_out_RF_clean, file_format)
        SeqIO.write(bad_FR, file_out_FR_bad, file_format)
        SeqIO.write(bad_RF, file_out_RF_bad, file_format)
        self.status_label.set_text('Idle')


    def clean_PE_reads(self, records_FR, records_RF, max_count=1e10):
        clean_FR = []
        clean_RF = []

        bad_FR = []
        bad_RF = []

        count = 0
        short_read_threshold = 10
        for record_FR, record_RF in zip(records_FR, records_RF):
            clean = False
            self.progressbar.set_fraction(count / max_count)
            while Gtk.events_pending():
                Gtk.main_iteration()

            clean_FR_rec, clean_RF_rec, bad_FR_rec, bad_RF_rec = clean_records(record_FR, record_RF,
                                                                               self.adapters_dict,
                                                                               self.adapter_min_length,
                                                                               self.adapter_similarity,
                                                                               short_read_threshold)
            if clean_FR_rec:
                clean_FR.append(clean_FR_rec)
            if clean_RF_rec:
                clean_RF.append(clean_RF_rec)
            if bad_FR_rec:
                bad_FR.append(bad_FR_rec)
            if bad_RF_rec:
                bad_RF.append(bad_RF_rec)

            count += 1

        return clean_FR, clean_RF, bad_FR, bad_RF
