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

from TrimMergeUI.Worker import count_length, compare_reads, init_trim_worker, clean_reads

from multiprocessing import Pool


class TrimPage(Gtk.Box):
    def __init__(self):
        super(TrimPage, self).__init__()

        self.run = False
        self.total_length_of_reads_count = 0
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=25)
        self.add(vbox)
        self.lock_label = Gtk.Label('Locked while running')
        vbox.pack_start(self.lock_label, True, True, 0)
        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=25)
        vbox.pack_start(hbox, True, True, 0)

        grid = Gtk.Grid(column_spacing=10, row_spacing=10)
        hbox.pack_start(grid, True, True, 0)

        statistics_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
        hbox.pack_start(statistics_box, True, True, 0)
        stat_title = Gtk.Label()
        stat_title.set_markup('_'*50 + '\n<big><b>Input reads statistics</b></big>')
        statistics_box.pack_start(stat_title, False, False, 0)
        input_stat_grid = Gtk.Grid(column_spacing=10, row_spacing=10)
        statistics_box.pack_start(input_stat_grid, True, True, 0)
        self.total_number_of_reads = Gtk.Label('n/a')
        self.total_length_of_reads = Gtk.Label('n/a')
        input_stat_grid.attach(Gtk.Label('Reads count'), 0, 0, 1, 1)
        input_stat_grid.attach(self.total_number_of_reads, 1, 0, 1, 1)
        input_stat_grid.attach(Gtk.Label('Total length'), 0, 1, 1, 1)
        input_stat_grid.attach(self.total_length_of_reads, 1, 1, 1, 1)
        stat_title = Gtk.Label()
        stat_title.set_markup('_'*50 + '\n<big><b>Output reads statistics</b></big>')
        statistics_box.pack_start(stat_title, False, False, 0)
        output_stat_grid = Gtk.Grid(column_spacing=10, row_spacing=10)
        statistics_box.pack_start(output_stat_grid, True, True, 0)
        self.total_number_of_reads_out = Gtk.Label('n/a')
        self.total_length_of_reads_out = Gtk.Label('n/a')
        self.total_number_of_reads_with_adapters = Gtk.Label('n/a')
        self.total_number_of_suspicious = Gtk.Label('n/a')
        self.max_number_of_adapters = Gtk.Label('n/a')
        self.total_number_of_short = Gtk.Label('n/a')
        output_stat_grid.attach(Gtk.Label('Reads count'), 0, 0, 1, 1)
        output_stat_grid.attach(self.total_number_of_reads_out, 1, 0, 1, 1)
        output_stat_grid.attach(Gtk.Label('Total length'), 0, 1, 1, 1)
        output_stat_grid.attach(self.total_length_of_reads_out, 1, 1, 1, 1)
        output_stat_grid.attach(Gtk.Label('Reads with adapters'), 0, 2, 1, 1)
        output_stat_grid.attach(self.total_number_of_reads_with_adapters, 1, 2, 1, 1)
        output_stat_grid.attach(Gtk.Label('Suspicious reads count'), 0, 3, 1, 1)
        output_stat_grid.attach(self.total_number_of_suspicious, 1, 3, 1, 1)
        output_stat_grid.attach(Gtk.Label('Max adapters in one read'), 0, 4, 1, 1)
        output_stat_grid.attach(self.max_number_of_adapters, 1, 4, 1, 1)
        output_stat_grid.attach(Gtk.Label('Short reads count'), 0, 5, 1, 1)
        output_stat_grid.attach(self.total_number_of_short, 1, 5, 1, 1)

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
        adapters = ["Nextera PE", ]
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

        self.run_button = Gtk.ToggleButton("Run")
        self.run_button.connect("toggled", self.on_run_toggled)

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

        self.file_format = None
        self.alphabet = None

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
        self.show_all()
        self.lock_label.hide()

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
                self._file_name_label_update(self.fr_file_label, self.fr_file_name)
            elif widget == self.select_rf_button:
                self.rf_file_name = dialog.get_filename()
                self._file_name_label_update(self.rf_file_label, self.rf_file_name)
        elif response == Gtk.ResponseType.CANCEL:
            print("Cancel clicked")

        dialog.destroy()

    def add_filters(self, dialog):
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
            self._file_name_label_update(self.output_dir_label, self.output_dir_name)
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
        max_len = min(len(self.adapters_dict['FR']), len(self.adapters_dict['RF']))
        print('max_len', max_len)
        if value > max_len:
            self.adapter_min_length_button.set_value(max_len)
        else:
            print("Adapter min length changed=%d" % value)
            self.adapter_min_length = value

    def on_adapter_similarity_changed(self, widget):
        value = int(widget.get_value())
        print("Adapter similarity changed=%d" % value)
        self.adapter_similarity = value

    def _file_name_label_update(self, label, file_name):
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

    def on_run_toggled(self, button):
        if button.get_active():
            if self.output_dir_name is None:
                print('Select output dir first!')
                self.status_label.set_text('Select output dir first!')
                self.run_button.set_active(False)
                return
            state = "on"
            self.run = True
        else:
            state = "off"
            self.run = False
        print("Run button was turned", state)
        self.lock_controls(lock=self.run)
        if self.run:
            while Gtk.events_pending():
                Gtk.main_iteration()
            self.prepare_to_run()

    def lock_controls(self, lock=True):
        if lock:
            release = False
            self.lock_label.show()
        else:
            release = True
            self.lock_label.hide()
        self.select_fr_button.set_sensitive(release)
        self.select_rf_button.set_sensitive(release)
        self.file_format_combo.set_sensitive(release)
        self.alphabet_combo.set_sensitive(release)
        self.select_out_dir_button.set_sensitive(release)
        self.adapter_combo.set_sensitive(release)
        self.adapter_min_length_button.set_sensitive(release)
        self.adapter_similarity_button.set_sensitive(release)
        self.get_parent().set_show_tabs(release)

    def prepare_to_run(self):
        self.status_label.set_text('Starting...')
        self.file_format = self.file_format_combo.get_active_text()
        self.alphabet = self.alphabet_combo.get_active_text()
        if self.alphabet == "Ambiguous DNA":
            self.alphabet = IUPAC.ambiguous_dna
        elif self.alphabet == "Unambiguous DNA":
            self.alphabet = IUPAC.unambiguous_dna
        elif self.alphabet == "Extended DNA":
            self.alphabet = IUPAC.extended_dna

        response, max_count = self.count_records()
        if not response:
            self.run_button.set_active(False)
            return

        response = self.compare_records()
        if not response:
            self.run_button.set_active(False)
            return

        clean_FR, clean_RF, bad_FR, bad_RF, short_FR, short_RF = self.clean_PE_reads(max_count)

        fr_base_name, fr_extension = splitext(basename(self.fr_file_name))
        rf_base_name, rf_extension = splitext(basename(self.rf_file_name))

        file_out_FR_clean = join(self.output_dir_name, fr_base_name + '_clean' + fr_extension)
        file_out_FR_bad = join(self.output_dir_name, fr_base_name + '_suspicious' + fr_extension)
        file_out_FR_short = join(self.output_dir_name, fr_base_name + '_short' + fr_extension)
        file_out_RF_clean = join(self.output_dir_name, rf_base_name + '_clean' + fr_extension)
        file_out_RF_bad = join(self.output_dir_name, rf_base_name + '_suspicious' + fr_extension)
        file_out_RF_short = join(self.output_dir_name, rf_base_name + '_short' + fr_extension)

        SeqIO.write(clean_FR, file_out_FR_clean, self.file_format)
        SeqIO.write(clean_RF, file_out_RF_clean, self.file_format)
        SeqIO.write(bad_FR, file_out_FR_bad, self.file_format)
        SeqIO.write(bad_RF, file_out_RF_bad, self.file_format)
        SeqIO.write(short_FR, file_out_FR_short, self.file_format)
        SeqIO.write(short_RF, file_out_RF_short, self.file_format)

        self.run_button.set_active(False)
        self.status_label.set_text('Idle')

    def count_records(self):
        self.status_label.set_text('Reading in sequences...')
        records_FR = SeqIO.parse(self.fr_file_name, self.file_format, alphabet=self.alphabet)
        records_RF = SeqIO.parse(self.rf_file_name, self.file_format, alphabet=self.alphabet)

        pool = Pool()
        len_fr_iter = pool.imap(count_length, records_FR, chunksize=100)
        len_rf_iter = pool.imap(count_length, records_RF, chunksize=100)
        fr_done = False
        rf_done = False
        max_count_fr = 0
        max_count_rf = 0
        len_fr = 0
        len_rf = 0
        while not (fr_done and rf_done):
            while Gtk.events_pending():
                Gtk.main_iteration()
            try:
                len_fr += next(len_fr_iter)
                max_count_fr +=1
            except StopIteration:
                fr_done = True
                pass
            try:
                len_rf += next(len_rf_iter)
                max_count_rf +=1
            except StopIteration:
                rf_done = True
                pass
            self.total_number_of_reads.set_text('%d' % max_count_fr)
            self.total_length_of_reads_count = len_fr + len_rf
            self.total_length_of_reads.set_text('%g' % self.total_length_of_reads_count)
        if max_count_fr != max_count_rf:
            dialog = Gtk.MessageDialog(self.get_toplevel(), 0, Gtk.MessageType.ERROR,
                                       Gtk.ButtonsType.CANCEL, "Different length of FR and RF files")
            dialog.format_secondary_text(
                "Number of forward and reverse reads must match")
            dialog.run()
            dialog.destroy()
            result = False, 0
        result = True
        self.status_label = Gtk.Label('Idle')
        return result, max_count_fr

    def compare_records(self):
        self.status_label.set_text('Checking sequences...')
        records_FR = SeqIO.parse(self.fr_file_name, self.file_format, alphabet=self.alphabet)
        records_RF = SeqIO.parse(self.rf_file_name, self.file_format, alphabet=self.alphabet)

        pool = Pool()
        len_fr_iter = pool.imap(compare_reads, zip(records_FR, records_RF), chunksize=100)
        cmp_done = False
        result = True
        while not cmp_done:
            while Gtk.events_pending():
                Gtk.main_iteration()
            try:
                result = next(len_fr_iter)
                if not result:
                    pool.terminate()
                    break
            except StopIteration:
                cmp_done = True
                pass

        if not result:
            dialog = Gtk.MessageDialog(self.get_toplevel(), 0, Gtk.MessageType.ERROR,
                                       Gtk.ButtonsType.CANCEL, "Different IDs in FR and RF files")
            dialog.format_secondary_text(
                "IDs of PE reads must match")
            dialog.run()
            dialog.destroy()
        self.status_label = Gtk.Label('Idle')
        return result

    def clean_PE_reads(self, max_count=1e10):
        self.status_label.set_text('Scanning for adapters...')
        records_FR = SeqIO.parse(self.fr_file_name, self.file_format, alphabet=self.alphabet)
        records_RF = SeqIO.parse(self.rf_file_name, self.file_format, alphabet=self.alphabet)
        short_read_threshold = self.adapter_min_length

        pool = Pool(initializer=init_trim_worker, initargs=(self.adapters_dict, self.adapter_min_length,
                                                            self.adapter_similarity,
                                                            short_read_threshold))
        results = pool.imap(clean_reads, zip(records_FR, records_RF), chunksize=1000)

        clean_FR = []
        clean_RF = []

        bad_FR = []
        bad_RF = []

        short_FR = []
        short_RF = []

        count = 0
        count_good = 0
        good_fr_len = 0
        good_rf_len = 0
        count_bad = 0
        count_short = 0
        count_adapters_fr = 0
        count_adapters_rf = 0
        max_adapters_per_read_fr = 0
        max_adapters_per_read_rf = 0
        max_similarity_fr = 0
        max_similarity_rf = 0

        working = True

        while working:
            if not self.run:
                pool.terminate()
                break
            while Gtk.events_pending():
                Gtk.main_iteration()
            try:
                clean_FR_rec, clean_RF_rec, bad_FR_rec, bad_RF_rec, \
                short_FR_rec, short_RF_rec, adapters_count, max_sim = next(results)
                count += 1
                self.progressbar.set_fraction(count / max_count)

                count_adapters_fr += adapters_count[0]
                count_adapters_rf += adapters_count[1]
                self.total_number_of_reads_with_adapters.set_text(
                    '%2.5g / %2.5g' % (count_adapters_fr, count_adapters_rf))
                if adapters_count[0] > max_adapters_per_read_fr:
                    max_adapters_per_read_fr = adapters_count[0]
                if adapters_count[1] > max_adapters_per_read_rf:
                    max_adapters_per_read_rf = adapters_count[1]
                self.max_number_of_adapters.set_text('%d, %d' % (max_adapters_per_read_fr, max_adapters_per_read_rf))
                if clean_FR_rec:
                    clean_FR.append(clean_FR_rec)
                    count_good += 1
                    self.total_number_of_reads_out.set_text('%d (%2.2f %%)' % (count_good, count_good / count * 100))
                    good_fr_len += len(clean_FR_rec.seq)
                if clean_RF_rec:
                    clean_RF.append(clean_RF_rec)
                    good_rf_len += len(clean_RF_rec.seq)
                    good_total_len = good_fr_len + good_rf_len
                    self.total_length_of_reads_out.set_text('%2.5g (%2.2f %%)' % (good_total_len,
                                                                                  good_total_len / self.total_length_of_reads_count * 100))
                if bad_FR_rec:
                    count_bad += 1
                    self.total_number_of_suspicious.set_text('%d (%2.2f %%)' % (count_bad, count_bad / count * 100))
                    bad_FR.append(bad_FR_rec)
                if bad_RF_rec:
                    bad_RF.append(bad_RF_rec)
                if short_FR_rec:
                    count_short += 1
                    self.total_number_of_short.set_text('%d (%2.2f %%)' % (count_short, count_short / count * 100))
                    short_FR.append(short_FR_rec)
                if short_RF_rec:
                    short_RF.append(short_RF_rec)
            except StopIteration:
                working = False
                pass

        return clean_FR, clean_RF, bad_FR, bad_RF, short_FR, short_RF
