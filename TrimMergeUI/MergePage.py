from __future__ import division, print_function

from os import remove
from os.path import join, basename, splitext

try:
    import itertools.izip as zip
except ImportError:
    pass

from Bio import SeqIO
from Bio.Alphabet import IUPAC

import gi

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas

from multiprocessing import Pool

from TrimMergeUI.Worker import partition_file, count_length, compare_reads, init_merge_worker, merge_overlaps
from TrimMergeUI.utils import batch_iterator


class MergePage(Gtk.Box):

    def __init__(self):
        super(MergePage, self).__init__()

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
        statistics_box.pack_start(input_stat_grid, False, False, 0)
        self.total_number_of_reads = Gtk.Label('n/a')
        self.total_length_of_reads = Gtk.Label('n/a')
        input_stat_grid.attach(Gtk.Label('Reads count'), 0, 0, 1, 1)
        input_stat_grid.attach(self.total_number_of_reads, 1, 0, 1, 1)
        input_stat_grid.attach(Gtk.Label('Total length'), 0, 1, 1, 1)
        input_stat_grid.attach(self.total_length_of_reads, 1, 1, 1, 1)
        stat_title = Gtk.Label()
        stat_title.set_markup('_'*50 + '\n<big><b>Contigs statistics</b></big>')
        statistics_box.pack_start(stat_title, False, False, 0)
        output_stat_grid = Gtk.Grid(column_spacing=10, row_spacing=10)
        statistics_box.pack_start(output_stat_grid, False, False, 0)
        self.total_number_of_reads_out = Gtk.Label('n/a')
        self.total_number_of_contigs_out = Gtk.Label('n/a')
        self.total_number_of_bad = Gtk.Label('n/a')
        output_stat_grid.attach(Gtk.Label('Reads count'), 0, 0, 1, 1)
        output_stat_grid.attach(self.total_number_of_reads_out, 1, 0, 1, 1)
        output_stat_grid.attach(Gtk.Label('Contigs built'), 0, 1, 1, 1)
        output_stat_grid.attach(self.total_number_of_contigs_out, 1, 1, 1, 1)
        output_stat_grid.attach(Gtk.Label('Not overlapping'), 0, 2, 1, 1)
        output_stat_grid.attach(self.total_number_of_bad, 1, 2, 1, 1)

        self.fig = Figure()
        self.ax_overlaps = self.fig.add_subplot(211)
        self.ax_overlaps.set_title('overlap len')
        self.ax_insert_len = self.fig.add_subplot(212)
        self.ax_insert_len.set_title('insert len')
        self.fig.subplots_adjust(left=0.2, hspace=0.3, top=0.95, bottom=0.1, right=0.9, wspace=0.0)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.set_size_request(250, 380)
        hbox.pack_start(self.canvas, True, True, 0)

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

        self.reverse_complement_rf = Gtk.CheckButton("Reverse complement RF file")
        self.reverse_complement_rf.set_active(True)

        self.select_out_dir_button = Gtk.Button("Choose output directory")
        self.select_out_dir_button.connect("clicked", self.on_folder_clicked)
        self.output_dir_name = None
        self.output_dir_label = Gtk.Label('Please choose output dir')

        self.overlap_min_length = 20
        adjustment = Gtk.Adjustment(self.overlap_min_length, 3, 100, 1, 10, 0)
        self.overlap_min_length_button = Gtk.SpinButton()
        self.overlap_min_length_button.set_adjustment(adjustment)
        self.overlap_min_length_button.set_numeric(True)
        self.overlap_min_length_button.set_update_policy(Gtk.SpinButtonUpdatePolicy.IF_VALID)
        self.overlap_min_length_button.connect("value-changed", self.on_overlap_min_length_changed)
        self.overlap_min_length_button.set_value(self.overlap_min_length)

        self.overlap_similarity = 60
        adjustment = Gtk.Adjustment(self.overlap_similarity, 10, 100, 1, 10, 0)
        self.overlap_similarity_button = Gtk.SpinButton()
        self.overlap_similarity_button.set_adjustment(adjustment)
        self.overlap_similarity_button.set_numeric(True)
        self.overlap_similarity_button.set_update_policy(Gtk.SpinButtonUpdatePolicy.IF_VALID)
        self.overlap_similarity_button.connect("value-changed", self.on_overlap_similarity_changed)
        self.overlap_similarity_button.set_value(self.overlap_similarity)

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

        grid.attach(self.reverse_complement_rf, 0, 4, 2, 1)

        grid.attach(self.select_out_dir_button, 0, 5, 1, 1)
        grid.attach(self.output_dir_label, 1, 5, 1, 1)
        grid.attach(Gtk.Label('Min overlap length:'), 0, 6, 1, 1)
        grid.attach(self.overlap_min_length_button, 1, 6, 1, 1)
        grid.attach(Gtk.Label('Overlap similarity %:'), 0, 7, 1, 1)
        grid.attach(self.overlap_similarity_button, 1, 7, 1, 1)
        grid.attach(self.run_button, 0, 8, 2, 1)
        grid.attach(self.progressbar, 0, 9, 2, 1)
        grid.attach(self.status_label, 0, 10, 2, 1)
        self.show_all()
        self.lock_label.hide()

        self.temp_files_FR = []
        self.temp_files_RF = []

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

    def on_overlap_min_length_changed(self, widget):
        value = int(widget.get_value())
        print("Overlap min length changed=%d" % value)
        self.overlap_min_length = value

    def on_overlap_similarity_changed(self, widget):
        value = int(widget.get_value())
        print("Overlap similarity changed=%d" % value)
        self.overlap_similarity = value

    def on_run_toggled(self, button):
        if button.get_active():
            state = "on"
            if self.output_dir_name is None:
                print('Select output dir first!')
                self.status_label.set_text('Select output dir first!')
                self.run_button.set_active(False)
                return
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

    def _file_name_label_update(self, label, file_name):
        if len(file_name) > 24:
            short_name = '...' + file_name[-21:]
        else:
            short_name = file_name
        label.set_text(short_name)

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
        self.reverse_complement_rf.set_sensitive(release)
        self.select_out_dir_button.set_sensitive(release)
        self.overlap_min_length_button.set_sensitive(release)
        self.overlap_similarity_button.set_sensitive(release)
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

        print('Partitioning files')
        response, max_count = self.partition_files()
        if not response:
            self.run_button.set_active(False)
            return

        print('Comparing FR and RF records')
        response = self.compare_records()
        if not response:
            self.run_button.set_active(False)
            return

        print('Merging')
        concat_FR, bad_seq_FR, bad_seq_RF, overlap_len, insert_len, matches = self.merge_PE_reads(max_count)

        fr_base_name, fr_extension = splitext(basename(self.fr_file_name))
        rf_base_name, rf_extension = splitext(basename(self.rf_file_name))

        file_out_FR_concat = join(self.output_dir_name, fr_base_name + '_contigs' + fr_extension)
        file_out_FR_bad = join(self.output_dir_name, fr_base_name + '_no_overlap' + fr_extension)
        file_out_RF_bad = join(self.output_dir_name, rf_base_name + '_no_overlap' + fr_extension)
        insert_len_file = join(self.output_dir_name, fr_base_name + '_insert_len.txt')
        overlap_len_file = join(self.output_dir_name, fr_base_name + '_overlap_len.txt')
        matches_file = join(self.output_dir_name, fr_base_name + '_matches.txt')
        f = open(insert_len_file, 'w')
        for item in insert_len:
            f.write("%s\n" % item)
        f.close()
        f = open(overlap_len_file, 'w')
        for item in overlap_len:
            f.write("%s\n" % item)
        f.close()
        f = open(matches_file, 'w')
        for item in matches:
            f.write("%s\n" % item)
        f.close()

        SeqIO.write(concat_FR, file_out_FR_concat, self.file_format)
        SeqIO.write(bad_seq_FR, file_out_FR_bad, self.file_format)
        SeqIO.write(bad_seq_RF, file_out_RF_bad, self.file_format)

        print('Removing temp files')
        self.delete_temp_files()

        self.run_button.set_active(False)
        self.status_label.set_text('Idle')

    def delete_temp_files(self):
        for file_descriptor in self.temp_files_FR + self.temp_files_RF:
            remove(file_descriptor['file_name'])
        self.temp_files_FR = []
        self.temp_files_RF = []

    def partition_files(self):
        self.status_label.set_text('Partitioning files...')
        chunk_size = 250
        records_FR = SeqIO.parse(self.fr_file_name, self.file_format, alphabet=self.alphabet)
        records_RF = SeqIO.parse(self.rf_file_name, self.file_format, alphabet=self.alphabet)
        records_FR_iterator = batch_iterator(records_FR, chunk_size)
        records_RF_iterator = batch_iterator(records_RF, chunk_size)
        prefix_FR, extension_FR = splitext(basename(self.fr_file_name))
        prefix_RF, extension_RF = splitext(basename(self.rf_file_name))
        count_FR_files = 0
        count_RF_files = 0
        count_FR = 0
        count_RF = 0
        len_fr = 0
        len_rf = 0
        fr_done = False
        rf_done = False
        while not (fr_done and rf_done):
            while Gtk.events_pending():
                Gtk.main_iteration()
            try:
                batch, len_fr_batch = next(records_FR_iterator)
                file_name = join(self.output_dir_name, prefix_FR + "_%04i" % (count_FR_files + 1) + extension_FR)
                self.temp_files_FR.append({
                    'file_name': file_name,
                    'format': self.file_format,
                    'alphabet': self.alphabet,
                    'out_dir': self.output_dir_name
                })
                with open(file_name, "w") as handle:
                    count = SeqIO.write(batch, handle, self.file_format)
                print("Wrote %i records to %s" % (count, file_name))
                count_FR_files += 1
                count_FR += len(batch)
                len_fr += len_fr_batch
            except StopIteration:
                fr_done = True
                pass
            try:
                batch, len_rf_batch = next(records_RF_iterator)
                file_name = join(self.output_dir_name, prefix_RF + "_%04i" % (count_RF_files + 1) + extension_RF)
                self.temp_files_RF.append({
                    'file_name': file_name,
                    'format': self.file_format,
                    'alphabet': self.alphabet,
                    'out_dir': self.output_dir_name
                })
                with open(file_name, "w") as handle:
                    count = SeqIO.write(batch, handle, self.file_format)
                print("Wrote %i records to %s" % (count, file_name))
                count_RF_files += 1
                count_RF += len(batch)
                len_rf += len_rf_batch
            except StopIteration:
                rf_done = True
                pass
            self.total_number_of_reads.set_text('%d' % count_FR)
            self.total_length_of_reads_count = len_fr + len_rf
            self.total_length_of_reads.set_text('%g' % self.total_length_of_reads_count)
        print('Split FR and RF into %d and %d parts' % (len(self.temp_files_FR), len(self.temp_files_RF)))
        if count_FR != count_RF:
            dialog = Gtk.MessageDialog(self.get_toplevel(), 0, Gtk.MessageType.ERROR,
                                       Gtk.ButtonsType.CANCEL, "Different length of FR and RF files")
            dialog.format_secondary_text("Number of forward and reverse reads must match")
            dialog.run()
            dialog.destroy()
            result = False, 0
        else:
            result = True, count_FR
        return result

    def compare_records(self):
        self.status_label.set_text('Checking sequences...')

        print('Creating pool')
        pool = Pool()
        print('Run FR-RF comparison')
        cmp_iter = pool.imap(compare_reads, list(zip(self.temp_files_FR, self.temp_files_RF)))
        cmp_done = False
        result = True
        while not cmp_done:
            while Gtk.events_pending():
                Gtk.main_iteration()
            try:
                result = next(cmp_iter)
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
        print('Closing pool')
        pool.close()
        pool.join()
        self.status_label = Gtk.Label('Idle')
        return result

    def merge_PE_reads(self, max_count=1e10):
        self.status_label.set_text('Scanning for overlaps...')
        print('Creating pool')
        pool = Pool(initializer=init_merge_worker, initargs=(self.overlap_min_length,
                                                             self.overlap_similarity,
                                                             True,
                                                             self.reverse_complement_rf.get_active()))

        results = pool.imap(merge_overlaps, list(zip(self.temp_files_FR, self.temp_files_RF)))

        concat_FR = []
        bad_seq_FR = []
        bad_seq_RF = []
        overlap_len = []
        matches = []
        insert_len = []

        count = 0
        count_overlapping = 0
        count_not_overlapping = 0
        last_plot_update = 0
        plot_update_interval = 1000

        working = True

        while working:
            if not self.run:
                pool.terminate()
                break
            while Gtk.events_pending():
                Gtk.main_iteration()
            try:
                concat_FR_i, bad_seq_FR_i, bad_seq_RF_i, overlap_len_i, matches_i, insert_len_i,\
                count_i, count_overlapping_i, count_not_overlapping_i = next(results)
                count += count_i
                count_overlapping += count_overlapping_i
                count_not_overlapping += count_not_overlapping_i
                concat_FR += concat_FR_i
                bad_seq_FR += bad_seq_FR_i
                bad_seq_RF += bad_seq_RF_i
                overlap_len += overlap_len_i
                matches += matches_i
                insert_len += insert_len_i
                self.progressbar.set_fraction(count / max_count)
                self.total_number_of_reads_out.set_text('%d (%2.2f %%)' % (count, count / max_count * 100))

                self.total_number_of_contigs_out.set_text('%d (%2.2f %%)' % (count_overlapping,
                                                                             count_overlapping / count * 100))
                if count - last_plot_update > plot_update_interval:
                    last_plot_update = count
                    self.ax_overlaps.cla()
                    [n, X, V] = self.ax_overlaps.hist(overlap_len, bins=50)#, normed=True)
                    self.ax_overlaps.set_title('overlap len')
                    self.ax_insert_len.cla()
                    [n, X, V] = self.ax_insert_len.hist(insert_len, bins=50)  # , normed=True)
                    self.ax_insert_len.set_title('insert len')
                    self.fig.canvas.draw()

                self.total_number_of_bad.set_text('%d (%2.2f %%)' % (count_not_overlapping,
                                                                     count_not_overlapping / count * 100))
            except StopIteration:
                working = False
                pass
        print('Closing pool')
        pool.close()
        pool.join()
        self.status_label = Gtk.Label('Idle')
        return concat_FR, bad_seq_FR, bad_seq_RF, overlap_len, insert_len, matches
