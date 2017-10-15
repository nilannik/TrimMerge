from __future__ import division, print_function

try:
    import itertools.izip as zip
except ImportError:
    pass

from Bio import SeqIO

from TrimMergeUI.utils import clean_records, find_overlaps, batch_iterator

from os.path import join, basename, splitext


def count_length(parameters):
    records = SeqIO.parse(parameters['file_name'], parameters['format'], alphabet=parameters['alphabet'])
    total_length = 0
    count = 0
    for record in records:
        total_length += len(record.seq)
        count += 1
    return total_length, count


def partition_file(parameters):
    records = SeqIO.parse(parameters['file_name'], parameters['format'], alphabet=parameters['alphabet'])
    out_dir = parameters['out_dir']
    prefix, extension = splitext(basename(parameters['file_name']))
    output = []
    for i, batch in enumerate(batch_iterator(records, 100)):
        file_name = join(out_dir, prefix + "_%04i" % (i + 1) + extension)
        output.append({
            'file_name': file_name,
            'format': parameters['format'],
            'alphabet': parameters['alphabet'],
            'out_dir': parameters['out_dir']
        })
        with open(file_name, "w") as handle:
            count = SeqIO.write(batch, handle, parameters['format'])
        print("Wrote %i records to %s" % (count, file_name))
    return output


def compare_reads(parameters):
    parameters_FR = parameters[0]
    parameters_RF = parameters[1]
    records_FR = SeqIO.parse(parameters_FR['file_name'], parameters_FR['format'], alphabet=parameters_FR['alphabet'])
    records_RF = SeqIO.parse(parameters_RF['file_name'], parameters_RF['format'], alphabet=parameters_RF['alphabet'])
    for record_FR, record_RF in list(zip(records_FR, records_RF)):
        if not record_FR.id == record_RF.id:
            print('ID mismatch:', records_FR.id, record_RF.id)
            return False
        try:
            if not record_FR.name == record_RF.name:
                print('Name mismatch:', records_FR.name, record_RF.name)
                return False
        except AttributeError:
            pass
        try:
            if not record_FR.description == record_RF.description:
                print('Description mismatch:', records_FR.description, record_RF.description)
                return False
        except AttributeError:
            pass
    return True


class TrimWorker(object):

    def __init__(self, adapters_dict, adapter_min_length, adapter_similarity, short_read_threshold):
        super(TrimWorker, self).__init__()
        self.adapters_dict = adapters_dict
        self.adapter_min_length = adapter_min_length
        self.adapter_similarity = adapter_similarity
        self.short_read_threshold = short_read_threshold


class MergeWorker(object):

    def __init__(self, overlap_min_length, overlap_similarity, correct_pq, reverse_complement):
        super(MergeWorker, self).__init__()
        self.overlap_min_length = overlap_min_length
        self.overlap_similarity = overlap_similarity
        self.correct_pq = correct_pq
        self.reverse_complement = reverse_complement


def init_trim_worker(adapters_dict, adapter_min_length, adapter_similarity, short_read_threshold):
    global trim_worker
    trim_worker = TrimWorker(adapters_dict, adapter_min_length, adapter_similarity, short_read_threshold)


def init_merge_worker(overlap_min_length, overlap_similarity, correct_pq, reverse_complement):
    global merge_worker
    merge_worker = MergeWorker(overlap_min_length, overlap_similarity, correct_pq, reverse_complement)


def clean_reads(parameters):
    parameters_FR = parameters[0]
    parameters_RF = parameters[1]
    records_FR = list(SeqIO.parse(parameters_FR['file_name'], parameters_FR['format'], alphabet=parameters_FR['alphabet']))
    records_RF = list(SeqIO.parse(parameters_RF['file_name'], parameters_RF['format'], alphabet=parameters_RF['alphabet']))
    clean_FR, clean_RF, bad_FR, bad_RF, short_FR, short_RF, num_found, max_similarity, \
    count, count_clean, count_bad, count_short,\
    clean_fr_len, clean_rf_len, clean_total_len = clean_records(list(records_FR), list(records_RF),
                                                                trim_worker.adapters_dict,
                                                                trim_worker.adapter_min_length,
                                                                trim_worker.adapter_similarity,
                                                                trim_worker.short_read_threshold)
    return clean_FR, clean_RF, bad_FR, bad_RF, short_FR, short_RF, num_found, max_similarity, \
        count, count_clean, count_bad, count_short,\
        clean_fr_len, clean_rf_len, clean_total_len


def merge_overlaps(parameters):
    parameters_FR = parameters[0]
    parameters_RF = parameters[1]
    records_FR = SeqIO.parse(parameters_FR['file_name'], parameters_FR['format'], alphabet=parameters_FR['alphabet'])
    records_RF = SeqIO.parse(parameters_RF['file_name'], parameters_RF['format'], alphabet=parameters_RF['alphabet'])
    concat_FR, bad_seq_FR, bad_seq_RF,\
    overlap_len, matches, insert_len,\
    count, count_overlapping, count_not_overlapping = find_overlaps(list(records_FR), list(records_RF),
                                                                    merge_worker.overlap_min_length,
                                                                    merge_worker.overlap_similarity,
                                                                    merge_worker.correct_pq,
                                                                    merge_worker.reverse_complement)
    return concat_FR, bad_seq_FR, bad_seq_RF, overlap_len, matches, insert_len,\
           count, count_overlapping, count_not_overlapping
