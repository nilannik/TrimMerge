from __future__ import division, print_function

try:
    import itertools.izip as zip
except ImportError:
    pass

from Bio import SeqIO

from TrimMergeUI.utils import clean_records, find_overlap, batch_iterator

from os import pardir
from os.path import join, dirname, basename, realpath, splitext


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
    for i, batch in enumerate(batch_iterator(records, 10000)):
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
    total_length = 0
    count = 0
    for record_FR, record_RF in zip(records_FR, records_RF):
        if not record_FR.id == record_RF.id:
            print('ID mismatch:', records_FR.id, record_RF.id)
            return False
    return True


class TrimWorker(object):

    def __init__(self, adapters_dict, adapter_min_length, adapter_similarity, short_read_threshold):
        super(TrimWorker, self).__init__()
        self.adapters_dict = adapters_dict
        self.adapter_min_length = adapter_min_length
        self.adapter_similarity = adapter_similarity
        self.short_read_threshold = short_read_threshold


class MergeWorker(object):

    def __init__(self, overlap_min_length, overlap_similarity, correct_pq):
        super(MergeWorker, self).__init__()
        self.overlap_min_length = overlap_min_length
        self.overlap_similarity = overlap_similarity
        self.correct_pq = correct_pq


def init_trim_worker(adapters_dict, adapter_min_length, adapter_similarity, short_read_threshold):
    global trim_worker
    trim_worker = TrimWorker(adapters_dict, adapter_min_length, adapter_similarity, short_read_threshold)


def init_merge_worker(overlap_min_length, overlap_similarity, correct_pq):
    global merge_worker
    merge_worker = MergeWorker(overlap_min_length, overlap_similarity, correct_pq)


def clean_reads(PE_records):
    record_FR = PE_records[0]
    record_RF = PE_records[1]
    clean_FR_rec, clean_RF_rec, bad_FR_rec, bad_RF_rec, \
    short_FR_rec, short_RF_rec, adapters_count, max_sim = clean_records(record_FR, record_RF,
                                                                        trim_worker.adapters_dict,
                                                                        trim_worker.adapter_min_length,
                                                                        trim_worker.adapter_similarity,
                                                                        trim_worker.short_read_threshold)
    return clean_FR_rec, clean_RF_rec, bad_FR_rec, bad_RF_rec, short_FR_rec, short_RF_rec, adapters_count, max_sim


def merge_overlaps(PE_records):
    record_FR = PE_records[0]
    record_RF = PE_records[1]
    ret_FR, ret_RF, overlap, max_match, \
    overlap_FR, overlap_RF, insert_Len, concatenated_seq = find_overlap(record_FR, record_RF,
                                                                        merge_worker.overlap_min_length,
                                                                        merge_worker.overlap_similarity,
                                                                        merge_worker.correct_pq,
                                                                        reverse_complement=True)
    return ret_FR, ret_RF, overlap, max_match, overlap_FR, overlap_RF, insert_Len, concatenated_seq


def merge_overlaps_no_reverse_complement(PE_records):
    record_FR = PE_records[0]
    record_RF = PE_records[1]
    ret_FR, ret_RF, overlap, max_match, \
    overlap_FR, overlap_RF, insert_Len, concatenated_seq = find_overlap(record_FR, record_RF,
                                                                        merge_worker.overlap_min_length,
                                                                        merge_worker.overlap_similarity,
                                                                        merge_worker.correct_pq,
                                                                        reverse_complement=False)
    return ret_FR, ret_RF, overlap, max_match, overlap_FR, overlap_RF, insert_Len, concatenated_seq
