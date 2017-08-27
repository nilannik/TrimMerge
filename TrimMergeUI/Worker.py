from __future__ import division, print_function

from TrimMergeUI.utils import clean_records, find_overlap


def count_length(record):
    return len(record.seq)


def compare_reads(PE_records):
    record_FR = PE_records[0]
    record_RF = PE_records[1]
    return record_FR.id == record_RF.id


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
                                                                        merge_worker.correct_pq)
    return ret_FR, ret_RF, overlap, max_match, overlap_FR, overlap_RF, insert_Len, concatenated_seq
