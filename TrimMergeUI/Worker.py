from __future__ import division, print_function

from TrimMergeUI.utils import clean_records


def count_length(record):
    return len(record.seq)


def compare_reads(PE_records):
    record_FR = PE_records[0]
    record_RF = PE_records[1]
    return record_FR.id == record_RF.id


class TrimmerWorker(object):

    def __init__(self, adapters_dict, adapter_min_length, adapter_similarity, short_read_threshold):
        super(TrimmerWorker, self).__init__()
        self.adapters_dict = adapters_dict
        self.adapter_min_length = adapter_min_length
        self.adapter_similarity = adapter_similarity
        self.short_read_threshold = short_read_threshold


def init_worker(adapters_dict, adapter_min_length, adapter_similarity, short_read_threshold):
    global worker
    worker = TrimmerWorker(adapters_dict, adapter_min_length, adapter_similarity, short_read_threshold)


def clean_reads(PE_records):
    record_FR = PE_records[0]
    record_RF = PE_records[1]
    clean_FR_rec, clean_RF_rec, bad_FR_rec, bad_RF_rec, \
    short_FR_rec, short_RF_rec, adapters_count, max_sim = clean_records(record_FR, record_RF,
                                                                        worker.adapters_dict,
                                                                        worker.adapter_min_length,
                                                                        worker.adapter_similarity,
                                                                        worker.short_read_threshold)
    return clean_FR_rec, clean_RF_rec, bad_FR_rec, bad_RF_rec, short_FR_rec, short_RF_rec, adapters_count, max_sim
