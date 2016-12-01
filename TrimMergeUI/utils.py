from __future__ import division, print_function

import numpy as np


def find_insert_in_record(insert, record, min_length, threshold, verbose=False):
    insert_array = np.array(insert.seq)
    insert_size = insert_array.size
    record_size = len(record.seq)
    start_idx = 0
    positions = []
    similarity = []
    inserts = []
    while start_idx < record_size - min_length:
        stop_idx = start_idx + insert_size
        if stop_idx - record_size > 0:
            offset = stop_idx - record_size
        else:
            offset = 0
        subseq = record.seq[start_idx:stop_idx]
        subseq_array = np.array(list(subseq))
        comparison = subseq_array == insert_array[:insert_size-offset]
        match = np.sum(comparison) / insert_size * 100
        if match >= threshold:
            if verbose: print('Match:', match, '% @', start_idx)
            if verbose: print('====> sequence', subseq)
            if verbose: print('====>  nextera', insert.seq, insert.id)
            inserted = False
            for i in range(len(positions)):
                if abs(positions[i] - start_idx) < insert_size - offset:
                    inserted = True
                    if similarity[i] < match:
                        positions[i] = start_idx
                        similarity[i] = match
                        inserts[i] = subseq
                        break
            if not inserted:
                positions.append(start_idx)
                similarity.append(match)
                inserts.append(subseq)
        start_idx += 1
    return positions, similarity, inserts


def clean_records(record_FR, record_RF, adapters_dict, min_length, similarity, short_read_threshold):
    clean_FR, clean_RF, bad_FR, bad_RF, short_FR, short_RF = None, None, None, None, None, None

    positions_FR_raw, similarity_FR_raw, inserts_FR_raw = find_insert_in_record(adapters_dict['FR'],
                                                                                record_FR,
                                                                                min_length,
                                                                                similarity,
                                                                                verbose=False)

    positions_RF_raw, similarity_RF_raw, inserts_RF_raw = find_insert_in_record(adapters_dict['RF'],
                                                                                record_RF,
                                                                                min_length,
                                                                                similarity,
                                                                                verbose=False)

    num_found = [len(positions_FR_raw), len(positions_FR_raw)]
    try:
        max_similarity = [max(similarity_FR_raw), max(similarity_RF_raw)]
    except ValueError:
        max_similarity = [0, 0]
    positions_FR = []
    positions_RF = []
    if len(positions_FR_raw) > 0:
        m = max(similarity_FR_raw)
        idx_FR = [i for i, j in enumerate(similarity_FR_raw) if j == m]
        positions_FR = [positions_FR_raw[i] for i in idx_FR]
        similarity_FR = [similarity_FR_raw[i] for i in idx_FR]
        inserts_FR = [inserts_FR_raw[i] for i in idx_FR]
    if len(positions_RF_raw) > 0:
        m = max(similarity_RF_raw)
        idx_RF = [i for i, j in enumerate(similarity_RF_raw) if j == m]
        positions_RF = [positions_RF_raw[i] for i in idx_RF]
        similarity_RF = [similarity_RF_raw[i] for i in idx_RF]
        inserts_RF = [inserts_RF_raw[i] for i in idx_RF]

    if len(positions_FR) > 0 and len(positions_FR) == len(positions_RF):
        phred_quality_FR = record_FR.letter_annotations['phred_quality'][:positions_FR[0]]
        phred_quality_RF = record_RF.letter_annotations['phred_quality'][:positions_RF[0]]
        record_FR.letter_annotations = {}
        record_RF.letter_annotations = {}
        record_FR.seq = record_FR.seq[:positions_FR[0]]
        record_RF.seq = record_RF.seq[:positions_RF[0]]
        record_FR.letter_annotations['phred_quality'] = phred_quality_FR
        record_RF.letter_annotations['phred_quality'] = phred_quality_RF
        clean = True
        short = False
    elif len(positions_FR) != len(positions_RF):
        print('\n\n*** Suspicious ***')
        print(positions_FR, positions_RF)
        print('*** Suspicious ***\n\n')
        clean = False
        short = False
    elif len(record_FR.seq) < short_read_threshold:
        clean = False
        short = True
        print('Short Sequence found: %d' % len(record_FR.seq))
        print('FR:', record_FR.seq)
        print('RF:', record_RF.seq)
    else:
        clean = True
        short = False
    if clean:
        clean_FR = record_FR
        clean_RF = record_RF
    elif short:
        short_FR = record_FR
        short_RF = record_RF
    else:
        bad_FR = record_FR
        bad_RF = record_RF
    return clean_FR, clean_RF, bad_FR, bad_RF, short_FR, short_RF, num_found, max_similarity
