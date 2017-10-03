from __future__ import division, print_function

import numpy as np


def batch_iterator(iterator, batch_size):
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        total_len = 0
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
            total_len += len(entry.seq)
        if batch:
            yield batch, total_len


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


def clean_records(records_FR, records_RF, adapters_dict, min_length, similarity, short_read_threshold):
    clean_FR = []
    count_clean = 0
    clean_fr_len = 0
    clean_RF = []
    clean_rf_len = 0
    clean_total_len = 0

    bad_FR = []
    count_bad = 0
    bad_RF = []

    short_FR = []
    count_short = 0
    short_RF = []

    num_found = [0, 0]
    max_similarity = [0, 0]

    count = 0

    for record_FR, record_RF in list(zip(records_FR, records_RF)):
        clean_FR_i, clean_RF_i, bad_FR_i, bad_RF_i,\
        short_FR_i, short_RF_i,\
        num_found_i, max_similarity_i = clean_record(record_FR, record_RF, adapters_dict, min_length,
                                                     similarity, short_read_threshold)
        count += 1
        num_found[0] += num_found_i[0]
        num_found[1] += num_found_i[1]
        max_similarity[0] += max_similarity_i[0]
        max_similarity[1] += max_similarity_i[1]

        if clean_FR_i:
            clean_FR.append(clean_FR_i)
            count_clean += 1
            clean_fr_len += len(clean_FR_i.seq)
        if clean_RF_i:
            clean_RF.append(clean_RF_i)
            clean_rf_len += len(clean_RF_i.seq)
        clean_total_len = clean_fr_len + clean_rf_len
        if bad_FR_i:
            count_bad += 1
            bad_FR.append(bad_FR_i)
        if bad_RF_i:
            bad_RF.append(bad_RF_i)
        if short_FR_i:
            count_short += 1
            short_FR.append(short_FR_i)
        if short_RF_i:
            short_RF.append(short_RF_i)

    return clean_FR, clean_RF, bad_FR, bad_RF, short_FR, short_RF, num_found, max_similarity,\
        count, count_clean, count_bad, count_short, clean_fr_len, clean_rf_len, clean_total_len


def clean_record(record_FR, record_RF, adapters_dict, min_length, similarity, short_read_threshold):
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
        #print('\n\n*** Suspicious ***')
        #print(positions_FR, positions_RF)
        #print('*** Suspicious ***\n\n')
        clean = False
        short = False
    elif len(record_FR.seq) < short_read_threshold:
        clean = False
        short = True
        #print('Short Sequence found: %d' % len(record_FR.seq))
        #print('FR:', record_FR.seq)
        #print('RF:', record_RF.seq)
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


def find_overlaps(records_FR, records_RF, min_length, threshold, correct_pq=True, reverse_complement=True):
    concat_FR = []
    bad_seq_FR = []
    bad_seq_RF = []
    overlap_len = []
    matches = []
    insert_len = []
    count = 0
    count_overlapping = 0
    count_not_overlapping = 0
    for record_FR, record_RF in list(zip(records_FR, records_RF)):
        ret_FR, ret_RF, overlap, max_match,\
        overlap_FR, overlap_RF, insert_Len, concatenated_seq = find_overlap(record_FR, record_RF,
                                                                            min_length, threshold,
                                                                            correct_pq,
                                                                            reverse_complement)
        if max_match >= threshold:
            concat_FR.append(concatenated_seq)
            insert_len.append(insert_Len)
            overlap_len.append(len(overlap.seq))
            matches.append(max_match)
            count_overlapping += 1
        else:
            bad_seq_FR.append(ret_FR)
            bad_seq_RF.append(ret_RF)
            count_not_overlapping += 1
        count += 1
    return concat_FR, bad_seq_FR, bad_seq_RF, overlap_len, matches, insert_len,\
           count, count_overlapping, count_not_overlapping


def find_overlap(seq_FR, seq_RF, min_length, threshold, correct_pq=True, reverse_complement=True):
    str_FR = np.array(list(seq_FR.seq))
    if reverse_complement:
        str_RF = np.array(list(seq_RF.reverse_complement().seq))
    else:
        str_RF = np.array(list(seq_RF.seq))
    insert_Len = str_FR.size + str_RF.size
    L_min = min(str_FR.size, str_RF.size)
    L_max = max(str_FR.size, str_RF.size)
    if str_FR.size >= str_RF.size:
        long_str = str_FR
        short_str = str_RF
        long_seq = seq_FR
        if reverse_complement:
            short_seq = seq_RF.reverse_complement()
        else:
            short_seq = seq_RF
        short_seq.id = seq_RF.id
    else:
        long_str = str_RF
        short_str = str_FR
        if reverse_complement:
            long_seq = seq_RF.reverse_complement()
        else:
            long_seq = seq_RF
        long_seq.id = seq_RF.id
        short_seq = seq_FR
    if reverse_complement:
        concatenated_seq = seq_FR + seq_RF.reverse_complement()
    else:
        concatenated_seq = seq_FR + seq_RF
    concatenated_seq.id = seq_FR.id
    start_long = 0
    start_short = max(0, L_min - min_length)
    L = min(L_min, min_length)
    max_match = 0
    max_match_start_long = start_long
    max_match_start_short = start_short
    max_match_L = L
    while start_long <= L_max - min_length:
        comparison = long_str[start_long:start_long + L] == short_str[start_short:start_short + L]
        match = np.sum(comparison) / L * 100
        if match > max_match:
            max_match = match
            max_match_start_long = start_long
            max_match_start_short = start_short
            max_match_L = L
        if start_short > 0:
            start_short -= 1
            L += 1
        else:
            start_long += 1
            L = min(L_max - start_long, L)
    if max_match >= threshold:
        # print(max_match_start_short, max_match_start_long, max_match_L, max_match)
        insert_Len = str_FR.size + str_RF.size - max_match_L
        long_substr = long_str[max_match_start_long:max_match_L + max_match_start_long]
        short_substr = short_str[max_match_start_short:max_match_start_short + max_match_L]
        comparison = long_substr == short_substr
        faults = np.where(comparison == 0)[0]
        if correct_pq:
            for idx in range(max_match_L):
                PQ_long = np.array(long_seq.letter_annotations['phred_quality'])[max_match_start_long + idx]
                PQ_short = np.array(short_seq.letter_annotations['phred_quality'])[max_match_start_short + idx]
                if long_str[max_match_start_long + idx] != short_str[max_match_start_short + idx]:
                    if PQ_long > PQ_short:
                        id_backup = short_seq.id
                        tmp_seq = short_seq.seq.tomutable()
                        tmp_seq[max_match_start_short + idx] = long_str[max_match_start_long + idx]
                        annot_backup = short_seq.letter_annotations
                        short_seq.letter_annotations = {}
                        short_seq.seq = tmp_seq.toseq()
                        short_seq.letter_annotations = annot_backup
                        short_seq.id = id_backup
                    elif PQ_long < PQ_short:
                        id_backup = long_seq.id
                        tmp_seq = long_seq.seq.tomutable()
                        tmp_seq[max_match_start_long + idx] = short_str[max_match_start_short + idx]
                        annot_backup = long_seq.letter_annotations
                        long_seq.letter_annotations = {}
                        long_seq.seq = tmp_seq.toseq()
                        long_seq.letter_annotations = annot_backup
                        long_seq.id = id_backup
                if PQ_long != PQ_short:
                    long_seq.letter_annotations['phred_quality'][max_match_start_long + idx] = max(PQ_long, PQ_short)
                    short_seq.letter_annotations['phred_quality'][max_match_start_short + idx] = max(PQ_long, PQ_short)

    overlap = long_seq[max_match_start_long:max_match_start_long + max_match_L]

    if str_FR.size >= str_RF.size:
        if reverse_complement:
            if max_match_start_short > 0:
                concatenated_seq = seq_RF.reverse_complement()[0:max_match_start_short] + seq_FR[:]
            else:
                concatenated_seq = seq_FR[0:max_match_start_long] + seq_RF.reverse_complement()[:]
        else:
            if max_match_start_short > 0:
                concatenated_seq = seq_RF[0:max_match_start_short] + seq_FR[:]
            else:
                concatenated_seq = seq_FR[0:max_match_start_long] + seq_RF[:]
        concatenated_seq.id = seq_FR.id
        overlap_FR = seq_FR[max_match_start_long:max_match_start_long + max_match_L]
        overlap_FR.id = seq_FR.id
        if reverse_complement:
            overlap_RF = seq_RF.reverse_complement()[max_match_start_short:max_match_start_short + max_match_L]\
                .reverse_complement()
        else:
            overlap_RF = seq_RF[max_match_start_short:max_match_start_short + max_match_L]
        overlap_RF.id = seq_RF.id
        ret_FR = long_seq
        if reverse_complement:
            ret_RF = short_seq.reverse_complement()
            ret_RF.id = seq_RF.id
        else:
            ret_RF = short_seq
    else:
        if reverse_complement:
            if max_match_start_short > 0:
                concatenated_seq = seq_FR[0:max_match_start_short] + seq_RF.reverse_complement()[:]
            else:
                concatenated_seq = seq_RF.reverse_complement()[0:max_match_start_long] + seq_FR[:]
        else:
            if max_match_start_short > 0:
                concatenated_seq = seq_FR[0:max_match_start_short] + seq_RF[:]
            else:
                concatenated_seq = seq_RF[0:max_match_start_long] + seq_FR[:]
        concatenated_seq.id = seq_FR.id
        if reverse_complement:
            overlap_RF = seq_RF.reverse_complement()[max_match_start_long:max_match_start_long + max_match_L]\
                .reverse_complement()
        else:
            overlap_RF = seq_RF[max_match_start_long:max_match_start_long + max_match_L]
        overlap_RF.id = seq_RF.id
        overlap_FR = seq_FR[max_match_start_short:max_match_start_short + max_match_L]
        overlap_FR.id = seq_FR.id
        ret_FR = short_seq
        if reverse_complement:
            ret_RF = long_seq.reverse_complement()
            ret_RF.id = seq_RF.id
        else:
            ret_RF = long_seq
    if L_max < min_length:
        insert_Len = 0
    return ret_FR, ret_RF, overlap, max_match, overlap_FR, overlap_RF, insert_Len, concatenated_seq
