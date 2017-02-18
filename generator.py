from __future__ import division, print_function
from os.path import join, dirname, realpath
import numpy as np

from matplotlib import pyplot as plt

from Bio import SeqIO
from Bio.Alphabet import IUPAC

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

alphabet = IUPAC.unambiguous_dna


def read_in_adapters(adapter_name):
    dir_path = join(dirname(realpath(__file__)))
    adapters_dir = join(dir_path, 'adapters')
    print(adapters_dir)
    adapters_dict = dict()
    if adapter_name == 'Nextera PE':
        adapters_file_name = join(adapters_dir, 'NexteraPE-PE_2.txt')
        print(adapters_file_name)
        nextera_records = SeqIO.parse(adapters_file_name, 'fasta', alphabet=IUPAC.ambiguous_dna)
        for record in nextera_records:
            if record.id == 'Trans2_rc':
                adapters_dict['FR'] = record
            elif record.id == 'Trans1_rc':
                adapters_dict['RF'] = record
    return adapters_dict


def generate_seq(length, gc_content=40):
    deck = [['G', 'C'][i] for i in np.random.randint(0, 2, size=int(gc_content))]
    deck += [['A', 'T'][i] for i in np.random.randint(0, 2, size=int(100 - gc_content))]
    np.random.shuffle(deck)
    return ''.join([deck[i] for i in np.random.randint(0, 100, size=int(length))])


def generate_genome(length):
    genome = SeqRecord(Seq(generate_seq(length)), id='Artificial genome', name='Generated genome',
                       description='Genome generated for testing purposes')
    genome.letter_annotations["phred_quality"] = [33] * int(length)
    return genome


def generate_seq_record(genome, adapters, read_length, insert_length, seq_id,
                        max_phred=33, min_phred=18, transition_length=30, change_threshold=25):
    start_position = np.random.randint(len(genome.seq)-insert_length)
    if insert_length > read_length:
        read = str(genome.seq[start_position: start_position+read_length])
        adapter_fr = ''
        noise_fr = ''
        adapter_rf = ''
        noise_rf = ''
    else:
        read = str(genome.seq[start_position: start_position + insert_length])
        bp_left = read_length - insert_length
        adapter_fr_len = len(adapters['FR'].seq)
        adapter_rf_len = len(adapters['RF'].seq)
        if bp_left >= adapter_fr_len:
            adapter_fr = str(adapters['FR'].seq)
            noise_fr = generate_seq(bp_left - adapter_fr_len)
        else:
            adapter_fr = str(adapters['FR'].seq[:read_length - insert_length])
            noise_fr = ''
        if bp_left >= adapter_rf_len:
            adapter_rf = str(adapters['RF'].seq)
            noise_rf = generate_seq(bp_left - adapter_rf_len)
        else:
            adapter_rf = str(adapters['RF'].seq[:read_length - insert_length])
            noise_rf = ''
    seq_fr_array = np.array(list(read + adapter_fr + noise_fr))
    _, pq = generate_phred_quality(read_length, max_phred, min_phred, transition_length)
    change_idx = np.where(pq < change_threshold)
    change_array = [alphabet.letters[i] for i in np.random.randint(0, 3, len(change_idx))]
    seq_fr_array[change_idx] = change_array
    seq_fr = Seq(''.join(seq_fr_array))
    seq_fr = SeqRecord(seq=seq_fr, id=seq_id, name='Generated sequence',
                       description='Generated read for algorithms testing')
    seq_fr.letter_annotations["phred_quality"] = pq

    seq_rf = Seq(noise_rf + adapter_rf + str(Seq(read).reverse_complement()))
    seq_rf_array = np.array(list(str(seq_rf)))
    _, pq = generate_phred_quality(read_length, max_phred, min_phred, transition_length)
    change_idx = np.where(pq < change_threshold)
    change_array = [alphabet.letters[i] for i in np.random.randint(0, 3, len(change_idx))]
    seq_rf_array[change_idx] = change_array
    seq_rf = Seq(''.join(seq_rf_array))
    seq_rf = SeqRecord(seq=seq_rf, id=seq_id, name='Generated sequence',
                       description='Generated read for algorithms testing')
    seq_rf.letter_annotations["phred_quality"] = pq

    insert = SeqRecord(Seq(read), id=seq_id, name='Generated sequence',
                       description='Clean insert of generated read')
    insert.letter_annotations["phred_quality"] = [33] * len(read)
    return seq_fr, seq_rf, insert


def generate_phred_quality(read_length, max_phred, min_phred, transition_length):
    x = np.arange(read_length)
    d = 1 - (max_phred - min_phred)/max_phred * np.exp((x - read_length) / transition_length)
    pq = max_phred - np.random.randint(min_phred/2, 2*max_phred, size=read_length) * (1 - d)
    return d * max_phred, pq


def plot_seq_diagram(read_length, max_phred, min_phred, transition_length):
    x = np.arange(read_length)
    change_threshold = 25
    change_threshold_array = np.ones_like(x) * change_threshold
    for i in range(1):
        d, pq = generate_phred_quality(read_length, max_phred, min_phred, transition_length)
        plt.plot(x, pq, 'bo')
    plt.plot(x, d, '-r')
    plt.plot(x, change_threshold_array, 'k--')
    change_idx = np.where(pq < change_threshold_array)
    markers = np.zeros_like(x)
    markers[change_idx] = change_threshold
    plt.plot(x, markers, 'go')
    plt.show()


if __name__ == '__main__':
    dir_path = join(dirname(realpath(__file__)))
    data_dir = join(dir_path, 'data')
    print(data_dir)
    fr_file_name = join(data_dir, 'reads_FR.fastq')
    rf_file_name = join(data_dir, 'reads_RF.fastq')
    inserts_file_name = join(data_dir, 'inserts.fastq')
    genome_file_name = join(data_dir, 'genome_contig.fastq')

    ref_genome = generate_genome(1e6)
    adapters = read_in_adapters('Nextera PE')
    reads_fr = []
    reads_rf = []
    clean_reads = []
    min_insert_length = 10
    max_insert_length = 300
    center_insert_length = 200
    deviation_insert_length = 30
    read_length = 250
    lib_size = 65*1e3
    insert_length = np.random.normal(center_insert_length, deviation_insert_length, size=lib_size)
    insert_length[np.where(insert_length < min_insert_length)] = min_insert_length
    insert_length[np.where(insert_length > max_insert_length)] = max_insert_length
    for seq_id in range(int(lib_size)):
        my_seq_fr, my_seq_rf, insert = generate_seq_record(ref_genome, adapters,
                                                           read_length, np.int(insert_length[seq_id]),
                                                           seq_id='TEST_READ_%04d' % (seq_id+1))
        reads_fr.append(my_seq_fr)
        reads_rf.append(my_seq_rf)
        clean_reads.append(insert)
    SeqIO.write(sequences=reads_fr, handle=fr_file_name, format='fastq')
    SeqIO.write(sequences=reads_rf, handle=rf_file_name, format='fastq')
    SeqIO.write(sequences=clean_reads, handle=inserts_file_name, format='fastq')
    SeqIO.write(sequences=ref_genome, handle=genome_file_name, format='fastq')
    #plot_seq_diagram(read_length=255, max_phred=33, min_phred=18, transition_length=50)
    plt.hist(insert_length)
    plt.show()
