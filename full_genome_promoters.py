from Bio import SeqIO, SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import nt_search
import csv
import sys

def search_for_promoters(promoter_list, promoter_type, window_start, window_end, ref_seq, strand):
    unique_search_results = set()
    if strand == 'F':
        seq_to_search = ref_seq[window_start:window_end]
    elif strand == 'R':
        seq_to_search = ref_seq[window_start:window_end].reverse_complement()
    else: print('strand error')
    for promoter in promoter_list:
        search_results = nt_search(str(seq_to_search), promoter)
        num_promoters_found = len(search_results) - 1
        if num_promoters_found > 0:
            for pos in search_results[1:]:
                if strand == 'F':
                    promoter_pos = window_start + pos
                    if promoter_type == -10:
                        promoter_seq = str(ref_seq[(promoter_pos - 3):(promoter_pos + 6)])
                    elif promoter_type == -35:
                        promoter_seq = str(ref_seq[promoter_pos:(promoter_pos + 6)])
                elif strand == 'R':
                    promoter_pos = window_end - pos
                    if promoter_type == -10:
                        promoter_seq = str(ref_seq[(promoter_pos - 6):(promoter_pos + 3)].reverse_complement())
                    elif promoter_type == -35:
                        promoter_seq = str(ref_seq[(promoter_pos - 6):promoter_pos].reverse_complement())
                else: print('strand error')
                unique_search_results.add((promoter_pos, promoter_seq))
    return sorted(unique_search_results)

def record_thirtyfive_ten_promoters(ref_seq, strand, min_len_gap, max_len_gap):
    # Promoter sequences to look for: 
    ### In -10 promoters: A2 and T6 always conserved plus 2 out of 4 variable sites (total of 4/6).
    ### In -35 promoters: T1 and T2 always conserved plus 2 out of 4 variable sites (total of 4/6).
    ten_promoters = ['TATNNT', 'TANANT', 'TANNAT', 'NATANT', 'NATNAT', 'NANAAT']
    thirtyfive_promoters = ['TTGANN', 'TTGNCN', 'TTGNNA', 'TTNACN', 'TTNANA', 'TTNNCA']
    
    full_results = []
    # Since searching the whole genome for -10 promoters, start of window = 0 and end of window = length of genome.
    ten_promoter_results = search_for_promoters(ten_promoters, -10, 0, len(ref_seq), ref_seq, strand)
    if len(ten_promoter_results) > 0:
        for promoter_record in ten_promoter_results:
            ten_promoter_pos = promoter_record[0]
            ten_promoter_seq = promoter_record[1]
            if strand == 'F':
                thirtyfive_win_start = ten_promoter_pos - (max_len_gap + 6)
                thirtyfive_win_end = ten_promoter_pos - min_len_gap
            elif strand == 'R':
                thirtyfive_win_start = ten_promoter_pos + min_len_gap
                thirtyfive_win_end = ten_promoter_pos + (max_len_gap + 6)
            else: print('strand error')
            thirtyfive_promoter_results = search_for_promoters(thirtyfive_promoters, -35, thirtyfive_win_start, thirtyfive_win_end, ref_seq, strand)
            if len(thirtyfive_promoter_results) > 0:
                if ten_promoter_seq.startswith('TG'):
                    full_results.append((ten_promoter_pos, ten_promoter_seq, 'extended', thirtyfive_promoter_results))
                else:
                    full_results.append((ten_promoter_pos, ten_promoter_seq, ' ', thirtyfive_promoter_results))
            else:
                if ten_promoter_seq.startswith('TG'):
                    full_results.append((ten_promoter_pos, ten_promoter_seq, 'extended'))
                else:
                    full_results.append((ten_promoter_pos, ten_promoter_seq, ' '))
    return (full_results)

def find_35_10_promoters(gb_path, outpath):
    ref_record = SeqIO.read(open(gb_path), 'genbank')
    ref_sequence = ref_record.seq
    print('reference length: ', len(ref_sequence))
    # Where to search for -35 promoters relative to -10 promoters
    min_len_gap = 15
    max_len_gap = 20

    with open(outpath, 'w') as outfile:
        writer = csv.writer(outfile, delimiter = ',')
        # Find and record promoters on the forward strand:
        forward_results = record_thirtyfive_ten_promoters(ref_sequence, 'F', min_len_gap, max_len_gap)
        if len(forward_results) > 0:
            for entry in forward_results:
                writer.writerow(['F', entry])
        else: print('no promoters found on forward strand')
        # Find and record promoters on the reverse strand:
        reverse_results = record_thirtyfive_ten_promoters(ref_sequence, 'R', min_len_gap, max_len_gap)
        if len(reverse_results) > 0:
            for entry in reverse_results:
                writer.writerow(['R', entry])
        else: print('no promoters found on forward strand')

if __name__ == '__main__':
    if len(sys.argv) == 3:
         find_35_10_promoters(sys.argv[1], sys.argv[2])
    else:
         print("Usage: full_genome_promoters.py reference_genbank.gb outfile.csv")
         sys.exit(0)