from functools import cmp_to_key
import re
import sys, threading
import math
import timeit


class BWT:

    def __init__(self):
        pass

    def comparator(self, index1, index2):
        seq_len = len(self.sequence)

        while (True):
            char1 = self.sequence[index1]
            char2 = self.sequence[index2]
            if (char2 < char1):
                return 1
            if (char1 < char2):
                return -1
            index1 = (index1 + 1) % seq_len
            index2 = (index2 + 1) % seq_len
        return 0

    def rotation_with_indices(self, sequence):

        """ Given a pattern/sequence, returns a vector of sequences.
            Other sequences are only rotations of the input sequence. """

        seq_len = len(sequence)
        index_vector = [""] * seq_len
        for i in range(len(index_vector)):
            index_vector[i] = i

        rot_vector = [""] * seq_len
        self.sequence = sequence
        self.vector = index_vector
        index_vector.sort(key=cmp_to_key(self.comparator))
        return index_vector

    def rotation(self, sequence):

        """ Given a pattern/sequence, returns a vector of sequences.
            Other sequences are only rotations of the input sequence. """

        seq_len = len(sequence)
        rot_vector = [""] * seq_len
        seq = sequence * 2

        for i in range(seq_len):
            rot_vector[i] = seq[i: seq_len + i]

        return rot_vector

    def sorting(self, rot_matrix):

        """ Given a vector with rotated sequences, returns a vector of sorted sequences
            which is a Burrows-Wheeler's matrix. Sequences are sorted in a lexicographic order. """

        bw_matrix = sorted(rot_matrix)
        return bw_matrix

    def first_last_column(self, num_array, seq):

        """ Given a numerical array (offsets) and the original sequence, returns the first column. """

        first = []
        last = []
        for i, num in enumerate(num_array):
            first.append(seq[num])
            last.append(seq[(num - 1) % len(seq)])

        first = ''.join(first)
        last = ''.join(last)
        return first, last

    def first_column(self, Dict):

        """ Given a map of symbols with number of appearance, returns the first column as a dictionary.
            The first column is a map that has a range of rows prefixed by the symbol . """

        first = {}
        totc = 0
        temp = sorted(Dict.items())

        for c, count in temp:
            first[c] = (totc, totc + count)
            totc += count

        return first

    def last_column(self, bw_matrix):
        bwt = ''
        for s in bw_matrix:
            bwt += s[-1]
        return bwt

    def B_rank(self, sequence):

        """ Given the first/last column, returns a rank for each element in that sequence.
            Also returns a dictionary with a number of appearing for each element. """

        Dict = dict()
        ranks = []

        for s in sequence:

            if s not in Dict:
                Dict[s] = 0

            ranks.append(Dict[s])
            Dict[s] += 1

        return ranks, Dict

    def LF_mapping(self, last):

        """ Given the last column, returns an original sequence.
            This is a revers transform! """

        ranks, dicts = self.B_rank(last)
        first = self.first_column(dicts)

        i = 0
        original = '$'

        while last[i] != '$':
            s = last[i]
            original = s + original
            i = first[s][0] + ranks[i]

        return original

def suffix_array(sequence, period=1):
    ''' Given a pattern/sequence, returns only suffix-array's offsets.
        Taking sub-sequences from an incremental position and sorting them in lexicographic order. '''

    s = sorted([(sequence[i:], i) for i in range(0, len(sequence))])
    offsets = list(map(lambda x: x[1], s))

    if period > 1:
        new_sa = []
        for i in range(0, len(offsets)):
            if i % period == 0:
                new_sa.append(offsets[i])
        return new_sa
    else:
        return offsets

def cat_suffix_array(sa, period):
    new_sa = []
    for i in range(0, len(sa)):
        if i % period == 0:
            new_sa.append(sa[i])
    return new_sa

def tally_matrix(last, period=1):
    ''' Given the last column of the BW matrix, returns tally matrix as a dictionary.
        Keys are symbols and values are columns of a matrix where each column coresponds
        to a specific symbol. Tally matrix has size (sequence_length, n_symbols). Removed $ sign. '''

    counter = {}
    tally = {}

    for sym in last:
        if sym != '$':
            if sym not in counter:
                counter[sym] = 0
                tally[sym] = []

    for i, symb in enumerate(last):
        if symb != '$':

            counter[symb] += 1
            if (i % period) == 0:

                for key in counter.keys():
                    tally[key].append(counter[key])

        else:

            if (i % period) == 0:

                for key in counter.keys():
                    tally[key].append(counter[key])

    return tally, counter

def count(BWT, bwt, pattern):
    """ Given the Burrow-Wheeler's transform sequence and a pattern, returns a the first and last
        location in the first column. Still, we don't know what are the indeces for alignment! """

    L_rank, L_tot = BWT.B_rank(bwt)
    F = BWT.first_column(L_tot)
    i = len(pattern) - 2

    lower_limit, upper_limit = F[pattern[-1]]
    while i >= 0 and upper_limit > 1:

        search_symbol = pattern[i]
        j = lower_limit

        while j < upper_limit:

            if bwt[j] == search_symbol:
                lower_limit = F[search_symbol][0] + L_rank[j]
                break

            j += 1

        if j == upper_limit:
            lower_limit = upper_limit
            break

        upper_limit -= 1

        while bwt[upper_limit] != search_symbol:
            upper_limit -= 1

        upper_limit = F[search_symbol][0] + L_rank[upper_limit] + 1
        i -= 1

    if upper_limit == 1:
        print('The pattern does not appear in the sequence!')

    return lower_limit, upper_limit - 1


def fast_count(L_rank, L_tot, F, bwt, pattern, tally, period):
    """ Given a BW transform, pattern that we search for and tally matrix, returns the first and last
        location in the first column. Still, we don't know what are the indeces for alignment! """

    i = len(pattern) - 2

    lower_limit, upper_limit = F[pattern[-1]]
    lower_limit -= 1
    upper_limit -= 1

    while i >= 0:

        search_symbol = pattern[i]
        if lower_limit % period != 0:
            low_rank = find_rank(bwt, lower_limit, search_symbol, period, tally)
        else:
            low_rank = tally[search_symbol][lower_limit // period]


        if upper_limit % period != 0:
            high_rank = find_rank(bwt, upper_limit, search_symbol, period, tally)
        else:
            high_rank = tally[search_symbol][upper_limit // period]


        diff = high_rank - low_rank  # The number of appearance for current symbol
        if diff == 0:
            print('The pattern does not appear in the sequence!')
            break

        i -= 1

        lower_limit = F[search_symbol][0] + low_rank - 1
        upper_limit = F[search_symbol][0] + low_rank + diff - 1


    return lower_limit + 1, upper_limit


def locate(L_rank, L_tot, F, bwt, f_range, suffix_array, period):
    """ Given a complete suffix_array and a range from a first column, returns a list of position
        indices where a sequence and a pattern are aligned. """

    if f_range[1] - f_range[0] < 2:
        return []

    # Do this only if the suffix array is not optimized
    if period == 1:
        low_limit = f_range[0]
        up_limit = f_range[1]
        indices = suffix_array[low_limit:up_limit + 1]
        return indices

    len_seq = len(bwt)

    indices = []
    # Suffix array is optimized with a period
    for i in range(f_range[0], f_range[1] + 1):
        if i % period == 0:
            indices.append(suffix_array[i // period])
        else:
            row, count = i, 0
            while row % period != 0:
                count += 1
                s = bwt[row]
                row = F[s][0] + L_rank[row]
            indices.append((suffix_array[row // period] + count) % len_seq)
    return indices

def find_rank(bwt, row, symbol, period, cp):
    n_sym_hidden = 0
    row_start = row

    while (row % period) != 0:

        if bwt[row] == symbol:
            n_sym_hidden += 1

        row -= 1

    cp_pos = cp[symbol][row // period] + n_sym_hidden
    return cp_pos

def read_file(fpath):
    """ Givem a .txt file return a list of integers """

    file = open(fpath, "r")
    data = []

    while True:
        line = file.readline()
        if not line:
            break
        data.append(int(line))

    return data

def read_seq(fpath):
    sys.setrecursionlimit(10 ** 7)  # max depth of recursion
    threading.stack_size(2 ** 27)  # new thread will get stack of such size

    ca_file = open(fpath)
    ca_sequence = ca_file.read()
    ca_sequence = ca_sequence + '$'
    ca_sequence = re.sub(">.*\n?", "", ca_sequence)
    ca_sequence = ca_sequence.replace("\n", '')
    ca_sequence = ca_sequence.replace("N", '')

    return ca_sequence

def convert_size(size_bytes):
    if size_bytes == 0:
        return "0B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return "%s %s" % (s, size_name[i])

def find_indices(data, seq, pattern, sa_factor=1, tally_factor=1):
    bw = BWT()
    f, l = bw.first_last_column(data, seq)
    L_rank, L_tot = bw.B_rank(l)
    F = bw.first_column(L_tot)

    # Suffix array
    offsets = cat_suffix_array(data, period=sa_factor)

    # Tally matrix
    cp, tally = tally_matrix(l, tally_factor)

    # Start counting time
    start_time = timeit.default_timer()

    # Range where pattern occurs
    m = fast_count(L_rank, L_tot, F, l, pattern, cp, tally_factor)

    # Locations
    loc = locate(L_rank, L_tot, F, l, m, offsets, sa_factor)

    end_time = timeit.default_timer()

    # Print time
    print('Processing time needed: ', (end_time - start_time)*1000, '(in miliseconds)')
    count_bytes_tally = 0
    for elem in cp:
        count_bytes_tally+=sys.getsizeof(cp[elem])

    print('Size of objects:')
    print('Tally matrix: ', convert_size(count_bytes_tally))
    print('Suffix array: ', convert_size(sys.getsizeof(offsets)))

def testBWTandFM():
    test_obj = BWT()
    test_sample_keyword = ['abaaba$', 'mississippi$']
    test_sample_rotation = [
        ['abaaba$',
         'baaba$a',
         'aaba$ab',
         'aba$aba',
         'ba$abaa',
         'a$abaab',
         '$abaaba'],
        ['mississippi$'
         'ississippi$m',
         'ssissippi$mi',
         'sissippi$mis',
         'issippi$miss',
         'ssippi$missi',
         'sippi$missis',
         'ippi$mississ',
         'ppi$mississi',
         'pi$mississip',
         'i$mississipp',
         '$mississippi']
    ]

    test_sample_sorted = [
        ['$abaaba',
         'a$abaab',
         'aaba$ab',
         'aba$aba',
         'abaaba$',
         'ba$abaa',
         'baaba$a', ],
        ['$mississippi',
         'i$mississipp',
         'ippi$mississ',
         'issippi$miss',
         'ississippi$m',
         'mississippi$',
         'pi$mississip',
         'ppi$mississi',
         'sippi$missis',
         'sissippi$mis',
         'ssippi$missi',
         'ssissippi$mi']
    ]
    test_sample_first_column = [
        '$aaaabb',
        '$iiiimppssss'
    ]
    test_sample_last_column = [
        'abba$aa',
        'ipssm$pissii'
    ]
    test_sample_B_ranking = [
        [[1, 1, 2, 3, 2, 4, 1], {'a': 4, 'b': 2, '$': 1}],
        [[1, 1, 1, 2, 2, 3, 4, 3, 1, 2, 4, 1], {'m': 1, 'i': 4, 's': 4, 'p': 2, '$': 1}]
    ]

    test_sample_suffix_array = [
        [6, 5, 2, 3, 0, 4, 1],
        [11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]
    ]

    test_sample_tally_matrix_1 = [
        {'a': [1, 1, 1, 2, 2, 3, 4], 'b': [0, 1, 2, 2, 2, 2, 2]},
        {'i': [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4], 'p': [0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2],
         's': [0, 0, 1, 2, 2, 2, 2, 2, 3, 4, 4, 4], 'm': [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]}
    ]

    test_sample_tally_matrix_2 = [
        {'a': [1, 1, 2, 4], 'b': [0, 2, 2, 2]},
        {'i': [1, 1, 1, 1, 2, 3], 'p': [0, 1, 1, 2, 2, 2], 's': [0, 1, 2, 2, 3, 4], 'm': [0, 0, 1, 1, 1, 1]}
    ]

    test_sample_tally_matrix_4 = [
        {'a': [1, 2], 'b': [0, 2]},
        {'i': [1, 1, 2], 'p': [0, 1, 2], 's': [0, 2, 3], 'm': [0, 1, 1]}
    ]

    test_sample_pattern = [
        'aba',
        'si'
    ]
    test_sample_pattern_lookup = [
        (3, 4),
        (8, 9)
    ]
    testFailed = False

    # testing function rotation:
    #
    for i in range(len(test_sample_keyword)):
        if (test_obj.rotation(test_sample_keyword[i]) != test_sample_rotation[i]):
            testFailed = True

    if (testFailed == False):
        print('Rotation failed')
    else:
        print('Rotation passed')

    # testing function sorting:
    #
    for i in range(len(test_sample_keyword)):
        if (test_obj.sorting(test_sample_keyword[i]) != test_sample_rotation[i]):
            testFailed = True

    if (testFailed == False):
        print('Sorting failed')
    else:
        print('Sorting passed')

    # testing function first column extraction:
    #
    for i in range(len(test_sample_keyword)):
        if (test_obj.first_column(test_obj.B_rank(test_sample_keyword[i])[1]) != test_sample_first_column[i]):
            testFailed = True

    if (testFailed == False):
        print('First column extraction failed')
    else:
        print('First column extraction passed')

    # testing function last column extraction:
    #
    for i in range(len(test_sample_keyword)):
        if (test_obj.last_column(test_obj.B_rank(test_sample_keyword[i])[1]) != test_sample_last_column[i]):
            testFailed = True

    if (testFailed == False):
        print('Last column extraction failed')
    else:
        print('Last column extraction passed')

    # testing function B ranking:
    #
    for i in range(len(test_sample_keyword)):
        if (test_obj.B_rank(test_sample_keyword[i]) != test_sample_B_ranking[i]):
            testFailed = True

    if (testFailed == False):
        print('B ranking test failed')
    else:
        print('B ranking test passed')

    # testing function LF mapping:
    #
    for i in range(len(test_sample_keyword)):
        if (test_obj.LF_mapping(test_sample_last_column[i]) != test_sample_keyword[i]):
            testFailed = True
    if (testFailed == False):
        print('LF mapping test failed')
    else:
        print('LF mapping test passed')

    # testing function suffix array:
    #
    for i in range(len(test_sample_keyword)):
        if (list(suffix_array(test_sample_keyword[i])) != test_sample_suffix_array[i]):
            testFailed = True
    if (testFailed == False):
        print('SA test failed')
    else:
        print('SA test passed')

    # testing tally matrix function output with ranking 1:
    #
    for i in range(len(test_sample_keyword)):
        if (tally_matrix(test_sample_last_column[i])[0] != test_sample_tally_matrix_1[i]):
            testFailed = True
    if (testFailed == False):
        print('tally matrix function with ranking 1 test failed')
    else:
        print('tally matrix function with ranking 1 test passed')

    # testing tally matrix function output with ranking 2:
    #
    for i in range(len(test_sample_keyword)):
        if (tally_matrix(test_sample_last_column[i])[0] != test_sample_tally_matrix_2[i]):
            testFailed = True
    if (testFailed == False):
        print('tally matrix function with ranking 2 test failed')
    else:
        print('tally matrix function with ranking 2 test passed')

    # testing tally matrix function output with ranking 4:
    #
    for i in range(len(test_sample_keyword)):
        if (tally_matrix(test_sample_last_column[i])[0] != test_sample_tally_matrix_4[i]):
            testFailed = True
    if (testFailed == False):
        print('Tally matrix function with ranking 4 test failed')
    else:
        print('Tally matrix function with ranking 4 test passed')

    # testing basic pattern lookup:
    #
    for i in range(len(test_sample_keyword)):
        if (count(test_obj, test_sample_last_column[i], test_sample_pattern[i]) != test_sample_pattern_lookup[i]):
            testFailed = True
    if (testFailed == False):
        print('Pattern lookup test failed')
    else:
        print('Pattern lookup test passed')

    coffea_arabica_path = 'data/13443_ref_Cara_1.0_chr1c.fa'
    mus_pahari_path = 'data/10093_ref_PAHARI_EIJ_v1.1_chrX.fa'
    alligator_path = 'data/ami_ref_ASM28112v4_chrMT.fa'

    coffea_arabica_patterns = ['ATGCATG', 'TCTCTCTA','TTCACTACTCTCA']
    mus_pahari_patterns = ['ATGATG', 'CTCTCTA', 'TCACTACTCTCA']
    alligator_patterns = ['AGTCA', 'AACTCA', 'GTGCTTAG']

    sa_factors = [1, 4, 16, 64, 254]
    tally_factors = [1, 8, 32, 128, 512]

    # coffea arabica time consumption
    #

    data = read_file("data/sorted_ca.txt")
    seq = read_seq(coffea_arabica_path)

    for pattern in coffea_arabica_patterns:
        for sa_factor in sa_factors:
            for tally_factor in tally_factors:
                print('SA factor: ', sa_factor,
                      ' Tally factor: ', tally_factor,
                      ' Pattern:', pattern)
                find_indices(data, seq, pattern, sa_factor, tally_factor)

    # mus pahari time consumption
    #
    data = read_file("data/sorted_mp.txt")
    seq = read_seq(mus_pahari_path)

    for pattern in mus_pahari_patterns:
        for sa_factor in sa_factors:
            for tally_factor in tally_factors:
                print('SA factor: ', sa_factor,
                      ' Tally factor: ', tally_factor,
                      ' Pattern:', pattern)
                find_indices(data, seq, pattern, sa_factor, tally_factor)

    # alligator time consumption
    #

    data = read_file("data/sorted_al.txt")
    seq = read_seq(alligator_path)

    for pattern in alligator_patterns:
        for sa_factor in sa_factors:
            for tally_factor in tally_factors:
                print('SA factor: ', sa_factor,
                      ' Tally factor: ', tally_factor,
                      ' Pattern:', pattern)
                find_indices(data, seq, pattern, sa_factor, tally_factor)

# currently available indexes
# to build upon this, one should run rotation_with_indices and write indexes to separate rows in .txt file
# After that, add new row to this matrix
# Program will fail gracefully if no sorted element is found
#
available_sequences = [
    ["13443_ref_Cara_1.0_chr1c.fa", "sorted_ca.txt"],
    ["10093_ref_PAHARI_EIJ_v1.1_chrX.fa", "sorted_mp.txt"],
    ["ami_ref_ASM28112v4_chrMT.fa", "sorted_al.txt"]
]


if __name__ == "__main__":
    file = sys.argv[1]
    pattern = sys.argv[2]
    sa_factor = int(sys.argv[3])
    tally_factor = int(sys.argv[4])
    indexes = ""

    sa_factors = [1, 4, 16, 64, 254]
    tally_factors = [1, 8, 32, 128, 512]

    if sa_factor not in sa_factors:
        sys.exit("Unappropriate SA offset")

    if tally_factor not in tally_factors:
        sys.exit("Unappropriate tally matrix offes")

    for elem in available_sequences:
        if (elem[0] == file):
            indexes = elem[1]
            break

    if (indexes == ""):
        sys.exit("Indexes for file not found!")

    data = read_file("data/" + indexes)
    seq = read_seq("data/" + file)

    find_indices(data, seq, pattern, sa_factor, tally_factor)