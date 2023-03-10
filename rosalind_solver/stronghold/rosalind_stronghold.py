from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio import SeqUtils as ut


class RosalindStronghold():
    """
    each function of this class solves a question lsited on RosalindStronghold\n
    imports are called within each function to avoid excessive import list
    """

    NUCLEOTIDE = ['A', 'T', 'C', 'G']
    DNA_COMPLEMENT_DICT = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C"
    }

    DNA_CODON_DICT = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }

    def read_input_content(self, input_file_path: str):
        # read in the file
        with open(input_file_path) as f:
            content = f.read()
            content = content.strip()
        return content

    def write_solution_into_output(self, content: str, output_file_path: str, appending_mode: bool = False):
        # write content into output file
        if not appending_mode:
            with open(output_file_path, "w") as f:
                f.write(content)
                f.close()
        if appending_mode:
            with open(output_file_path, "a+") as f:
                f.write(content)
                f.close()

    def read_seq_from_inputfile(self, input_file_path: str) -> Seq:
        return Seq(self.read_input_content(input_file_path).upper())

    def solve_DNA(self, input_file_path: str):
        """
        Given: A DNA string s of length at most 1000 bp.
        Return: Four integers (separated by spaces) representing the respective number of times that the symbols
        'A', 'C', 'G', and 'T' occur in s. Note: You must provide your answer in the format shown in the sample output below.
        """
        import collections

        # read in the file
        seq = self.read_input_content(input_file_path)
        # first validate the input
        for nuc in seq:
            if nuc.upper() not in self.NUCLEOTIDE:
                raise ValueError('sequence contain none nucleotide string')
        # count the nucleotide frequency
        return dict(collections.Counter(seq))

    def solve_RNA(self, input_file_path: str):
        """
        Given: A DNA string t having length at most 1000 nt.
        Return: The transcribed RNA string of t.
        """
        dna_seq = self.read_input_content(input_file_path)
        # first validate the input
        for nuc in dna_seq:
            if nuc.upper() not in self.NUCLEOTIDE:
                raise ValueError('sequence contain none nucleotide string')
        # transcribe into RNA string
        return dna_seq.replace('T', 'U')

    def solve_REVC(self, input_file_path: str):
        """
        Given: A DNA string s of length at most 1000 bp.
        Return: The reverse complement sc of s.
        """
        # if not using biopython just use python's string.translate() with dict (complement)
        # and string[::-1] (reverse)
        dna_seq = Seq(self.read_input_content(input_file_path).upper())
        revc = dna_seq.reverse_complement()
        self.write_solution_into_output(
            str(revc), "solution/revc_solution.txt")

    def solve_FIB(self, input_file_path: str):
        """
        Given: Positive integers n and k.
        Return: The total number of rabbit pairs that will be present after n months,
        if we begin with 1 pair and in each generation,
        every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
        """
        content = self.read_input_content(input_file_path)
        try:
            n, k = content.split(' ')
        except:
            raise ValueError('input file content is not in expected format')

        # in the problem's demo:
        # F0 rabbit count = 1 pair (immature)
        # F1 rabbit count = 1 pair (able to reproduce)
        # F2 rabbit count = 1 pair (the pair from F1) + 3*1 pair (reproduced by F1, immature) = 4 pair
        # F3 rabbit count = 4 pair (mature) + 3 pair(immature) = 7 pair
        # F4 rabbit count = 7 + 3*4 = 19 pair
        # F5 rabbit count = 19 + 3*7 = 40 pair

        def fib(n, k):
            if n < 3:  # Rosalind test case treats F5 as month4, F4 as month 3, etc
                return 1
            else:
                return k*fib(n-2, k) + fib(n-1, k)
        solution = fib(int(n), int(k))
        self.write_solution_into_output(
            str(solution), 'solution/fib_solution.txt')

    def solve_GC(self, input_file_path: str):
        """
        Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
        Return: The ID of the string having the highest GC-content, followed by the GC-content of that string.
                Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated;
                please see the note on absolute error below.
        """
        # short cutting by using GC, if not simply use collection.Counter with iterator

        # setup counter variables
        highest_recorded_gc = 0
        seq_id = None
        # take a shortcut by using biopython's fasta iterator
        with open(input_file_path) as fasta_file:
            # FastaIterator reads each Fasta item and return individual item as SeqRecord class object
            for seq_record in FastaIO.FastaIterator(fasta_file):
                if ut.gc_fraction(seq_record.seq)*100 > highest_recorded_gc:
                    highest_recorded_gc = ut.gc_fraction(seq_record.seq)*100
                    seq_id = seq_record.id
        self.write_solution_into_output(
            f"{seq_id}\n{highest_recorded_gc}", "solution/gc_solution.txt")

    def solve_HAMM(self, input_file_path: str):
        content = self.read_input_content(input_file_path)
        s1 = content.split("\n")[0]
        s2 = content.split("\n")[1]
        # validate len s1 = len s2
        if len(s1) != len(s2):
            raise ValueError('input strings have different len')
        distance = 0
        for i in range(0, len(s1)):
            if s1[i] != s2[i]:
                distance += 1
        self.write_solution_into_output(
            f"{distance}", "solution/hamm_solution.txt")

    def solve_IPRB(self, input_file_path: str):
        """
        Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: k
               individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.
        Return: The probability that two randomly selected mating organisms will produce an individual possessing
                a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
        """
        content = self.read_input_content(input_file_path)
        k = int(content.split(' ')[0])
        m = int(content.split(' ')[1])
        n = int(content.split(' ')[2])
        total = k+m+n
        domdom = (k/total)*((k-1)/(total-1))
        hethet = (m/total)*((m-1)/(total-1)) * 0.75
        domrec = ((k/total)*(n/(total-1))) + ((n/total)*(k/(total-1)))
        domhet = ((k/total)*(m/(total-1))) + ((m/total)*(k/(total-1)))
        hetrec = (((m/total)*(n/(total-1))) + ((n/total)*(m/(total-1)))) * 0.5
        chance = domdom + hethet + domrec + domhet + hetrec
        self.write_solution_into_output(
            f'{chance}', 'solution/iprb_solution.txt')

    def solve_PROT(self, input_file_path: str):
        """translate rna into protein"""
        dna_seq = self.read_seq_from_inputfile(input_file_path)
        rna_seq = dna_seq.translate(to_stop=True)
        self.write_solution_into_output(
            f'{rna_seq}', 'solution/prot_solution.txt')

    def solve_SUBS(self, input_file_path: str):
        """
        Given: Two DNA strings s and t (each of length at most 1 kbp).
        Return: All locations of t as a substring of s.
        """
        content = self.read_input_content(input_file_path)
        seq = content.split('\n')[0].upper()
        sub_seq = content.split('\n')[1].upper()
        sub_len = len(sub_seq)
        solution = []
        for i in range(len(seq)-len(sub_seq)+1):
            if seq[i:i+sub_len] == sub_seq:
                solution.append(i+1)
        result = (' ').join([str(x) for x in solution])
        self.write_solution_into_output(result, 'solution/subs_solution.txt')

    def solve_CONS(self, input_file_path: str):
        """
        Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.
        Return: A consensus string and profile matrix for the collection.
                (If several possible consensus strings exist, then you may return any one of them.)
        """
        # initiate the dna string matrix as a dict
        matrix = dict()
        # read sequences into matrix
        with open(input_file_path) as fasta_file:
            for seq_record in FastaIO.FastaIterator(fasta_file):
                matrix[seq_record.id] = seq_record.seq
        # set up profile
        profile = dict()
        # seq length
        seq_len = len(list(matrix.values())[0])
        for nuc in ['A', 'T', 'C', 'G']:
            profile[nuc] = [0]*seq_len
        # calculate profile matrix
        for seq in matrix.values():
            for position, nuc in enumerate(seq):
                profile[nuc][position] += 1
        # calculate consensus
        consensus_list = list()
        for position in range(seq_len):
            count = 0
            for nuc in ['A', 'T', 'C', 'G']:
                if profile[nuc][position] >= count:
                    consensus = nuc
                    count = profile[nuc][position]
            consensus_list.append(consensus)
        # output result
        profile_result = str()
        for nuc in ['A', 'C', 'G', 'T']:
            profile_result += f"{nuc}: {' '.join([str(x) for x in profile[nuc]])}\n"
        consensus_result = f"{''.join(consensus_list)}\n"
        self.write_solution_into_output(
            f"{consensus_result}{profile_result}", 'solution/cons_solution.txt')

    def solve_FIBD(self, input_file_path: str):
        """
        Given: Positive integers n???100 and m???20.
        Return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.
        start with 1 pair of rabbits, each pair takes 1 month to grow and give birth of 1 pair of off-springs 1 month after grown
        """
        content = self.read_input_content(input_file_path)
        n = content.split(' ')[0]
        m = content.split(' ')[1]

        # this is too memory heavy
        # unfinished
        def fibd(n, m):
            if n in [1, 2]:
                return 1
            elif n < m:
                return fibd(n-1, m) + fibd(n-2, m)
            elif n == m:
                return fibd(n-1, m) + fibd(n-2, m)
            elif n > m:
                return fibd(n-1, m) + fibd(n-2, m) - fibd(n-m-2, m)

        def fibd2(n, m):
            if n in [1, 2]:
                return 1
            rabbit_count_list = list()
            rabbit_count_list.append(1)  # F0 or Month1
            rabbit_count_list.append(1)  # F1 or Month2
            for i in range(2, n):
                if i < m:
                    rabbit_count_list.append(
                        rabbit_count_list[i-1] + rabbit_count_list[i-2])
                if i == m:
                    rabbit_count_list.append(
                        rabbit_count_list[i-1] + rabbit_count_list[i-2] - rabbit_count_list[i-m])
                if i > m:
                    rabbit_count_list.append(
                        # rabbit born this month + rabbit we already have - rabbit died(which = rabbit born m-1 month ago)
                        # all rabbit we have m month ago = list[i-m-1] + list[i-m-2]
                        # since list[i-m-2] is the rabbit lived from previous months, we don't care about them
                        # list[i-m-1] is the rabbit born that month
                        rabbit_count_list[i-1] + rabbit_count_list[i-2] - rabbit_count_list[i-m-1])
            return rabbit_count_list[n-1]
        # output result
        result = fibd2(int(n), int(m))
        self.write_solution_into_output(
            f"{result}", 'solution/fibd_solution.txt')

    def solve_GRPH(self, input_file_path: str):
        """
        Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.
        Return: The adjacency list corresponding to O3. You may return edges in any order.
        """
        # set up a dict where key = seq.id, values = [head(k),tail(k)]. k=3.
        ligation_dict = dict()
        k = 3
        with open(input_file_path) as fasta_file:
            for seq_record in FastaIO.FastaIterator(fasta_file):
                ligation_dict[seq_record.id] = [
                    seq_record.seq[:k], seq_record.seq[-k:]]
        # alignment
        alignment = list()
        for key, value in ligation_dict.items():
            for key2, value2 in ligation_dict.items():
                if value[1] == value2[0] and value != value2:
                    alignment.append((key, key2))
        for item in alignment:
            print(f"{item[0]} {item[1]}")
            self.write_solution_into_output(
                f"{item[0]} {item[1]}\n", 'solution/grph_solution.txt', appending_mode=True)

    def solve_IEV(self, input_file_path: str):
        """
        Given: Six nonnegative integers, each of which does not exceed 20,000.
        The integers correspond to the number of couples in a population possessing each genotype pairing for a given factor.
        In order, the six given integers represent the number of couples having the following genotypes:
        AA-AA
        AA-Aa
        AA-aa
        Aa-Aa
        Aa-aa
        aa-aa
        Return: The expected number of offspring displaying the dominant phenotype in the next generation, u
        nder the assumption that every couple has exactly two offspring.
        """
        content = self.read_input_content(input_file_path)
        l = [int(x) for x in content.split(' ')]
        pr = (l[0] + l[1] + l[2] + l[3]*0.75 + l[4]*0.5)*2
        print(pr)
        self.write_solution_into_output(f"pr", "solution/iev_solution.txt")

    def solve_LCSM(self, input_file_path: str):
        """
        Given: A collection of k (k???100) DNA strings of length at most 1 kbp each in FASTA format.
        Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)
        """
        def find_shortest_seq(input_file_path):
            seq_len = 0
            with open(input_file_path) as fasta_file:
                for fasta in FastaIO.FastaIterator(fasta_file):
                    # first find the shortes tring
                    if seq_len == 0:
                        seq_len = len(fasta.seq)
                        seq = fasta.seq
                    elif len(fasta.seq) < seq_len:
                        seq_len = len(fasta.seq)
                        seq = fasta.seq
            return seq

        def list_all_possible_motif(seq):
            all_possible_motif_set = set()  # set is faster in searching b/c hash table
            for i in range(len(seq)):
                for x in range(i+1, len(seq)+1):
                    all_possible_motif_set.add(seq[i:x])
            return all_possible_motif_set

        def identify_longest_motif(input_file_path: str, all_possible_motif_set: set):
            # remove false motifs
            updated_motif_set = all_possible_motif_set.copy()
            with open(input_file_path) as fasta_file:
                for fasta in FastaIO.FastaIterator(fasta_file):
                    seq = fasta.seq  # so that fasta.seq is only calcuated 1 time per fasta
                    for motif in all_possible_motif_set:
                        if motif not in seq and motif in all_possible_motif_set:
                            try:
                                updated_motif_set.remove(motif)
                            except:
                                pass
                    # each iteration we will get a smaller possible motif set to search from
                    all_possible_motif_set = updated_motif_set.copy()
            # find the longest motif in motif list
            length = 0
            for motif in updated_motif_set:
                if len(motif) >= length:
                    longest_motif = motif
                    length = len(motif)
            return longest_motif
        shortest_seq = find_shortest_seq(input_file_path)
        all_possible_motif = list_all_possible_motif(shortest_seq)
        result = identify_longest_motif(input_file_path, all_possible_motif)
        self.write_solution_into_output(
            f"{result}", "solution/lcsm_solution.txt")
        print(result)

    def solve_LIA(self, input_file_path: str):
        """
        Given: Two positive integers k(k???7) and N(N???2k).
            In this problem, we begin with Tom, who in the 0th generation has genotype Aa Bb.
            Tom has two children in the 1st generation, each of whom has two children,
            and so on. Each organism always mates with an organism having genotype Aa Bb.
        Return: The probability that at least N
            Aa Bb organisms will belong to the k-th generation of Tom's family tree (don't count the Aa Bb mates at each level).
            Assume that Mendel's second law holds for the factors.
        """
        import math
        # AaBb x any genotype will have a 0.25 chance to produce AaBb offspring
        content = self.read_input_content(input_file_path)
        k = int(content.split(' ')[0])
        n = int(content.split(' ')[1])
        p = 2**k
        probability = 0
        # binomial distribution with cumulative distribution (because problem asks for 'at least N')
        for i in range(n, p + 1):
            prob = (math.factorial(p) / (math.factorial(i) *
                    math.factorial(p - i))) * (0.25**i) * (0.75**(p - i))
            probability += prob
        print(probability)
        self.write_solution_into_output(
            f"{probability}", 'solution/lia_solution.txt')

    def solve_MPRT(self, input_file_path: str):
        """
        Given: At most 15 UniProt Protein Database access IDs.
        Return: For each protein possessing the N-glycosylation motif,
                output its given access ID followed by a list of locations in the protein string where the motif can be found.
        """
        import urllib.request
        # use urlllib to extract given protein's sequence
        content = self.read_input_content(input_file_path)
        seq_dict = dict()
        # retrive protein sequence in fasta format
        for acc_id in content.split('\n'):
            if '_' in acc_id:
                protein_id = acc_id.split('_')[0]
                with urllib.request.urlopen(f"http://www.uniprot.org/uniprot/{protein_id}.fasta") as response:
                    res = response.read()
                seq_dict[acc_id] = str(res).split('\\n', 1)[
                    1].replace(r'\n', '')
            else:
                protein_id = acc_id
                with urllib.request.urlopen(f"http://www.uniprot.org/uniprot/{protein_id}.fasta") as response:
                    res = response.read()
                seq_dict[acc_id] = str(res).split('\\n', 1)[
                    1].replace(r'\n', '')
        # identify motif location
        motif_dict = dict()
        for acc_id, seq in seq_dict.items():
            # check if motif in seq
            # motif window = 4
            position_list = list()
            for i in range(len(seq)-4):
                # check if P in seq since N-glycosylation motif does not have P
                if 'P' in seq[i:i+4]:
                    pass
                elif seq[i] == 'N' and seq[i+2] in ['S', 'T']:
                    # rosalind takes position starting with 1
                    position_list.append(i+1)
            if len(position_list) != 0:
                motif_dict[acc_id] = position_list
        result = f''
        for acc_id, position in motif_dict.items():
            result = result + \
                f"{acc_id}\n{' '.join([str(x) for x in position])}\n"
        # output
        self.write_solution_into_output(result, 'solution/MPRT_solution.txt')

    def solve_MRNA(self, input_file_path: str):
        seq = self.read_input_content(input_file_path)
        possible_dna_seq = 0
        for aa in seq:
            possible_condon = len(
                [x for x in self.DNA_CODON_DICT.values() if x == aa])
            if possible_condon != 0:
                if possible_dna_seq == 0:
                    possible_dna_seq = possible_condon
                else:
                    possible_dna_seq = possible_dna_seq*possible_condon
        possible_dna_seq = possible_dna_seq*len(
            [x for x in self.DNA_CODON_DICT.values() if x == '_'])
        mod = possible_dna_seq % 1000000
        print(mod)

    def solve_ORF(self, input_file_path: str):
        # a reading frame must contain a stop codon
        # e.g. ATG-without a stop codon is not a reading frame
        def get_orf(dna):
            orf_list = list()
            dna = dna.upper()
            for frame in (0, 1, 2):
                for i in range(frame, len(dna), 3):
                    if dna[i:i+3] == 'ATG':  # i will be orf start position
                        for j in range(i, len(dna), 3):
                            if dna[j:j+3] in ['TAA', 'TAG', 'TGA']:  # j+3 is the end position
                                orf_list.append(
                                    dna[i:j+3])
                                break  # stop at the first stop codon
            return orf_list

        with open(input_file_path) as fasta_file:
            for seq_record in FastaIO.FastaIterator(fasta_file):
                input_sequence = Seq(seq_record.seq)
            # createa a list to store all final translated results
            protein_result_list = list()
            # original reading
            seq_list = get_orf(input_sequence)
            for seq in seq_list:
                protein_result_list.append(seq.translate(to_stop=True))
            # reverse complement reading
            seq_list2 = get_orf(
                input_sequence.reverse_complement())
            for seq in seq_list2:
                protein_result_list.append(seq.translate(to_stop=True))
        # drop duplicates in protein_result_list
        result = '\n'.join([str(x) for x in list(set(protein_result_list))])
        # output
        self.write_solution_into_output(result, "solution/orf_solution.txt")
        print(result)

    def solve_PERM(self, input_file_path: str):
        import itertools
        import math
        n = int(self.read_input_content(input_file_path))
        # calculate total possibility
        pr = math.factorial(n)

        def permutations(iterable, r=None):
            # this is the permutation function grabbed from itertools
            # there are 2 method provided by itertools, this one uses product() and modified to return all result in list
            result = list()
            pool = tuple(iterable)
            n = len(pool)
            r = n if r is None else r
            for indices in itertools.product(range(n), repeat=r):
                if len(set(indices)) == r:
                    result.append(tuple(pool[i] for i in indices))
            return result

        result = '\n'.join([str(x)
                           for x in permutations(range(1, n+1)) if x != None])
        for symbol in ['(', ')', ',']:
            result = result.replace(symbol, '')
        self.write_solution_into_output(
            f"{pr}\n{result}", 'solution/perm_solution.txt')
        print(f"{pr}\n{result}")
