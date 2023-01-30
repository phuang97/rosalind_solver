from Bio.Seq import Seq
from Bio.SeqIO import FastaIO


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

    def read_input_content(self, input_file_path: str):
        # read in the file
        with open(input_file_path) as f:
            content = f.read()
            content = content.strip()
        return content

    def write_solution_into_output(self, content: str, output_file_path: str):
        # write content into output file
        with open(output_file_path, "w") as f:
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
        from Bio.SeqUtils import gc_fraction

        # setup counter variables
        highest_recorded_gc = 0
        seq_id = None
        # take a shortcut by using biopython's fasta iterator
        with open(input_file_path) as fasta_file:
            # FastaIterator reads each Fasta item and return individual item as SeqRecord class object
            for seq_record in FastaIO.FastaIterator(fasta_file):
                if gc_fraction(seq_record.seq)*100 > highest_recorded_gc:
                    highest_recorded_gc = gc_fraction(seq_record.seq)*100
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
