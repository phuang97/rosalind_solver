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
