# this is the rosalind_solver

class RosalindSolver():
    """
    each function of this class solves a question lsited on RosalindSolver\n
    imports are called within each function to avoid excessive import list
    """

    NUCLEOTIDE = ['A', 'T', 'C', 'G']

    def solve_introduction_to_the_bioinformatics_armory(self, input_file: str):
        """
        Given: A DNA string s of length at most 1000 bp.
        Return: Four integers (separated by spaces) representing the respective number of times that the symbols 
        'A', 'C', 'G', and 'T' occur in s. Note: You must provide your answer in the format shown in the sample output below.
        """
        import collections
        # read in the file
        with open(input_file) as f:
            seq = f.read()
            seq = seq.strip()
        # first validate the input
        for nuc in seq:
            if nuc.upper() not in self.NUCLEOTIDE:
                raise ValueError('sequence contain none nucleotide string')
        # count the nucleotide frequency
        return dict(collections.Counter(seq))
