# this is the rosalind_solver

class RosalindArmory():
    """
    each function of this class solves a question lsited on RosalindArmory\n
    imports are called within each function to avoid excessive import list
    """

    NUCLEOTIDE = ['A', 'T', 'C', 'G']

    def solve_introduction_to_the_bioinformatics_armory(self, input_file_path: str):
        """
        Given: A DNA string s of length at most 1000 bp.
        Return: Four integers (separated by spaces) representing the respective number of times that the symbols 
        'A', 'C', 'G', and 'T' occur in s. Note: You must provide your answer in the format shown in the sample output below.
        """
        import collections
        # read in the file
        with open(input_file_path) as f:
            seq = f.read()
            seq = seq.strip()
        # first validate the input
        for nuc in seq:
            if nuc.upper() not in self.NUCLEOTIDE:
                raise ValueError('sequence contain none nucleotide string')
        # count the nucleotide frequency
        return dict(collections.Counter(seq))

    def solve_genbank_introduction(self, genus_name: str, date1: str, date2: str):
        from Bio import Entrez
        search_term = f'{genus_name}[All Fields]'
        Entrez.email = "phuang0114@gmail.com"
        ncbi_search_engine = Entrez.esearch(db='nucleotide',
                                            term=f'{search_term}',
                                            retmax='40',
                                            mindate=date1,
                                            maxdate=date2)
        result = Entrez.read(ncbi_search_engine)
        return result['Count']
