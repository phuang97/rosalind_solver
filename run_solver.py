import pandas as pd
import os
from rosalind_solver.armory.rosalind_armory import RosalindArmory
from Bio import Entrez

ra = RosalindArmory()

# solve question 1
# print(ra.solve_introduction_to_the_bioinformatics_armory(
#    "data/q1_input.txt"))


# solve question 2
result = ra.solve_genbank_introduction(
    'Anthoxanthum', '2003/7/25', '2005/12/27')
print(result)
