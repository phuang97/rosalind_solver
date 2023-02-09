import pandas as pd
import os
from rosalind_solver.armory.rosalind_armory import RosalindArmory
from rosalind_solver.stronghold.rosalind_stronghold import RosalindStronghold
from Bio import Entrez

ra = RosalindArmory()
rs = RosalindStronghold()

# solve question 1
# print(ra.solve_introduction_to_the_bioinformatics_armory(
#    "data/q1_input.txt"))

# solve question 2
# result = ra.solve_genbank_introduction(
#    'Anthoxanthum', '2003/7/25', '2005/12/27')
# print(result)

# solve RNA
# solution = rs.solve_RNA("data/rosalind_rna.txt")
# rs.write_solution_into_output(solution, 'solution/rna_solution.txt')

# solve REVC
# rs.solve_REVC("data/rosalind_revc.txt")

# solve FIB
# rs.solve_FIB("data/rosalind_fib.txt")

# solve GC
# rs.solve_GC('data/rosalind_gc.txt')

# solve HAMM
# rs.solve_HAMM("data/rosalind_hamm.txt")

# solve IPRb
# rs.solve_IPRB('data/rosalind_iprb.txt')

# solve PORT
# rs.solve_PROT('data/rosalind_prot.txt')

# solve SUBS
# rs.solve_SUBS('data/rosalind_subs.txt')

# solve CONS
# rs.solve_CONS('data/rosalind_cons.txt')

# solve FIBD
# rs.solve_FIBD('data/rosalind_fibd.txt')

# solve GRPH
# rs.solve_GRPH('data/rosalind_grph.txt')

# solve IEV
# rs.solve_IEV('data/rosalind_iev.txt')

# solve LCSM
# rs.solve_LCSM('data/rosalind_lcsm.txt')

# solve LIA
# rs.solve_LIA('data/rosalind_lia.txt')

# solve MPRT
rs.solve_MPRT('data/rosalind_mprt.txt')
