import os
import shutil

ofile = open('input_dir_list.txt', 'r')
dir_list = ofile.readlines()
line_list = [1113,1808,1829,1830,1831,1832,1837,1839,1842,2339,2340,4125,4174,4177,4179,
             4180,4183,4194,4198,4213,4234,4245,4249,4266,4275,4277,4279,4289,4309,4316,
             4317,4318,4319,4320,4321,4322,4323,4324,4325,4326,4327,4328,4329,4330,4331,
             4690,4705,4879,4885,4888,4889,4911,4912,4913,4938,4939,4940,4942,4947,4948,
             4956,4962]
for iline in line_list:
    copy_to_dir = dir_list[iline-1].strip('\n') #These are 1 indexed for slurm
    copy_from_dir = copy_to_dir.replace('SNEC', 'SNEC2/SNEC')
    shutil.copy2(copy_from_dir, copy_to_dir)
    print('copied {} to {}'.format(copy_from_dir, copy_to_dir))

