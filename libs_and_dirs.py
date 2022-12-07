# General
import os
from shutil import rmtree
import logging
import time
import multiprocessing as mp
from enum import Enum

# Sequence comparison
from subprocess import run as subproc_run
from collections import Counter

# Archaic (Intended for deletion)
from Bio.SeqRecord import SeqRecord
from datetime import datetime

# Analysis
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


unallowed_filenames = ("con", "com", "prn", "aux", "nul", "lpt")

# Main IO
const_input_dir = "Main_Input\\"
file_1_dir = const_input_dir + "File 1\\"
file_2_dir = const_input_dir + "File 2\\"
const_output_dir = "Main_Output\\"
# Analysis IO
anlys_inp_dir = r"Analysis Input\\"
anlys_op_dir = r"Analysis Results\\"
anlys_memr_dir = r"Memory\\Analysis\\"
# Memory
memory_dir = "Memory\\MainProgram\\"
memr_inp_dir = "Memory\\MainProgram\\Input\\"
memr_op_dir = "Memory\\MainProgram\\Output\\"
temp_dir = "Memory\\DefaultTemp\\"
gene_lib_dir = temp_dir + "GeneLib\\"
blast_dir = "Memory\\blast\\bin\\"
# Logs
log_dir = "log\\"


# Setup files
# Main IO
if not os.path.exists(const_input_dir):
    os.makedirs(const_input_dir)
if not os.path.exists(file_1_dir):
    os.makedirs(file_1_dir)
if not os.path.exists(file_2_dir):
    os.makedirs(file_2_dir)
if not os.path.exists(const_output_dir):
    os.makedirs(const_output_dir)
# Analysis IO
if not os.path.exists(anlys_inp_dir):
    os.makedirs(anlys_inp_dir)
if not os.path.exists(anlys_op_dir):
    os.makedirs(anlys_op_dir)
if not os.path.exists(anlys_memr_dir):
    os.makedirs(anlys_memr_dir)
# Memory
if not os.path.exists(memory_dir):
    os.makedirs(memory_dir)
if not os.path.exists(memr_inp_dir):
    os.makedirs(memr_inp_dir)
if not os.path.exists(memr_op_dir):
    os.makedirs(memr_op_dir)

# Logs
if not os.path.exists(log_dir):
    os.makedirs(log_dir)


logging.basicConfig(filename=log_dir + str(time.time()) + ".log", filemode='w', level=logging.DEBUG,
                    format="(%(asctime)s.%(msecs)03d): %(message)s", datefmt='%M:%S')


class Gene:
    def __init__(self, name, id, start_pos, end_pos):
        self.name = name
        self.id = id
        self.start_pos = start_pos
        self.end_pos = end_pos


class Allele:
    def __init__(self, name: str, locations: tuple[str | int, str | int]):
        """

        :param name: Gene name.
        :param locations: Locations of gene in genomes. (<location in genome 1>, <location in genome 2>)
        """

        self.name = name
        self.locations = (None, int(locations[0]), int(locations[1]))
