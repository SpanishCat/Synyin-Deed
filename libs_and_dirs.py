# General
import os
import logging
import time
import multiprocessing as mp
import re

# Specific
# OS
from shutil import rmtree
from subprocess import run as subproc_run
# Math
from copy import deepcopy
from itertools import groupby
from math import floor
from collections import Counter
from datetime import datetime
# Others
from enum import Enum
from Bio.SeqRecord import SeqRecord
from functools import partial

# Analysis
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


# Configurations
max_threads = 8

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
memory_dir = "Memory\\"
blast_dir = "Memory\\blast\\bin\\"
# Logs
log_dir = "log\\"

# Constants
no_gene_symbol = '*'
unallowed_filenames = {"con", "com", "prn", "aux", "nul", "lpt"}
gff_formats = ("gff", "gff3")
fasta_formats = ("fasta", "fa", "fna", "ffn", "faa", "frn")

# Config
conversion_1_output_format = ".genelist1"
"""Output format for the GenomeToGenelist conversion, i.e GFF3 to single genelist """


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
# Logs
if not os.path.exists(log_dir):
    os.makedirs(log_dir)


logging.basicConfig(filename=log_dir + str(time.time()) + ".log", filemode='w', level=logging.DEBUG,
                    format="(%(asctime)s.%(msecs)03d): %(message)s", datefmt='%M:%S')


def create_genelist_tuple(gff_text: str | tuple | list | set):
    """
    Creates a genelist.

    :param gff_text:
    :return: genelist: [[<genes in chrom. 1>], [<genes in chrom. 2>]...]
    """

    if type(gff_text) == str:
        gff_text = tuple(gff_text.splitlines())
    else:
        gff_text = tuple(gff_text)

    genelist = tuple(Gene(
        line.split(";Name=")[1][:line.split(";Name=")[1].find(';')],  # Name
        # line.split("GeneID:")[1].split(';')[0].split(',')[0],   # ID
        "",
        int(line.split()[3]),   # Start
        int(line.split()[4])    # End
    )
        for line in gff_text
        if not line.startswith('#') and line.split()[2] == "gene")
    # genelist.sort(key=lambda gene: gene.start_pos)

    # return tuple(genelist)
    return tuple(sorted(genelist, key=lambda gene: gene.start_pos))


def sort_key(e):
    if not e[0]:
        return int(e[1].start_pos) ** 2 + 1
    elif not e[1]:
        return int(e[0].start_pos) ** 2 + 2
    else:
        return int(e[0].start_pos) ** 2 + 3


# Enums DON'T DELETE!
class MatchBy(Enum):
    Name = 1
    ID = 2


class Order(Enum):
    List = 1
    Location = 2


# Genomics
class Gene:
    def __init__(self, name, id, start_pos, end_pos):
        self.name = name
        self.id = id
        self.start_pos = start_pos
        self.end_pos = end_pos


class Genome:
    def __init__(self, annotations: tuple | list | set | str, sequence: tuple | list | set | str,
                 place_in_list=None):
        """
        A genome object.

        :param annotations: GFF text.
        :param sequence: FASTA text.
        """

        # Adjust input
        if type(sequence) == str:
            sequence = sequence.splitlines()
        else:
            sequence = tuple(sequence)

        # Sequence
        name = " ".join(sequence[0].split()[1:3]).title()
        seq_rec = SeqRecord(seq="".join([line for line in sequence if not line.startswith('>')]))

        self.name = name
        self.annotations = annotations
        self.sequence = seq_rec.seq
        self.genelist = create_genelist_tuple(annotations)

    def replace(self, gene_dict: tuple[tuple[str, str]] | dict):
        if False:  # Remove if statement if you want to change the gene names from the gff itself and not just in the final genelist
            if type(self.annotations) == str:
                edited_annots = tuple(line for line in self.annotations.splitlines()
                                      if not line.startswith('#') and
                                      len(splits := line.split('\t')) == 9 and
                                      splits[2] == "gene"
                                      )
            else:
                edited_annots = tuple(line for line in self.annotations
                                      if not line.startswith('#') and
                                      len(splits := line.split('\t')) == 9 and
                                      splits[2] == "gene"
                                      )

            edited_annots = "\n".join(edited_annots)

        gene_dict_extracted = {gene for tup in gene_dict for gene in tup}
        edited_genelist = []

        for gene in self.genelist:
            if gene.name in gene_dict_extracted:
                dict_tup = next((tup for tup in gene_dict if gene.name in tup), None)
                edited_genelist.append(Gene(f"{dict_tup[0]}|{dict_tup[1]}", gene.id, gene.start_pos, gene.end_pos))

            else:
                edited_genelist.append(gene)

        self.genelist = edited_genelist


class PossibleMatch:
    def __init__(self, main_gene: Gene, candidates: Counter[str]):
        self.main_gene = main_gene
        self.candidates = candidates


class DefiniteMatch:
    def __init__(self, gene_1: Gene | str, gene_2: Gene | str):
        self.gene_1 = gene_1
        self.gene_2 = gene_2


class DoubleGenelist:

    def __init__(self, genome_1: Genome, genome_2: Genome):
        self.genomes = (genome_1.name, genome_2.name)

        temp_unmatched_1 = tuple((gene_1, None) for gene_1 in genome_1.genelist
                                 if len([gene_2 for gene_2 in genome_2.genelist if gene_2.name == gene_1.name]) == 0)

        temp_unmatched_2 = tuple((None, gene_2) for gene_2 in genome_2.genelist
                                 if len([gene_1 for gene_1 in genome_1.genelist if gene_2.name == gene_1.name]) == 0)

        temp_matched = tuple((gene_1, [gene_2 for gene_2 in genome_2.genelist if gene_2.name == gene_1.name][0])
                             for gene_1 in genome_1.genelist
                             if len([gene_2 for gene_2 in genome_2.genelist if gene_2.name == gene_1.name]) > 0)

        temp_genelist = (*temp_matched, *temp_unmatched_2, *temp_unmatched_1)

        self.list = tuple(sorted(temp_genelist, key=sort_key))


# Synteny
class Allele:
    def __init__(self, name: str, locations: tuple[str | int, str | int]):
        """
        A gene that exists in both genomes.

        :param name: Gene name.
        :param locations: Locations of gene in genomes. (<location in genome 1>, <location in genome 2>)
        """

        self.name = name

        # In case input is not numeric
        try:
            loc_1 = int(locations[0])
        except:
            loc_1 = None

        try:
            loc_2 = int(locations[1])
        except:
            loc_2 = None
        self.locations = (None, int(loc_1), int(loc_2))


# Misc.
class TextFile:
    def __init__(self, dir: str):
        """
        Gets the content of a text file.

        :param dir: Directory of file (including file name & suffix)
        """
        with open(dir, 'r') as inp_file:
            self.name = inp_file.name
            self.text = inp_file.read()
