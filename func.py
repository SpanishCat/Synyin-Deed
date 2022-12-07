import logging
from inspect import currentframe

from libs_and_dirs import *


# Up to user configuration
results_file_name = "SI Results.txt"
conversion_1_output_format = ".genelist1"
"""Output format for the GenomeToGenelist conversion, i.e GFF3 to single genelist """
fasta_formats = ("fasta", "fa", "fna", "ffn", "faa", "frn")


""" Dictionaries Structure:
    familiar_genes[<Gene Name>] = [<Location in Genome 1>, <Location in Genome 2>]
    unfamiliar_genes[<Gene Name>] = [<Genome Number (1 / 2)>, <Location in Genome>]
"""
files_dict = dict()
familiar_genes = dict()
unfamiliar_genes = dict()
counter_list = []


class Order(Enum):
    List = 1
    Location = 2


# Analysis Functions
def apply_algorithm(input_dir=anlys_inp_dir, output_dir=anlys_op_dir, acceptable_formats=("genelist2",)):
    """
    Applies Synteny Index algorithm.

    Input:
    Genelist2 file in format[
    1:genome_1,2:genome_2
    gene_name location_in_genome_1 location_in_genome_2
    gene_name location_in_genome_1 location_in_genome_2
    gene_name location_in_genome_1 location_in_genome_2...]

    Output:
    Synteny Index of the two genomes

    :arg input_dir: Input directory
    :arg output_dir: Output directory
    :arg acceptable_formats: File formats the function accepts. Formats must be lowercase
    :return: None
    """

    synteny_indexes_of_genes = []

    # Get input
    filename = show_and_get_files(1,
                                  "Choose a genelist2 file to compare using SI:\n",
                                  input_dir_=anlys_inp_dir, acceptable_formats=acceptable_formats)

    while True:
        try:
            k = int(input("k = ?\n"))
            break
        except:
            print("Invalid Input! k should be a natural number")

    # Open file
    with open(input_dir + filename, 'r') as inp_file:
        inp_text = inp_file.read().splitlines()
    logging.debug(f'Text file: "{filename}"')

    genome_1_name = inp_text[0].split(':')[1].split(',')[0]
    genome_2_name = inp_text[0].split(':')[2]

    # Remove comments
    # todo: Put genome names in ## in genelist2 format
    inp_text.pop(0)
    inp_text = [line for line in inp_text if not line.startswith('#')]

    # fixme: every neighborhood has a different length
    def make_nbrhood(center_gene: Allele, nbrhood_num: int) -> list:
        """
        Create neighborhood for given gene.

        :param center_gene: The gene around which the neighborhood is built. center_gene = [<gene name>, <location in neighborhood 1>, <location in neighborhood 2>].
        :param nbrhood_num: 1 or 2. Number of neighborhood. 1: Left in genelist, 2: Right in genelist.
        :return: Neighborhood.
        """

        center_gene_num = center_gene.locations[nbrhood_num]   # Just for readability
        logging.debug(f"Center Gene: {center_gene.name} ({center_gene_num})\n\n")

        # Set neighborhood
        nbrhood = [line.split()[0] for line in inp_text if
                   # 1. If gene is homozygous (exists in both genomes)
                   line.split()[nbrhood_num] != "X" and

                   # 2. If gene is either a:
                   # - Direct Neighbor (e.g 2 and 3 in [1,2,3,4,5])
                   (center_gene_num - k <= int(line.split()[nbrhood_num]) <= center_gene_num + k or
                    # - Cyclical Neighbor (e.g 1 and 5 in [1,2,3,4,5])
                    center_gene_num + k - len(inp_text) >= int(line.split()[nbrhood_num]) or
                    center_gene_num - k + len(inp_text) <= int(line.split()[nbrhood_num]))
                   ]

        for line in inp_text:
            try:
                logging.debug(f"\n\n\n{line.split()[0]}:\n line.split()[{nbrhood_num}] != 'X': {line.split()[nbrhood_num] != 'X'}")
                logging.debug(f"{center_gene_num} - k({k}) <= int(line.split()[{nbrhood_num}]) <= {center_gene_num} + k({k}): {(center_gene_num - k <= int(line.split()[nbrhood_num]) <= center_gene_num + k)}")
                logging.debug(f"{center_gene_num} + k({k}) - {len(inp_text)} >= int(line.split()[{nbrhood_num}]): {center_gene_num + k - len(inp_text) >= int(line.split()[nbrhood_num])}")
                logging.debug(f"{center_gene_num} - k({k}) + {len(inp_text)} <= int(line.split()[{nbrhood_num}]): {center_gene_num - k + len(inp_text) <= int(line.split()[nbrhood_num])}")

            except: pass

        if center_gene.name in nbrhood:
            nbrhood.remove(center_gene.name)

        # logging.debug(f"\nnbrhood: " + str(nbrhood).replace("',", f"': {}").replace(', ', ',\n'))
        logging.debug(f"\n\nNeighborhood {nbrhood_num}: ")
        for i, gene in enumerate(nbrhood):
            logging.debug(f"{gene}: {i+1}")
        return nbrhood

    # Apply SI to each gene
    text_no_asterisks = [line for line in inp_text if "*" not in line]
    try:
        for line_num, current_gene_line in enumerate(text_no_asterisks):
            if line_num % 100 == 0:
                print(f"\r{int(100 * ((line_num + 1) / len(text_no_asterisks)))}%", end="")

            # Define current gene
            gene_in_check = Allele(current_gene_line.split()[0],  # Name
                                   (current_gene_line.split()[1],  # Location in genome 1
                                   current_gene_line.split()[2]))  # Location in genome 2

            # Create neighborhoods for current gene
            neighborhoods = [None, make_nbrhood(gene_in_check, 1), make_nbrhood(gene_in_check, 2)]
            if line_num % 1000:
                logging.debug(f"Neighborhoods for gene {gene_in_check}:\n"
                              f"Neighborhood 1: {neighborhoods[1]}\n Neighborhood 2: {neighborhoods[2]}")

            # Find common genes
            intersection = []
            for gene_1, gene_2 in zip(neighborhoods[1], neighborhoods[2], strict=True):
                print(f"Gene: {gene_1}; {len(neighborhoods[1])} == {len(neighborhoods[2])}")
                if (gene_1 in neighborhoods[2] and not gene_1.startswith("*")) or \
                   gene_1.startswith("*") and (gene_2.startswith("*")):
                    intersection.append(gene_1)
            logging.debug("Intersection: " + str(intersection))

            # SI formula
            x = (1 / (k * 2)) * len(intersection)
            synteny_indexes_of_genes.append(x)
            logging.debug(f"\nSI for gene {gene_in_check.name}: {x}\n")

    except ValueError:
        raise Exception(f"Neighborhoods are of different lengths ({len(neighborhoods[1])}, {len(neighborhoods[2])})")

    # Average out gene SI's to find genome's SI
    if len(synteny_indexes_of_genes) > 0:
        synteny_index_of_genome = sum(synteny_indexes_of_genes) / len(synteny_indexes_of_genes)
    else:
        synteny_index_of_genome = None
        print("\nERROR: No recognized genes found. SI is NONE")
    print(f"SI for genomes {genome_1_name} and {genome_2_name} = {synteny_index_of_genome}\n")

    # Write results in results file
    if not os.path.exists(output_dir + results_file_name):
        results_file = open(output_dir + results_file_name, 'w')
    else:
        results_file = open(output_dir + results_file_name, 'a')
        results_file.write("\n\n")

    results_file.write("Comparison results for genomes " + genome_1_name + " and " + genome_2_name + ": " +
                       "SI = " + str(synteny_index_of_genome))
    results_file.close()


def find_genome(dir):
    """
    Find GFF3 and Fasta in specified directory
    :param dir:
    :return:
    """
    gff3, fasta = None, None
    for filename in [file for (path, dirs, files) in os.walk(dir) for file in files]:

        if filename.split('.')[-1].lower() in fasta_formats:
            fasta = filename
            continue

        else:
            with open(dir + filename, 'r') as file:
                text = file.readlines()

            if text[0] == "##gff-version 3\n":
                gff3 = filename
                continue

        if gff3 and fasta: break

    return gff3, fasta


def create_genelist1(text_name, gene_type="Name", input_dir=const_input_dir, output_dir=const_output_dir, order=Order.Location):
    """
    Yes INDEED
    :param order: Order.Location/Order.List.
    List: number the genes according to location in genelist, e.g:
    PAU8 1 1
    YAL067W-A 1 2
    SEO1 1 3

    Location: number the genes according to location in genome, e.g:
    PAU8 1 2169
    YAL067W-A 1 2707
    SEO1 1 9016
    :return:
    """

    # Open file
    with open(input_dir + text_name, 'r') as input_file:
        input_text = input_file.read().splitlines()
    genes = []
    chromosomes = []

    # Make GENELIST file
    if input_text[0] == "##gff-version 3":
        logging.debug("File is a GFF3 file. Proceeding")
        output_file = open(output_dir + text_name.removesuffix("." + text_name.split(".")[-1]) + conversion_1_output_format, 'w')

        # Read data file
        for line in input_text:
            # Chromosomes
            if len(line.split()) >= 8 and line.split()[2] == "region":
                chr = line[line.find("chromosome=") + len("chromosome=")]

                if not chromosomes.__contains__(chr):
                    logging.debug("New chromosome found: " + str(chr))
                    chromosomes.append(chr)

            # Setup Gene
            if len(line.split()) >= 8 and line.split()[2] == "gene":
                this_gene = Gene("", "", 0, 0)
                this_gene.name = line[line.find(";Name=") + 6:(line.find(";", line.find(";Name=") + 6, line.find("\n")))]
                this_gene.id = line[line.find("ID=") + 3:(line.find(";", line.find("ID=") + 3, line.find("\n")))]
                # this_gene.position = gene_num
                this_gene.start_pos = int(line.split()[3])

                genes.append(this_gene)

            # logging.debug("Genes: " + str(len(genes)))

        # Upon finishing:
        # i = 1
        # logging.debug("Chromosomes in genome: " + str(chromosomes))

        if order == Order.List:
            if gene_type == "Name":
                for num, gene in enumerate(genes):
                    output_file.write(f"{gene.name} {num}\n")
                    # i += 1

            elif gene_type == "ID":
                for num, gene in enumerate(genes):
                    output_file.write(f"{gene.id} {num} \n")
                    # i += 1

        elif order == Order.Location:
            if gene_type == "Name":
                for num, gene in enumerate(genes):
                    output_file.write(f"{gene.name} {gene.start_pos}\n")

            elif gene_type == "ID":
                for num, gene in enumerate(genes):
                    output_file.write(f"{gene.id} {gene.start_pos} \n")

        output_file.close()

    elif input_text[0].startswith("##gff-version"):
        logging.error("Unsupported format version (GFF " + str(input_text[0].split()[1]) + ")")

    else:
        logging.error("Unsupported file format")


def merge_genelist(filename_1, filename_2, input_dir=const_input_dir, output_dir=const_output_dir, asterisks=False):
    def sort_key(e):
        if e[2] == "X":
            return int(e[1]) ** 2 + 1
        elif e[1] == "X":
            return int(e[2]) ** 2 + 2
        else:
            return int(e[1]) ** 2 + 3

    # Apply Search
    # Open files
    file_1 = open(input_dir + str(filename_1), 'r')
    file_2 = open(input_dir + str(filename_2), 'r')
    text_1 = file_1.read().splitlines()
    text_2 = file_2.read().splitlines()
    # Get filenames without format type
    genelist_name_1 = filename_1.removesuffix("." + filename_1.split('.')[-1])
    genelist_name_2 = filename_2.removesuffix("." + filename_2.split('.')[-1])
    output_file = open(output_dir + "Merged-" + genelist_name_1 + "-" + genelist_name_2 + ".genelist2", 'w')
    output_txt = ""
    unfamiliar_genes_file = open(output_dir + "Unfamiliar Genes.genelist2", 'w')
    unfamiliar_genes_txt = ""

    output_txt += f"##Genomes: {genelist_name_1} {genelist_name_2}\n"
    # processed = 0
    seqs: list = []
    wrong_asters = 0

    if not asterisks:
        # for line_num, line_1 in enumerate(text_1):
        #     genename_1 = line_1.split()[0]
        #     logging.debug("Processing gene " + genename_1)
        #     is_homozygous = False
        #
        #     for line_2 in text_2:
        #         if f"{genename_1} " in line_2:
        #
        #             familiar_genes[genename_1] = [line_1.split()[1], line_2.split()[1]]
        #             is_homozygous = True
        #             break
        #
        #     # Check Unidentified Genes from Genome 1
        #     if not is_homozygous:
        #         unfamiliar_genes[genename_1] = [1, line_1.split()[1]]
        #
        #     # processed += 1
        #     print("\rProcessing " + str(int(line_num / len(text_1) * 100)) + "%", end="")

        unfamiliar_genes: dict = {line_1.split()[0]: [1, line_1.split()[1]] for line_1 in text_1
                                  if f"{line_1.split()[0]} " not in text_2}

        familiar_genes: dict = {line_1.split()[0]: [line_1.split()[1],
                                [line_2 for line_2 in text_2
                                if f"{line_1.split()[0]} " in line_2][0]]

                                for line_1 in text_1 if f"{line_1.split()[0]} " in text_2
                                }

    # If asterisks:
    else:
        familiar_genes: dict = {line_1.split()[0]: [line_1.split()[1],
                                [line_2 for line_2 in text_2
                                if f"{line_1.split()[0]} " in line_2][0]]

                                for line_1 in text_1 if f"{line_1.split()[0]} " in text_2}

        print('Non-asterisks 1: ' + str([(line_1.split()[0], line_1.split()[1], [line_2 for line_2 in text_2
              if f"{line_1.split()[0]} " in line_2][0])

             for line_1 in text_1 if f"{line_1.split()[0]} " in text_2]))

        seqs = [
            # Matched genes
            *[(line_1.split()[0], line_1.split()[1], [line_2 for line_2 in text_2
              if f"{line_1.split()[0]} " in line_2][0])

             for line_1 in text_1 if f"{line_1.split()[0]} " in text_2],

            # Unmatched genes
            # Genome 1
            *[(f"*{line_1.split()[0]}", line_1.split()[1], "X") for line_1 in text_1
            if line_1.split()[0] not in text_2],
            # Genome 2
            *[(f"*{line_2.split()[0]}", "X", line_2.split()[1]) for line_2 in text_2
             if line_2.split()[0] not in text_1]
        ]

        print(f"Non-asterisks: {[x for x in seqs if not '*' in x[0]]}")
        wrong_asters = len([aster for aster in seqs if '*' in aster[0] and '|' in aster[0]])

        # # Add unmatched genes
        # seqs.append((f"*{line_1.split()[0]}", line_1.split()[1], "X") for line_1 in text_1
        #             if line_1.split()[0] not in text_2)
        # seqs.append((f"*{line_2.split()[0]}", "X", line_2.split()[1]) for line_2 in text_2
        #             if line_2.split()[0] not in text_1)


        # for line_1 in text_1:
        #     genename_1 = line_1.split()[0]
        #     logging.debug("Processing gene " + genename_1)
        #
        #     for line_2 in text_2:
        #         if f"{genename_1} " in line_2:
        #
        #             familiar_genes[genename_1] = [line_1.split()[1], line_2.split()[1]]
        #             seqs.append((genename_1, line_1.split()[1], line_2.split()[1]))
        #             break
        #
        # # Check Unidentified Genes from Genome 1
        # else:
        #     seqs.append((f"*{genename_1}", line_1.split()[1], "X"))
        #     if "|" in genename_1:
        #         print(f"Wrong asterisk1: {genename_1}")
        #         wrong_asters += 1
        #
        # processed += 1
        # print("\rProcessing " + str(int(processed / len(text_1) * 100)) + "%", end="")

    # Check Unidentified Genes from Genome 2
    # for line_2 in text_2:
    #     genename_2 = line_2.split()[0]
    #     if not asterisks and genename_2 not in unfamiliar_genes and genename_2 not in familiar_genes:
    #         unfamiliar_genes[genename_2] = [2, line_2.split()[1]]
    #
    #     if asterisks and genename_2 not in [x[0] for x in seqs]:
    #         seqs.append((f"*{genename_2}", "X", line_2.split()[1]))
    #         if "|" in genename_2:
    #             print(f"Wrong asterisk2: {genename_2}")
    #             wrong_asters += 1

    if not asterisks:
        output_txt += str(familiar_genes) \
                          .removeprefix('{')\
                          .removesuffix('}')\
                          .replace(":", "")\
                          .replace("'", "")\
                          .replace("[", "")\
                          .replace("], ", "\n")\
                          .replace("]", "")\
                          .replace(",", "")

        unfamiliar_genes_txt += (str(unfamiliar_genes)
                                    .removeprefix('{')
                                    .removesuffix('}')
                                    .replace("[", "")
                                    .replace("], ", "\n")
                                    .replace("'", "")
                                    .replace(",", "")
                                    )

    elif asterisks:
        seqs.sort(key=sort_key)
        output_txt += (str(seqs)
                          .removeprefix('[')
                          .removesuffix(']')
                          .replace("), ", "\n")
                          .replace("'", "")
                          .replace(",", "")
                          .replace("(", "")
                          .replace(")", "")
                          )

    output_file.write(output_txt)
    unfamiliar_genes_file.write(unfamiliar_genes_txt)

    # Close files
    file_1.close()
    file_2.close()
    output_file.close()
    unfamiliar_genes_file.close()

    if asterisks:
        print(f"Wrote {wrong_asters} wrong asterisks ({int(100 * (wrong_asters / (len(text_1) + len(text_2))))}%)")

    if len(unfamiliar_genes_txt) > 0:
        # print("do_blast = True")
        return True
    else:
        # print("do_blast = False")
        return False


# "FastaSeperatorByGene" Definitions
def separate_by_gene(gff3, fasta, genelist_output_format=".genelist1", output_filename="", output_format=".fasta",
                     one_outp_file=True, specific_genes_only=False, specific_genes=(),
                     input_dir=const_input_dir, output_dir=const_output_dir,
                     memory_dir=r"Memory\\DefaultTemp\\"):

    """

    :param gff3: Given GFF3 file
    :param fasta: Given FASTA file
    :param one_outp_file: If True: Output all sequences as different sections in a single FASTA file
                          If False: Output every sequence as its own FASTA file
    :param genelist_output_format: What GENELIST format to use for the separation
    :param output_filename:
    :param output_format:
    :param specific_genes_only:
    :param specific_genes:
    :param input_dir: Input directory
    :param output_dir: Output directory
    :param memory_dir: (Temporary) Memory directory
    :return:
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    def write_genelist(text_name):
        input_file = open(input_dir + text_name, 'r')
        output_file = open(memory_dir + text_name.removesuffix("." + text_name.split(".")[-1]) + genelist_output_format, 'w')
        genes = []

        # Read the data file
        while True:
            line = input_file.readline()
            if line:

                # GFF3
                if text_name.split(".")[-1].lower() == "gff3":
                    # Create Gene
                    if len(line.split()) >= 8 and line.split()[2] == "gene":
                        this_gene = Gene("", "", 0, 0)

                        # Write Gene Data
                        this_gene.name = line[line.find(";Name=") + 6:
                                              (line.find(";", line.find(";Name=") + 6, line.find("\n")))]
                        this_gene.id = line[line.find("ID=") + 3:
                                              (line.find(";", line.find("ID=") + 3, line.find("\n")))]
                        this_gene.start_pos = line.split()[3]
                        this_gene.end_pos = line.split()[4]

                        # Exceptions
                        for exception in unallowed_filenames:
                            if this_gene.name.lower().startswith(exception):
                                this_gene.name = this_gene.name + "_"

                        # Add to list
                        if not specific_genes_only or \
                                (specific_genes_only and specific_genes.__contains__(this_gene.name)):
                            genes.append(this_gene)

            # Upon finishing:
            else:
                # GFF3
                if text_name.split(".")[-1].lower() == "gff3":
                    i = 1
                    for gene in genes:
                        output_file.write(gene.name + " " + str(gene.start_pos) + " " + str(gene.end_pos) + "\n")
                        i += 1
                break

    def write_fasta(genelist_file_name):
        genelist_file = open(memory_dir + genelist_file_name, 'r')
        genelist_text = genelist_file.read().splitlines()
        fasta_file = open(input_dir + fasta, 'r')
        fasta_text = fasta_file.read()
        # fasta_text = fasta_text.removeprefix(fasta_text.splitlines()[0] + "\n")
        # logging.debug("Fasta text length: " + str(len(fasta_text)))

        # Remove info lines
        for line in fasta_text:
            if line.startswith(">"):
                fasta_text.replace(line, "")
        fasta_text = fasta_text.replace("\n", "").replace("\r", "")
        # logging.debug("Fasta text length: " + str(len(fasta_text)))

        # Cleanup
        for paths, dirs, files in os.walk(output_dir):
            for file in files:
                os.remove(output_dir + file)

        for line in genelist_text:
            this_gene = Gene(line.split()[0], None, int(line.split()[1]), int(line.split()[2]))

            # In case windows can't handle the file name
            try:
                gene_fa_file = open(output_dir + this_gene.name + ".FASTA", 'w')

            except:
                logging.debug("Failed to open " + this_gene.name)

                for exception in unallowed_filenames:
                    if this_gene.name.lower().startswith(exception):
                        this_gene.name = str(this_gene.name) + "_"
                    gene_fa_file = open(output_dir + this_gene.name + ".FASTA", 'w')

            gene_fa_file.write(fasta_text[this_gene.start_pos - 1:this_gene.end_pos])
            # gene_fa_file.write(fasta_text[this_gene.end_pos + 1:
            #                               this_gene.end_pos + len(fasta_text[
            #                                     this_gene.start_pos - 1:this_gene.end_pos].splitlines()) - 1])

            logging.debug(this_gene.name + " contains: " +
                          str(len(fasta_text[this_gene.start_pos - 1:this_gene.end_pos].splitlines())) + " lines")
            gene_fa_file.close()

        # Close files
        genelist_file.close()
        fasta_file.close()

    def write_as_1_fa(genelist_file_name):
        genelist_file = open(memory_dir + genelist_file_name, 'r')
        genelist_text = genelist_file.read().splitlines()
        inp_fa_file = open(input_dir + fasta, 'r')
        inp_fa_text = inp_fa_file.read()
        outp_txt =""

        if output_filename == "":
            outp_file = open(f"{output_dir}Genes-{fasta}", 'w')
        else:
            outp_file = open(f"{output_dir}{output_filename}{output_format}", 'w')

        # Remove info lines
        for line in inp_fa_text.splitlines():
            if line.startswith(">") or line.startswith(";"):
                inp_fa_text = inp_fa_text.replace(line, "")
        inp_fa_text = inp_fa_text.replace("\n", "").replace("\r", "")

        # Write sequence
        for line in genelist_text:
            this_gene = Gene(line.split()[0], None, int(line.split()[1]), int(line.split()[2]))
            outp_txt += f">{this_gene.name}\n" \
                        f"{inp_fa_text[this_gene.start_pos:this_gene.end_pos]}\n"

        outp_file.write(outp_txt)

        # Close files
        genelist_file.close()
        inp_fa_file.close()
        outp_file.close()

    write_genelist(gff3)

    if one_outp_file:
        write_as_1_fa(gff3.lower().removesuffix(".gff3").capitalize() + genelist_output_format)
    else:
        write_fasta(gff3.lower().removesuffix(".gff3").capitalize() + genelist_output_format)


def merge_txt_files(input_dir=const_input_dir, output_dir=const_output_dir,
                    list_filename="list", list_file_suffix=".txt"):
    print("Add spaces between files? (y/n)")
    bool = str(input())
    if bool.lower() == "y" or bool.lower() == "yes" or bool.lower() == "t" or bool.lower() == "true":
        add_space = True
    elif bool.lower() == "n" or bool.lower() == "no" or bool.lower() == "f" or bool.lower() == "false":
        add_space = False

    list_read = open(list_filename + list_file_suffix, 'r')

    list_text = list_read.read()
    print(list_text.find("\n"))
    list_lines = list_text.splitlines()
    output_file_name = list_lines[0]
    output_file = open(output_dir + output_file_name, 'w')
    print(list_lines)

    print(input_dir[0:-2])
    for path, dirs, files in os.walk(input_dir[0:-2]):
        remaining_files = files
        i = 1
        while len(remaining_files) > 0:
            for file in files:
                print("files: " + file[0:file.rfind('.')])
                # print("line: " + list_lines[i][0:list_lines[i].rfind('.')])
                if file[0:file.rfind('.')] == list_lines[i][0:list_lines[i].rfind('.')]:
                    current_file = open(input_dir + file, 'r')
                    output_file.write(current_file.read())

                    if add_space:
                        output_file.write("\n")

                    current_file.close()
                    print("File " + file + " Done")
                    remaining_files.remove(file)
                    print(remaining_files)
                    if i < len(list_lines):
                        i += 1
    # output_file = output_file.read().removesuffix("10\n")
    # output_file.write("mdaskl;wq")
    output_file.close()
    list_read.close()


def blast(is_protein=False, evalue=1e-50,
          input_dir=r"Main_Input\\File 1\\", output_dir=r"Main_Output\\",
          blast_dir="Memory\\blast\\bin\\") -> None:
    """
    Run BLAST from on given genomes.
    \nis_protein: True = Match protein sequences. False = Match nucleotide sequences.
    """

    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)

    path = os.path.dirname(os.path.abspath(__file__)) + '\\' + blast_dir

    # Create database (CMD)
    if is_protein:
        print(subproc_run(["makeblastdb",
                              "-in", "Genome1.fasta",
                              "-out", "Genome1",
                              "-dbtype", "prot"],
                             text=True, shell=True, cwd=path))
    else:
        print(subproc_run(["makeblastdb",
                              "-in", "Genome1.fasta",
                              "-out", "Genome1",
                              "-dbtype", "nucl"],
                             text=True, shell=True, cwd=path))

    # Align (CMD)
    if is_protein:
        print(subproc_run(["blastp",
                              "-db", "Genome1",
                              "-query", "Genome2.fasta",
                              "-out", "results.txt",
                              "-outfmt", "6",
                              "-evalue", str(evalue)],
                             text=True, shell=True, cwd=path))
    else:
        print(subproc_run(["blastn",
                              "-db", "Genome1",
                              "-query", "Genome2.fasta",
                              "-out", "results.txt",
                              "-outfmt", "6",
                              "-evalue", str(evalue),
                              # "-subject_besthit",
                              "-max_target_seqs", "5",
                              "-task", "megablast"],
                             text=True, shell=True, cwd=path))

    # todo: Get read of shell=True
    # todo: Add to blast alignment parameters the following:
    """ - Open gap
        - Gap penalty
        - nucl match reward
        - nucl mismatch penalty
        - megablast V
        - limit output V
    """


def fix_parallels(num, parallel, length) -> tuple:
    """
    This function's code runs for every parallel
    :return: temp_list, to_delete
    """
    if mp.current_process().name.split('-')[-1] == '1' and num % 3 == 0:
        print(f"\rFixing parallels: {num} ({int((num / length) * 100)}%)", end="")

    temp_list_obj = [{dict_item: dic[1][dict_item]} for dic in counter_list for dict_item in dic[1]
                     if dic[1][dict_item] > 1 and dict_item == list(parallel.keys())[0]]

    # fixme: to_delete isn't working
    to_delete_obj = [
        f"{dict_item}-{counter_list.index(dic)}" for dic in counter_list for dict_item in dic[1]
        if dic[1][dict_item] <= 1 and dict_item == list(parallel.keys())[0]]

    return temp_list_obj, to_delete_obj


def process_blast_results_line(file_line, f_line_num, f_length, memory_dir) -> list | None:
    """
    Processes blast results file. This function is being run for each line in results file.
    The results are written as a counter list in counter_list.txt output file.

    :param file_line:
    :param f_line_num:
    :param f_length:
    :param memory_dir:
    :return: None
    """

    # fixme: this function doesn't read Counter List.txt correctly
    # FOR EACH LINE IN RESULTS FILE: ↓ ↓ ↓
    # print(f"\rProcessing results line: {f_line_num}", end="")
    already_exists = False

    # Fill counter list up
    # fixme: all first 435 lines are don't appear in counter list. Starts erasing from line 6250.
    if f_line_num < 10:
        print(f"[x[0] for x in counter_list: {[item[0] for item in counter_list][0:10]}")
        print(f"file_line: {file_line}")
    if file_line.split()[0] not in [item[0] for item in counter_list]:
        if f_line_num < 10:
            print(f"{file_line.split()[0]} not in counter list")
        already_exists = False
        counter_list.append([str(file_line.split()[0]), Counter()])
        # if "PAU8" not in [x[0] for x in counter_list]:
        #     print(f"PAU8 is not in counter list. Line: {f_line_num}")
    else:
        if f_line_num < 10:
            print(f"{file_line.split()[0]} already in counter list")
        already_exists = True
    # for item in counter_list:
    #     if file_line.split()[0] in item:
    #         break
    # else:
    #     counter_list.append([str(file_line.split()[0]), Counter()])


    # Add possible matches
    for item in counter_list:
        if item[0] == file_line.split()[0]:
            item[1][file_line.split()[1]] += float(file_line.split()[11])

    # When done processing file
    if f_line_num + 1 == f_length:
        print("Finished file")
        # print(f"\nOriginal counter List: {len(counter_list)}")

        # Save to file
        with open(memory_dir + "Counter List.txt", 'w') as counter_list_f:
            counter_list_f.write(str(counter_list).removeprefix('[').removesuffix(']')
                                 .replace("', ", "' ").replace("})], ", "\n")
                                 .replace('[', "").replace("Counter({", "")
                                 .replace("'", "").replace(':', '').replace("})]", ""))
        print(f"\n\ncounter list (func): {[x[0] for x in counter_list][0:15]}\n\n")

        # return counter_list

    if mp.current_process().name.split('-')[-1] == '1':
        print(f"\rLines processed: {f_line_num} ({int(100 * ((f_line_num + 1) / f_length))}%),"
              f" genes matched: {len(counter_list)}",
              end="")
    # return counter_list[f_line_num]

    for item in counter_list:
        if item[0] == file_line.split()[0] and not already_exists:
            return item


# def blast_results_to_dict(input_dir="Memory\\blast\\bin\\", max_threads=os.cpu_count()):
#     # proc = mp.Process()
#     # Matching
#     results_file = open(input_dir + "results.txt", 'r')
#     results_txt = results_file.readlines()
#
#     lines = []
#     counter_list = []
#     procs = []
#
#     def match_unrecognized_seqs(line, line_num):
#         for item in counter_list:
#             if item[0] == line.split()[0]:
#                 break
#         else:
#             counter_list.append([str(line.split()[0]), Counter()])
#             # print(str(line.split()[0] + " does not exist. Adding it to the list"))
#
#         for item in counter_list:
#             if item[0] == line.split()[0]:
#                 item[1][line.split()[1]] += float(line.split()[11])
#
#             print("\rLines processed: " + str(line_num) + " (" + str(int(100 * ((line_num + 1) / len(results_txt)))) +
#                   "%), genes matched: " + str(len(counter_list)), end="")
#     # return lines
#
#     for line_num, line in enumerate(results_txt):
#         lines.append((line, line_num))
#     print("For loop 1 done")
#
#     # for i in range(max_threads):
#         # procs.append(mp.Process(target=match_unrecognized_seqs, args=(lines[i])))
#     print("For loop 2 done")
#     return procs, counter_list


def make_gene_dict(counter_list):
    """ Counter List:
        :List = [[<query gene 1>, <Counter>],[<query gene 2>, <Counter>],[<query gene 3>, <Counter>]...]
        :Counter = [<target gene 1>, <target gene 2>, <target gene 3>...]

        :counter_list = [[<query gene 1>, [<target gene 1>, <target gene 2>, <target gene 3>...]],
                         [<query gene 1>, [<target gene 1>, <target gene 2>, <target gene 3>...]],
                         [<query gene 1>, [<target gene 1>, <target gene 2>, <target gene 3>...]]]
        :return: dict()
    """

    # input_dir = "C:\\Users\\Yoni\\PycharmProjects\\pythonProject\\Memory\\blast\\bin\\"

    print(f"\n\n\nCounter list from make_gene_dict: {len(counter_list)}")

    # Matching
    # results_file = open(input_dir + "results.txt", 'r')
    # results_txt = results_file.readlines()

    # counter_list = list()
    gene_dict = dict()
    # proc = mp.Process()

    # for line in results_txt:

    # procs[0].start()
    # print("started process")

    for item in [x for x in counter_list if len(x[1]) > 0]:
        gene_dict[item[0]] = item[1].most_common(1)[0][0]
    logging.debug("\nGene dictionary: " + str(gene_dict))

    # results_file.close()
    return gene_dict


# "GeneIdentifier" Definitions
# def identify_genes(input_dir=r"MainInput\\", output_dir=r"MainOutput\\", temp_dir=r"Memory\\DefaultTemp\\",
#                    memory_dir=r"Memory\\DefaultTemp\\GeneLib\\", threshold=0.7,
#                    algorithm_type="Pairwise2-ScoreOnly", seq_len_max_difference=15,
#                    temp_list_filename="percentages.txt", temp_list_format="Normal"):
#
#     genes_list = list()
#
#     """
#     - results_str: Results of matching between sequence in check and every other sequence in
#     memory. From that string the program decides later which of the matches is the best one and assigns it.
#     - final_results_str: Best matches list of every sequence inputted.
#     """
#
#     def open_gene_files(dir):
#         # Open all genes that are stored in memory and prepare them for matching
#         for paths, dirs, files in os.walk(dir):
#             for file in files:
#                 if file.split('.')[1].lower() == "fasta":
#                     print("Opening file: " + file)
#                     gene_seq = open(dir + file, 'r')
#                     genes_list.append(SeqRecord(Seq(gene_seq.read()), name=file))
#                     # genes_dict[Seq(gene_seq.read())] = file
#                     gene_seq.close()
#
#     def apply_sbm_algorithm(input_gene_seq, memory_gene_seq, memory_gene_name, results_str_=""):
#         # Match given sequence (input_gene_seq) to a sequence from memory (memory_gene_sequence)
#         # using a Sequence Based Method.
#
#         """
#         Algorithm Types (algorithm_type):
#         - Hamming: Hamming Distance
#         - Pairwise2: Pairwise2 algorithm from Biopython library.
#         - Pairwise2-ScoreOnly: Pairwise2 algorithm from Biopython library. score_only set to True (shorter runtime).
#         """
#
#         if algorithm_type == "Hamming":
#             mutations = 0
#
#             length = min(len(input_gene_seq), len(memory_gene_seq))
#
#             for i in range(0, length):
#                 char1 = input_gene_seq[i]
#                 char2 = memory_gene_seq[i]
#
#                 if not char2.lower() == char1.lower():
#                     mutations += 1
#
#             hum_dist = mutations / length
#
#             # Save Results in Temp Memory
#             results_str_ += memory_gene_name.split('\\')[-1] + " " + str(mutations) + "/" + str(length) + " " + str(hum_dist) + "\n"
#
#         elif algorithm_type == "Pairwise2":
#             seq1 = Seq(input_gene_seq)
#             seq2 = Seq(memory_gene_seq)
#
#             matches = pairwise2.align.globalxx(seq1, seq2, score_only=False)
#             logging.debug("Match results for gene " + memory_gene_name + ": " +
#                           str(matches[0][2] / (matches[0][4] - matches[0][3])))
#
#             # Save Results in Temp Memory
#             if temp_list_format == "Normal":
#                 results_str_ += memory_gene_name.split('\\')[-1] + " " + str(matches[0][2]) + "/" + \
#                                 str(matches[0][4] - matches[0][3]) + " " + \
#                                 str(matches[0][2] / (matches[0][4] - matches[0][3])) + "\n"
#             elif temp_list_format == "PercentageOnly":
#                 results_str_ += str(matches[0][2] / (matches[0][4] - matches[0][3])) + "\n"
#
#         elif algorithm_type == "Pairwise2-ScoreOnly":
#             seq1 = Seq(input_gene_seq)
#             seq2 = Seq(memory_gene_seq)
#             bigger_seq = max(len(seq1), len(seq2))
#
#             # Apply algorithm
#             matches = pairwise2.align.globalxx(seq1, seq2, score_only=True)
#             logging.debug("Match results for gene " + memory_gene_name + ": " +
#                           str(matches / bigger_seq))
#
#             # Save Results in Temp Memory
#             if temp_list_format == "Normal":
#                 results_str_ += memory_gene_name.split('\\')[-1] + " " + str(matches) + "/" + \
#                                 str(bigger_seq) + " " + str(matches / bigger_seq) + "\n"
#             elif temp_list_format == "PercentageOnly":
#                 results_str_ += str(matches / bigger_seq) + "\n"
#
#         return results_str_
#
#     def check_gene_in_memory(input_gene_file__):
#         # Match given gene (input_gene_file__) to all genes stored in memory.
#
#         results_str__ = ""
#
#         if not os.path.exists(temp_dir):
#             os.makedirs(temp_dir)
#
#         logging.debug("\n\n-----Checking Genes-----\n")
#         input_gene_seq = input_gene_file__.read()
#         logging.debug("From input file: " + input_gene_file__.name)
#
#         for i in range(0, len(genes_list)):
#             if len(input_gene_seq) - seq_len_max_difference <= len(genes_list[i]) \
#                <= len(input_gene_seq) + seq_len_max_difference:
#                 memory_gene_seq = genes_list[i]
#                 results_str__ = apply_sbm_algorithm(
#                     input_gene_seq, memory_gene_seq.seq, memory_gene_seq.name, results_str_=results_str__
#                 )  # genes_dict[genes_list[i]]
#
#                 # print("\n" + str(len(genes_list) - i) + " Files remaining...")
#
#                 # if len(genes_list) > 0:
#                 #     print("\rProcessing " + input_gene_file__.name.split(r"\\")[-1] + ": " +
#                 #           str(100 - int(((len(genes_list) - i) / len(genes_list)) * 100)) + "%", end="")
#
#         return results_str__
#
#     def determine_best_match(seq_in_check, temp_list_txt, final_results_str_=""):
#         # Decide which (if any) of the sequences stored in memory is the same one as the given sequence (seq_in_check).
#
#         logging.debug("\n\n-----Determining Best Match-----\n")
#         top_match_gene = "None"
#         gene_dist_equivalents = {}
#         # _temp_list = open_1_file(temp_dir + temp_list_filename, 'r')
#         # temp_list_txt = _temp_list.readlines()
#         lowest_num = 1.0  # For Hamming Distance etc., where the lower the better
#         highest_num = 0.0  # For pairwise2 alignment etc., where the higher the better
#
#         # Find Lowest/Highest Distance
#         # Hamming Distance (Lowest number)
#         if algorithm_type == "Hamming":
#             # print("Hamming")
#             for line in temp_list_txt:
#                 if line.__contains__(" "):
#                     num = float(line.split()[2])
#                     if num < lowest_num:
#                         lowest_num = num
#                         top_match_gene = str(line.split('.')[0])
#                         logging.debug("Record broken with " + str(num) + " (" + top_match_gene + ")")
#
#         # Pairwise2 (Highest number)
#         elif algorithm_type == "Pairwise2":
#             # print("Pairwise2")
#             for line in temp_list_txt:
#                 if line.__contains__(" "):
#                     if temp_list_format == "Normal":
#                         num = float(line.split()[2])
#                     elif temp_list_format == "PercentageOnly":
#                         num = float(line)
#
#                     if num > highest_num:
#                         highest_num = num
#                         top_match_gene = str(line.split('.')[0])
#                         logging.debug("Record broken with " + str(num) + " (" + top_match_gene + ")")
#
#         # Pairwise2 Score Only (Highest number)
#         elif algorithm_type == "Pairwise2-ScoreOnly":
#             # print("Pairwise2-ScoreOnly")
#             for line in temp_list_txt:
#                 if line.__contains__(" "):
#                     if temp_list_format == "Normal":
#                         num = float(line.split()[2])
#                     elif temp_list_format == "PercentageOnly":
#                         num = float(line)
#                     else:
#                         logging.error("List format '" + temp_list_format + "' is not a valid format."
#                                                                            " \nSetting format to 'Normal'")
#                         print("List format '" + temp_list_format + "' is not a valid format."
#                                                                    " \nSetting format to 'Normal'")
#                         num = float(line.split()[2])
#
#                     if num > highest_num:
#                         highest_num = num
#                         top_match_gene = str(line.split('.')[0])
#                         logging.debug("Record broken with " + str(num) + " (" + top_match_gene + ")")
#                     else:
#                         logging.debug("num: " + str(num) + " < highest_num: " + str(highest_num))
#
#         else:
#             print("Else")
#
#         # Check For Equivalents
#         for line in temp_list_txt:
#             if line.__contains__(" "):
#                 if temp_list_format == "Normal":
#                     num = float(line.split()[2])
#                 elif temp_list_format == "PercentageOnly":
#                     num = float(line)
#
#                 gene_name = str(line.split('.')[0])
#
#                 if algorithm_type == "Hamming" and num == lowest_num or \
#                         algorithm_type == "Pairwise2" and num == highest_num or \
#                         algorithm_type == "Pairwise2-ScoreOnly" and num == highest_num:
#                     gene_dist_equivalents[gene_name] = num
#
#         # Write Output (Results)
#         if len(gene_dist_equivalents) > 0:
#             final_results_str_ += seq_in_check + " = "
#
#             # for i in equal_dist_genes:
#             final_results_str_ += str(gene_dist_equivalents).removeprefix("{").removesuffix("}").replace("'", "") + "\n"
#             return final_results_str_, gene_dist_equivalents
#
#         else:
#             # final_results_str_ += seq_in_check + " = " + top_match_gene + " (" + str(lowest_num) + ")\n"
#             final_results_str_ += seq_in_check + " = None : 1.0\n"
#             return final_results_str_, [top_match_gene]
#
#     def apply_threshold(threshold=0.7):
#         print("\nApplying Threshold...\n")
#         logging.debug("Applying Threshold...")
#         final_results_above_threshold = list()
#
#         for line in final_results_str.splitlines():
#             match_percentage = float(line.split(':')[-1])
#             if match_percentage >= threshold and line.split()[2] != "None":
#                 final_results_above_threshold.append(line)
#
#         results_above_threshold_file = open_1_file(output_dir + "Results.txt", 'w')
#         for i in range(len(final_results_above_threshold)):
#             results_above_threshold_file.write(str(final_results_above_threshold[i]))
#
#     # -----Function Codes End Here-----
#
#     final_results_str = ""
#
#     open_gene_files(memory_dir)
#     current_file = 1
#
#     for paths, dirs, files in os.walk(input_dir):
#         for filename in files:
#             # Open Files
#             temp_list_str = ""
#             # temp_list = open_1_file(temp_dir + temp_list_filename, 'w')
#             # results = open_1_file(output_dir + "Results " + filename + ".txt", 'w')
#
#             logging.debug("Processing " + filename)
#             # start = time.time()
#             input_gene_file_ = open_1_file(input_dir + filename, 'r')
#             results_str = check_gene_in_memory(input_gene_file_)
#             # end = time.time()
#             # print("\nTime: " + str(end - start))
#
#             # Estimate time to finish operation
#             estimated_time = 0
#             if "%.1f" % (current_file / len(files) * 100) == 0.1:
#                 start = time.time()
#             elif "%.1f" % (current_file / len(files) * 100) == 0.2:
#                 end = time.time()
#                 estimated_time = (end - start) * 1000
#                 estimated_time = time.gmtime(estimated_time)
#                 estimated_time = time.strftime("%H Hours, %M Minutes, %S Seconds", estimated_time)
#
#             if estimated_time == 0:
#                 print("\r" + str(len(files) - current_file) + " Files remaining... (" +
#                       str("%.1f" % (current_file / len(files) * 100)) + "%, Currently processing: " +
#                       str(filename) + ")", end="")
#
#             else:
#                 print("\r" + str(len(files) - current_file) + " Files remaining... (" +
#                       str("%.1f" % (current_file / len(files) * 100)) + "%, Currently processing: " +
#                       str(filename) + ") Estimated time: " + str(estimated_time), end="")
#
#             logging.debug("results_str before writing = " + results_str)
#             # temp_list.write(results_str)
#             # temp_list.close()
#             temp_list_str += results_str + "\n"
#             # logging.debug("temp_list file includes: " + open(temp_dir + temp_list_filename, 'r').read())
#             results_str = ""
#             logging.debug("Results String has been transported to file")
#
#             final_results_str += determine_best_match(filename, temp_list_str.splitlines())[0]
#             # print("final_results_str = " + final_results_str)
#             close_all_files()
#             logging.debug(filename + " has been successfully processed")
#
#             current_file += 1
#
#     results_file = open_1_file(output_dir + "Results Before Threshold.txt", 'w')
#     results_file.write(final_results_str)
#     apply_threshold(threshold=threshold)
#
#     close_all_files()  # Just in case


# def open_1_file(path, mode):
#     file = open(path, mode)
#     active_files.append(file)
#     return file


# def close_all_files():
#     logging.debug("\n\n-----Closing Files-----\n")
#     for file in active_files:
#         file.close()
#         active_files.remove(file)


# "CleanFolder" Definitions
def delete_files(_dir):
    for path, dirs, files in os.walk(_dir):
        for folder in dirs:
            print("Removing Folder " + folder)
            rmtree(_dir + folder)
        for file in files:
            print("Removing File " + file)
            os.remove(_dir + file)
    print("All files successfully deleted in " + _dir + "\n")


def delete_old_logs(log_dir, limit=10):
    def date(filename):
        return datetime.fromtimestamp(os.stat(log_dir + filename).st_mtime)

    files_list = list()
    do_print = False

    for path, dirs, files in os.walk(log_dir):
        for file_ in files:
            files_list.append(file_)

    files_list.sort(key=date, reverse=True)

    for i in range(len(files_list)):
        if i >= limit:
            do_print = True
            os.remove(log_dir + str(files_list[i]))

    if do_print:
        print("Deleting old log files\n")


def show_and_get_files(file_num, instruction, input_dir_=const_input_dir, acceptable_formats=()):
    check_acceptable_formats = True
    files_dict = dict()

    if not acceptable_formats:
        check_acceptable_formats = False

    print(instruction)
    current_gene_num = 1
    for paths, dirs, files in os.walk(input_dir_):
        for file in files:
            if not check_acceptable_formats or\
                    check_acceptable_formats and acceptable_formats.__contains__(file.split('.')[-1].lower()):
                print(str(current_gene_num) + ". " + file)
                files_dict[str(current_gene_num)] = file
                current_gene_num += 1

    # Get input
    is_file = False
    while not is_file:
        print("\nFile " + str(file_num) + ":")
        filename = input()
        if os.path.exists(input_dir_ + filename):
            print("File 1 is " + filename)
            is_file = True

        elif files_dict.__contains__(filename):
            if os.path.exists(input_dir_ + files_dict[filename]):
                print("File is valid!")
                filename = files_dict[filename]
                return filename

        else:
            print("File not found in directory\n")
            is_file = False
    return filename


def check_similarity(format="genelist2", input_dir=anlys_inp_dir, output_dir=anlys_op_dir):

    if format == "genelist1":
        # Enter files
        filename_1 = show_and_get_files(1, "Choose 2 files you wish to check:\n", input_dir_=anlys_inp_dir, acceptable_formats=[format])
        filename_2 = show_and_get_files(2, "Choose 2 files you wish to check:\n", input_dir_=anlys_inp_dir, acceptable_formats=[format])

        # Apply Search
        file_1 = open(input_dir + str(filename_1), 'r')
        file_2 = open(input_dir + str(filename_2), 'r')
        text_1 = file_1.read().splitlines()
        text_2 = file_2.read()
        file_1.close()
        file_2.close()
        matches = []

        for line in text_1:
            occurrence = text_2.find(line.split()[0])
            # print(text_2[occurrence:text_2.find()])

            if occurrence > -1:
                print("Occurrence found: " + line + "\n")
                matches.append(line)

        percentage = (len(matches) * 2) / (len(text_1) + len(text_2.splitlines())) * 100

    elif format == "genelist2":
        # Open & read file
        filename_1 = show_and_get_files(1, "Choose 2 files you wish to check:\n", input_dir_=anlys_inp_dir, acceptable_formats=[format])
        file_1 = open(input_dir + str(filename_1), 'r')
        text_1 = file_1.read().splitlines()
        file_1.close()

        filename_1 = f"{text_1[0].split(':')[1].removesuffix(', ')}.genelist1"
        filename_2 = f"{text_1[0].split(':')[2]}.genelist1"

        text_1 = [line for line in text_1 if not line.startswith("1") and not line.startswith("#")]
        matches = []

        # Search
        for line in [l for l in text_1
                     if not "*" in l and
                     not "|" in l]:
            matches.append(line)

        non_asters = [match for match in text_1 if "*" not in match]
        percentage = (len(matches) * 2) / ((len(text_1) - len(non_asters)) + (2 * len(non_asters))) * 100

    else:
        raise ValueError("Incompatible input format. format needs to be 'genelist1' or 'genelist2'")

    # Show results
    percentage = int(percentage)
    print(matches)
    print("\n" + str(len(matches)) + " matches where found (" + str(percentage) + "%)!\n")

    # Output Matching Genes List
    results_file = open(output_dir + " " +
                        filename_1.split('.')[0].capitalize() + " & " +
                        filename_2.split('.')[0].capitalize() + " " +
                        "Match Results.genelist", 'w')
    for match in matches:
        results_file.write(match + "\n")
    results_file.close()

    # Diagram
    if format == "genelist1":
        diagram = venn2(subsets=(
            (len(text_1) - len(matches)),
            len(text_2.splitlines()) - len(matches),
            len(matches)),
            set_labels=(filename_1, filename_2))

    elif format == "genelist2":
        unmatched_1 = [line for line in text_1 if line.split()[2] == "X"]
        unmatched_2 = [line for line in text_1 if line.split()[1] == "X"]

        diagram = venn2(subsets=(
            len(unmatched_1),
            len(unmatched_2),
            len(matches)),
            set_labels=(filename_1, filename_2))

        print(f"Unmatched 1: {len(unmatched_1)}")
        print(f"Unmatched 2: {len(unmatched_2)}")
        print(f"Matched: {len(matches)}")

    # print("len(text_1): " + str(len(text_1)))
    # print("len(matches): " + str(len(matches)))
    # print("len(text_2): " + str(len(text_2.splitlines())))

    plt.title("Similarity Check Results\n" + str(percentage) + "% of lines match")
    # plt.annotate(str(percentage) + "% of lines match", xy=(0,0))
    plt.show()

    # Close files
    # file_1.close()
    # file_2.close()


def check_distribution(input_dir=anlys_inp_dir):
    in_array = []

    # Add genes to in_array
    i = 1
    for paths, dirs, files in os.walk(input_dir):
        for file in files:
            print("File: ", file)
            print("Number of files: ", i)
            gene_file = open(input_dir + file, 'r')
            in_array.append(len(gene_file.read().replace("\n", "").replace("\r", "")))
            gene_file.close()
            i += 1

    # Show graph
    print("Array contains " + str(len(in_array)) + " items")
    plt.hist(in_array, bins=len(in_array))
    plt.ylabel('Probability')
    plt.show()


def check_asterisks_percentage(input_dir=anlys_inp_dir):
    inp_f = open(
        input_dir +
        show_and_get_files(1, "Choose GENELIST2 file:\n", acceptable_formats=["genelist2"], input_dir_=input_dir),
        'r')
    inp_txt = [x for x in inp_f.readlines() if len(x.split()) == 3]
    aster_num = 0

    for line in inp_txt:
        if line.split()[1] == "X" or line.split()[2] == "X":
            aster_num += 1

    genes_num = ((len(inp_txt) - aster_num) * 2) + aster_num
    percentage = (aster_num / genes_num) * 100
    hist_total, hist_aster = [], []

    print(f"Asterisks: {aster_num}; Total: {genes_num}")

    for i in range(genes_num):
        hist_total.append("")

    for i in range(aster_num):
        hist_aster.append("")

    # Show Diagram
    plt.hist(hist_total, bins="auto", edgecolor='black', color="blue", alpha=0.5)
    plt.hist(hist_aster, bins="auto", edgecolor='black', color="yellow", alpha=0.5)

    plt.title("Percentage of Unmatched Genes")
    plt.ylabel("Number of Genes")
    plt.text(-0.25, 2 * 1000, "Unmatched: " + "%.2f" % percentage + "%")
    plt.show()
