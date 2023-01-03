from libs_and_dirs import *


""" Redundant: Dictionaries Structure:
    familiar_genes[<Gene Name>] = [<Location in Genome 1>, <Location in Genome 2>]
    unfamiliar_genes[<Gene Name>] = [<Genome Number (1 / 2)>, <Location in Genome>]
"""


# Analysis Functions
def apply_algorithm(inp_text, k):
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

    :param inp_text: Double genelist as a string.
    :return: None
    """

    # Adjust input text
    if type(inp_text) == str:
        inp_text = inp_text.splitlines()

    elif '\n' in next(line for line in inp_text):
        inp_text = "".join(inp_text).splitlines()

    # logging.debug(f'Text file: "{filename}"')

    species = next((line for line in inp_text if line.lower().startswith("##species")))
    species = species.removeprefix(species.split(':')[0] + ':').split('-')
    species_1_name = species[0]
    species_2_name = species[1]

    # Remove comments
    inp_text = tuple(line for line in inp_text if not line.startswith('#'))

    # region Add fillers to missing genes
    pos_1s = tuple(int(line.split()[1]) for line in inp_text if '*' not in line.split()[1])
    pos_2s = tuple(int(line.split()[2]) for line in inp_text if '*' not in line.split()[2])

    len_genome_1 = sorted(pos_1s)[-1]
    len_genome_2 = sorted(pos_2s)[-1]

    lines_to_add_1 = tuple(f"FILLER {pos} *" for pos in range(len_genome_1)
                           if pos not in pos_1s and pos != 0)
    lines_to_add_2 = tuple(f"FILLER * {pos}" for pos in range(len_genome_2)
                           if pos not in pos_2s and pos != 0)

    lines_to_add = lines_to_add_1 + lines_to_add_2
    total_num_of_genes = (None, len_genome_1, len_genome_2)

    print(f"\n\nFillers: {lines_to_add}")
    inp_text = inp_text + lines_to_add
    # endregion

    def make_nbrhood(center_gene: Allele, nbrhood_num: int, len_genome: tuple) -> tuple:
        """
        Create neighborhood for given gene.

        :param center_gene: The gene around which the neighborhood is built. center_gene = [<gene name>, <location in neighborhood 1>, <location in neighborhood 2>].
        :param nbrhood_num: 1 or 2. Number of neighborhood. 1: Left in genelist, 2: Right in genelist.
        :param len_genome: Tuple of both genome lengths.
        :return: Neighborhood.
        """

        center_gene_num = center_gene.locations[nbrhood_num]   # Just for readability
        logging.debug(f"Center Gene: {center_gene.name} ({center_gene_num})\n\n")

        # Set neighborhood
        nbrhood = [
            line.split()[0] for line in inp_text if

            # 1. If gene is homozygous (exists in both genomes)
            line.split()[nbrhood_num] != no_gene_symbol and

            # 2. If gene is either a:
            # - Direct Neighbor (e.g 2 and 3 in [1,2,3,4,5])
            (center_gene_num - k <= int(line.split()[nbrhood_num]) <= center_gene_num + k or
             # - Cyclical Neighbor (e.g 1 and 5 in [1,2,3,4,5])
             center_gene_num + k - len_genome[nbrhood_num] >= int(line.split()[nbrhood_num]) or
             center_gene_num - k + len_genome[nbrhood_num] <= int(line.split()[nbrhood_num])) and

            # 3. Gene is not central gene (main gene is not supposed to be in neighborhood)
            line.split()[0] != center_gene.name
        ]

        tuple(
            nbrhood.remove(item)
            for item in nbrhood
            if nbrhood.count(item) > 1 and item != "FILLER"
        )

        nbrhood = tuple(nbrhood)

        # If there was an error
        if len(nbrhood) != k * 2:
            print(f"\nBroken neighborhood for gene {center_gene.name}; {center_gene_num} in genome {nbrhood_num}:\n"
                  f"{nbrhood}")

        return nbrhood

    # Apply SI to each gene
    text_no_asterisks = tuple(line for line in inp_text if no_gene_symbol not in line)
    try:
        synteny_indexes_of_genes = []

        for line_num, current_gene_line in enumerate(text_no_asterisks):
            if line_num % 10 == 0:
                print(f"\r{int(100 * ((line_num + 1) / len(text_no_asterisks)))}%", end="")

            # Define current gene
            gene_in_check = Allele(current_gene_line.split()[0],  # Name
                                   (current_gene_line.split()[1],  # Location in genome 1
                                   current_gene_line.split()[2]))  # Location in genome 2

            # Create neighborhoods for current gene
            neighborhoods = (None,
                             make_nbrhood(gene_in_check, 1, total_num_of_genes),
                             make_nbrhood(gene_in_check, 2, total_num_of_genes))

            # Find common genes
            intersection = tuple(gene_1

                                 for gene_1, gene_2 in zip(neighborhoods[1], neighborhoods[2], strict=True)

                                 if gene_1 in neighborhoods[2] and no_gene_symbol not in gene_1
                                 or
                                 no_gene_symbol in gene_1 and no_gene_symbol in gene_2
                                 )

            # logging.debug("Intersection: " + str(intersection))

            # SI formula
            x = (1 / (k * 2)) * len(intersection)
            synteny_indexes_of_genes.append(x)
            logging.debug(f"\nSI for gene {gene_in_check.name}: {x}\n")

    except ValueError as exception:
        try:
            raise Exception(f"Neighborhoods are of different lengths ({len(neighborhoods[1])}, {len(neighborhoods[2])})")
        except:
            raise Exception(exception)

    # Average out gene SI's to find genome's SI
    if len(synteny_indexes_of_genes) > 0:
        synteny_index_of_genome = sum(synteny_indexes_of_genes) / len(synteny_indexes_of_genes)
    else:
        synteny_index_of_genome = None

    return (species_1_name, species_2_name), synteny_index_of_genome


def find_genome(dir):
    """
    Find GFF3 and Fasta in specified directory

    :param dir: Genome directory.
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

        if gff3 and fasta:
            break

    return gff3, fasta


def separate_by_gene(genelist: (), sequence: str | SeqRecord):
    """
    Separate genomic sequence by gene.

    :param genelist: Single-gene-list.
    :param sequence: Whole genomic sequence.
    :return: Str in a FASTA format.
    """

    sequence = str(sequence)
    modified_seq = ""

    if '>' in sequence:
        sequence = "".join([line for line in sequence.splitlines() if '>' not in line])

    for gene in genelist:
        modified_seq += f">{gene.name}\n{sequence[gene.start_pos-1:gene.end_pos]}\n"

    return modified_seq


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


def blast(is_protein=False, evalue=1e-50, blast_dir_=blast_dir) -> None:
    """
    Run BLAST from on given genomes.
    \nis_protein: True = Match protein sequences. False = Match nucleotide sequences.
    """

    if not os.path.exists(blast_dir_):
        os.makedirs(blast_dir_)

    path = os.path.dirname(os.path.abspath(__file__)) + '\\' + blast_dir_

    # Align (CMD)
    try:
        if is_protein:
            print(subproc_run(("blastp",
                               "-subject", "Genome1.fasta",
                               "-query", "Genome2.fasta",
                               "-out", "results.txt",
                               "-outfmt", "6",
                               "-evalue", str(evalue),
                               "-max_target_seqs", "5",
                               "-task", "megablast"),
                              text=True, shell=True, cwd=path))
        else:
            print(subproc_run(("blastn",
                               "-subject", "Genome1.fasta",
                               "-query", "Genome2.fasta",
                               "-out", "results.txt",
                               "-outfmt", "6",
                               "-evalue", str(evalue),
                               "-max_target_seqs", "5",
                               "-task", "megablast"),
                              text=True, shell=True, cwd=path))

    except Exception as exception:
        raise Exception(f"There has been a problem with BLAST: {exception}")

    # todo: Get read of shell=True
    # todo: Add to blast alignment parameters the following:
    """ - Open gap
        - Gap penalty
        - nucl match reward
        - nucl mismatch penalty
        - megablast V
        - limit output V
    """


def get_best_match(candidates: iter) -> str:
    """
    Sorts a given group of matches by their bitscore.

    :param candidates: A parallel candidate
    :return: Best match
    """
    edited_group: list = list(candidates)
    edited_group.sort(key=lambda match_: candidates[match_], reverse=True)
    return edited_group[0]


def get_parallels_to_remove(parallel_real_matches, parallels, i: int, blast_match: PossibleMatch, candidate):
    if mp.current_process().name.split('-')[-1] == '1' and i % 25 == 0:
        pass
    print(f"\r{100 * (i / len(parallel_real_matches))}%", end="")

    # corres_paral_couple = next((paral for paral in classed_parallels
    #                            if paral.gene_2 == candidate), None)


    # if corres_paral_couple and blast_match.main_gene != corres_paral_couple.gene_1:
    if any(parallel.main_gene == blast_match.main_gene and candidate in parallel.candidates
           for parallel in parallels) and \
       DefiniteMatch(blast_match.main_gene, candidate) not in parallel_real_matches:
        return i, candidate
    # if match.main_gene != getattr(next((definite_match for definite_match in classed_parallels if definite_match.gene_2 == candidate), None), "gene_1", None):
        # <Corresponding parallel group>.<best match of parallel>


def process_blast_results_line(file_line) -> PossibleMatch:
    """
    Processes blast results file. This function is being run for each line in results file.
    The results are written as a counter list in counter_list.txt output file.

    :param file_line:
    :return: Result line as a PossibleMatch object
    """

    # FOR EACH LINE IN RESULTS FILE: ↓ ↓ ↓
    return PossibleMatch(file_line.split()[0], Counter({file_line.split()[1]: int(file_line.split()[11])}))


def create_parallel(match, candidate):
    """
    Runs this function for every parallel-to-be.

    :return:
    """

    return PossibleMatch(match.main_gene, Counter({candidate: match.candidates[candidate]}))


# DON'T DELETE- for folder cleaning function
def delete_files(_dir):
    for path, dirs, files in os.walk(_dir):
        for folder in dirs:
            rmtree(_dir + folder)
        for file in files:
            os.remove(_dir + file)


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
        unmatched_1 = [line for line in text_1 if line.split()[2] == no_gene_symbol]
        unmatched_2 = [line for line in text_1 if line.split()[1] == no_gene_symbol]

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
        if line.split()[1] == no_gene_symbol or line.split()[2] == no_gene_symbol:
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
