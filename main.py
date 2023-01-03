from libs_and_dirs import *
import func


def main_menu():
    print(
        "Choose an action to apply on the files in the input folder:\n"
        "0. Close program\n"
        "1. Multiples\n"
        "2. Convert\n"
        "3. Analyze\n"
        "4. Text Editor\n"
        "5. Clean Folder"
    )
    command = input()
    if command.isnumeric():
        command = int(command)

        if command == 0:
            close_program()
        elif command == 1:
            multiple_menu()
        elif command == 2:
            convert_menu()
        elif command == 3:
            analyze_menu()
        else:
            print("Invalid command!")
            main_menu()
    else:
        print("Invalid command!")
        main_menu()


def convert_menu():
    print(
        "Convert:\n"
        "0. Return to main menu\n"
        "\nGFF3 & FASTA:\n\n"
        "1. -> GENELIST2\n"
        "2. Compare a group of genomes\n"
        "\n GENELIST1:\n\n"
        "3. -> GENELIST2\n"
        "4. -> GENELIST1 numbered according to location in list"
    )
    command = input()
    if command.isnumeric():
        command = int(command)

        if command == 0:
            main_menu()
        elif command == 1:
            start = time.perf_counter()
            try:
                # Validate input files
                chosen_gff3_name_1 = func.find_genome(file_1_dir)[0]
                chosen_fasta_name_1 = func.find_genome(file_1_dir)[1]
                chosen_gff3_name_2 = func.find_genome(file_2_dir)[0]
                chosen_fasta_name_2 = func.find_genome(file_2_dir)[1]

                print(f"Species 1 GFF3: {func.find_genome(file_1_dir)[0]}\n"
                      f"Species 1 Fasta: {func.find_genome(file_1_dir)[1]}\n"
                      f"Species 2 GFF3: {func.find_genome(file_2_dir)[0]}\n"
                      f"Species 2 Fasta:  {func.find_genome(file_2_dir)[1]}")
                logging.debug(f"Species 1 GFF3: {func.find_genome(file_1_dir)[0]}\n"
                              f"Species 1 Fasta: {func.find_genome(file_1_dir)[1]}\n"
                              f"Species 2 GFF3: {func.find_genome(file_2_dir)[0]}\n"
                              f"Species 2 Fasta:  {func.find_genome(file_2_dir)[1]}")

                species, outp_genelist = convert_genomes_to_genelist2(
                    gffs=(TextFile(file_1_dir + chosen_gff3_name_1).text,
                          TextFile(file_2_dir + chosen_gff3_name_2).text),
                    fastas=(TextFile(file_1_dir + chosen_fasta_name_1).text,
                            TextFile(file_2_dir + chosen_fasta_name_2).text),
                                             )

                print("\nDouble genelist has been successfully made!"
                      f"\nSaving it as a file in {const_output_dir}...")

                # Save as file
                outp_genelist_txt = "\n".join(outp_genelist)

                with open(f"{const_output_dir}{species[0]}-{species[1]}.genelist2", 'w') as out_file:
                    out_file.write(f"##Species:{species[0]}-{species[1]}\n"
                                   f"##Order:List\n")
                    out_file.write(outp_genelist_txt)

                # Stop timer
                end = time.perf_counter()
                time_passed = end - start
                if time_passed > 60:
                    minutes = floor(time_passed / 60)
                    seconds = time_passed % 60
                    logging.debug(f"\nTook {minutes} minutes and {seconds} seconds to finish operation")
                    print(f"\nTook {minutes} minutes and {seconds} seconds to finish operation")
                else:
                    logging.debug(f"\nTook {time_passed} seconds to finish operation")
                    print(f"\nTook {time_passed} seconds to finish operation")

                main_menu()

            except Exception as ex:
                # Stop timer
                end = time.perf_counter()
                time_passed = end - start
                if time_passed > 60:
                    minutes = floor(time_passed / 60)
                    seconds = time_passed % 60
                    logging.debug(f"\nTook {minutes} minutes and {seconds} seconds to finish operation")
                    print(f"\nTook {minutes} minutes and {seconds} seconds to finish operation")
                else:
                    logging.debug(f"\nTook {time_passed} seconds to finish operation")
                    print(f"\nTook {time_passed} seconds to finish operation")

                raise ex


        elif command == 2:
            print("Applying comparison between:")
            for i, file in tuple(enumerate(files) for (path, dirs, files) in os.walk(const_input_dir)
                              for _, _, files in (path, dirs, files)):
                print(f"{i + 1}. {file}")

            # fixme Change to folders(?) (.gff, .fasta)
            print("Sorting input files... ", end="")
            for _, _, files in os.walk(const_input_dir):
                gffs = tuple(file for file in files if file.lower().split('.')[-1] in gff_formats)
                fastas = tuple(file for file in files if file.lower().split('.')[-1] in fasta_formats)
                genome_files = dict()

                # Get gffs & fastas paired
                print("Done\nPairing gffs & fastas... ", end="")
                for gff_filename in gffs:
                    fasta_file = tuple((filename_fa, fasta_file_.close()) for filename_fa in fastas
                                       if gff_filename in (fasta_file_ := open(const_input_dir + filename_fa, 'r')).read())

                    fasta_file = tuple(file for file in fasta_file if file)[0]
                    genome_files[gff_filename] = fasta_file
                    continue

                # Run for each genome
                print(f"Done\nStarting long run for genomes: {gffs}")
                for gff_name_1 in gffs:
                    gff_1 = TextFile(const_input_dir + gff_name_1)
                    fasta_1 = TextFile(const_input_dir + genome_files[gff_1.name])

                    for gff_name_2 in tuple(gff_name_2 for gff_name_2 in gffs if gff_name_2 != gff_1.name):
                        gff_2 = TextFile(const_input_dir + gff_name_2)
                        fasta_2 = TextFile(const_input_dir + genome_files[gff_2.name])

                        # region Run
                        start = time.perf_counter()
                        try:
                            # Validate input files

                            print(f"Species 1 GFF3: {gff_1.name}\n"
                                  f"Species 1 Fasta: {fasta_1.name}\n"
                                  f"Species 2 GFF3: {gff_2.name}\n"
                                  f"Species 2 Fasta:  {fasta_2.name}")

                            logging.debug(f"Species 1 GFF3: {gff_1.name}\n"
                                          f"Species 1 Fasta: {fasta_1.name}\n"
                                          f"Species 2 GFF3: {gff_2.name}\n"
                                          f"Species 2 Fasta:  {fasta_2.name}")

                            species, outp_genelist = convert_genomes_to_genelist2(
                                gffs=(gff_1.text, gff_2.text),
                                fastas=(fasta_1.text, fasta_2.text),
                            )

                            print("\nDouble genelist has been successfully made!"
                                  f"\nSaving it as a file in {const_output_dir}...")

                            # Save as file
                            outp_genelist_txt = "\n".join(outp_genelist)

                            with open(f"{const_output_dir}{species[0]}-{species[1]}.genelist2", 'w') as out_file:
                                out_file.write(f"##Species: {species[0]}-{species[1]}\n"
                                               f"##Order: List\n")
                                out_file.write(outp_genelist_txt)

                            # Stop timer
                            end = time.perf_counter()
                            time_passed = end - start
                            if time_passed > 60:
                                minutes = floor(time_passed / 60)
                                seconds = time_passed % 60
                                logging.debug(f"\nTook {minutes} minutes and {seconds} seconds to finish operation")
                                print(f"\nTook {minutes} minutes and {seconds} seconds to finish operation")
                            else:
                                logging.debug(f"\nTook {time_passed} seconds to finish operation")
                                print(f"\nTook {time_passed} seconds to finish operation")

                            main_menu()

                        except Exception as ex:
                            # Stop timer
                            end = time.perf_counter()
                            time_passed = end - start
                            if time_passed > 60:
                                minutes = floor(time_passed / 60)
                                seconds = time_passed % 60
                                logging.debug(f"\nTook {minutes} minutes and {seconds} seconds to finish operation")
                                print(f"\nTook {minutes} minutes and {seconds} seconds to finish operation")
                            else:
                                logging.debug(f"\nTook {time_passed} seconds to finish operation")
                                print(f"\nTook {time_passed} seconds to finish operation")

                            raise ex
                        # endregion

        elif command == 3:
            print("This command does not exist yet!")
            convert_menu()
        elif command == 4:
            # Get input
            file = func.show_and_get_files(1, "Choose a genelist(1/2) file:\n",
                                           acceptable_formats=("genelist1", "genelist2"))
            order_seqs_in_genelist(file)
            main_menu()
        else:
            print("Invalid command!")
            convert_menu()
    else:
        print("Invalid command!")
        convert_menu()


def analyze_menu():
    print(
        "Analyze:\n"
        "0. Return to main menu\n"
        "1. Apply SI algorithm\n"
        "2. Check similarity of genomes\n"
        "3. Check files length distribution\n"
        "4. Check percentage of unmatched genes"
    )
    command = input()
    if command.isnumeric():
        command = int(command)

        if command == 0:
            main_menu()
        elif command == 1:
            acceptable_formats = "genelist2",
            results_filename = "SI Results.txt"

            # Get input
            filename = func.show_and_get_files(1,
                                               "Choose a genelist2 file to compare using SI:\n",
                                               input_dir_=anlys_inp_dir, acceptable_formats=acceptable_formats)

            with open(anlys_inp_dir + filename, 'r') as inp_f:
                inp_text = inp_f.read().splitlines()

            while True:
                try:
                    k = int(input("k = ?\n"))
                    break
                except:
                    print("Invalid Input! k should be a natural number")

            species, si = func.apply_algorithm(inp_text, k)
            if si:
                print(f"SI for genomes {species[0]} and {species[1]} = {si}\n")
            else:
                print("\nERROR: No recognized genes found. SI is NONE")

            # Write results in results file
            if not os.path.exists(anlys_op_dir + results_filename):
                results_file = open(anlys_op_dir + results_filename, 'w')
            else:
                results_file = open(anlys_op_dir + results_filename, 'a')
                results_file.write("\n\n")

            results_file.write(f"Comparison results for genomes {species[0]} and {species[1]}: SI = {si}")
            results_file.close()
            main_menu()

        elif command == 2:
            func.check_similarity(format="genelist1")
            main_menu()
        elif command == 3:
            func.check_distribution()
            main_menu()
        elif command == 4:
            func.check_asterisks_percentage()
            main_menu()
        else:
            print("Invalid command!")
            analyze_menu()
    else:
        print("Invalid command!")
        analyze_menu()


def multiple_menu():
    print("Do the following to all files in input at once:\n"
          "0. Return to main menu\n"
          "1. Get double genelist\n"
          "2. Apply SI-algorithm")

    command = input()
    if command.isnumeric():
        command = int(command)

        if command == 0:
            main_menu()
        # Get double genelist
        elif command == 1:
            pass
        # Apply SI
        elif command == 2:
            acceptable_formats = "genelist2",
            # results_filename = "SI Results.txt"
            pool = mp.Pool()

            # region Get input
            results_filename = input("What to call the final results file?\n")

            while True:
                try:
                    k = int(input("k = ?\n"))
                    break
                except:
                    print("Invalid Input! k should be a natural number")

            files_to_analyze = ()
            for _, _, files in os.walk(anlys_inp_dir):
                files_to_analyze = tuple(TextFile(anlys_inp_dir + filename) for filename in files
                                         if filename.split('.')[-1].lower() in acceptable_formats)
            # endregion

            # region Apply algorithm
            genelist_txts = tuple(file.text for file in files_to_analyze)

            results = tuple(pool.imap_unordered(partial(func.apply_algorithm, k=k), genelist_txts, max_threads))
            species_names, genomes_si = tuple(x[0] for x in results), tuple(x[1] for x in results)
            # endregion

            # region Save to file
            with open(anlys_op_dir + results_filename + ".txt", 'w') as op_file:
                op_file.write("SI results:")

                for species_name, result in tuple(zip(species_names, genomes_si)):
                    print(species_name[0], species_name[1], result)
                    op_file.write(f"\n\n{species_name[0]}-{species_name[1]} = {result}")
            # endregion
            main_menu()


def text_editor_menu():
    print("Text Editor:\n"
          "0. Return to main menu\n"
          "1. ")

    command = input()
    if command.isnumeric():
        command = int(command)

        if command == 0:
            main_menu()
        elif command == 1:
            func.merge_txt_files()
            main_menu()
        else:
            print("Invalid command!")
            analyze_menu()
    else:
        print("Invalid command!")
        analyze_menu()


def convert_genomes_to_genelist2(gffs: tuple[str, str], fastas: tuple[str, str]) -> tuple:
    """
    Converts GFF & FASTA texts of both genomes to a single GENELIST2 file ready for SI comparison.

    :return: Double genelist of both genomes as line-splitted str.
    """

    pool = mp.Pool()

    # Build Genomes
    print("\nBuilding genomes... ", end="")
    genome_1 = Genome(annotations=gffs[0], sequence=fastas[0])
    genome_2 = Genome(annotations=gffs[1], sequence=fastas[1])
    print("Done")

    # Create double genelist
    print("Creating a double gene-list... ", end="")
    double_genelist = DoubleGenelist(genome_1, genome_2)
    print("Done\n"
          f"Double gene-list has been created for species: {' & '.join(double_genelist.genomes)}")

    # Check if blast is required
    print("Is BLAST needed? ", end="")
    do_blast = any(None in tuple for tuple in double_genelist.list)
    print(do_blast)

    if do_blast:
        print("\n\033[4m\033[1m" + "Blast:" + "\033[0m")

        # Prepare sequence as FASTA for BLAST
        print("Saving genome as a fasta file... ", end="")
        fasta_by_gene_1 = func.separate_by_gene(genome_1.genelist, genome_1.sequence)
        fasta_by_gene_2 = func.separate_by_gene(genome_2.genelist, genome_2.sequence)

        with open(blast_dir + "\\Genome1.fasta", 'w') as blast_fasta_1:
            blast_fasta_1.write(fasta_by_gene_1)
        with open(blast_dir + "\\Genome2.fasta", 'w') as blast_fasta_2:
            blast_fasta_2.write(fasta_by_gene_2)
        print("Done")

        # region Blast
        print("Blasting...")
        func.blast(blast_dir_=blast_dir, is_protein=False)

        # Analyze BLAST results
        with open(blast_dir + "results.txt", 'r') as results_file:
            results_txt = tuple(results_file.read().splitlines())
        blast_matches: tuple[PossibleMatch] | tuple[tuple] = \
            tuple(pool.map(func.process_blast_results_line, results_txt, max_threads))

        print(f"Blast matches: {tuple((x.main_gene, x.candidates) for x in blast_matches)[:20]}")
        # endregion

        # Adjusting blast results:
        # region Remove parallels
        # Parallel = A possible match that appears more than once
        print("\nLooking for parallels... ", end="")

        # Find parallels
        print("\nStep 1")
        total_blast_cands = tuple(candidate for match in blast_matches for candidate in match.candidates)

        print("Step 2")
        # starmap_args = pool.starmap(func.get_starmap_arguments_parallel, (range(len(total_blast_cands)), blast_matches), max_threads)
        starmap_args = tuple((match, candidate) for match in blast_matches for candidate in match.candidates
                             if total_blast_cands.count(candidate) > 1)

        print("Step 3")
        parallels = tuple(pool.starmap(func.create_parallel, starmap_args, chunksize=max_threads))
        print(f"Done\nFound {len(parallels)} parallels")

        # If there are parallels
        if len(parallels) > 0:
            print(f"Parallels: {[paral.main_gene for paral in parallels][0:20]}")
            print("Grouping parallels... ", end="")

            # Group parallels by the paralleled gene & sort each type of parallel by bitscore
            grouped_parallels = tuple(list(x) for _, x in groupby(parallels, key=lambda paral: paral.main_gene))
            grouped_parallels = tuple({match.main_gene: sum([tuple(match.candidates.values())[0] for group in grouped_parallels if match in group])}
                                      for match in {counter for group in grouped_parallels for counter in group})
            print(f"Done\ngrouped_parallels: {grouped_parallels[:20]}")

            parallel_real_matches: tuple[DefiniteMatch] = tuple(
                DefiniteMatch(tuple(group)[0], func.get_best_match(group))
                for group in grouped_parallels)

            # Remove parallels in a new tuple
            print("Done\n\nRemoving parallels... ", end="")
            edited_blast_matches = deepcopy(blast_matches)
            tuple(
                (edited_blast_matches[i].candidates.pop(candidate))
                for i, match in enumerate(blast_matches) for candidate in match.candidates

                if any(parallel.main_gene == match.main_gene and candidate in parallel.candidates
                       for parallel in parallels) and \
                DefiniteMatch(match.main_gene, candidate) not in parallel_real_matches
            )

            # todo: Multiprocessing for faster runtime
            # starmap_args = tuple((i, match, candidate)
            #                      for i, match in enumerate(blast_matches) for candidate in match.candidates)
            #
            # print("Args list is finished")
            #
            # to_pop = tuple(pool.starmap(
            #     partial(func.get_parallels_to_remove, parallels_real_matches=set(parallel_real_matches), parallels=set(parallels)),
            #     starmap_args, chunksize=max_threads))
            # print("Multiprocessing is finished")
            # to_pop = tuple(item for item in to_pop if item)
            # print(to_pop)
            #
            # tuple(temp_blast_matches[tuple[0]].candidates.pop(tuple[1]) for tuple in to_pop)

            blast_matches = edited_blast_matches
            print("Done")

            # Test if all parallels have been actually removed
            total_blast_cands = tuple(candidate for match in edited_blast_matches for candidate in match.candidates)

            starmap_args = tuple((match, candidate) for match in edited_blast_matches for candidate in match.candidates
                                 if total_blast_cands.count(candidate) > 1)
            print(f"Remaining parallels (Should be 0): {len(tuple(pool.starmap(func.create_parallel, starmap_args, max_threads)))}")

        # endregion
        # region Remove multiple matches
        print("\nLooking for multiple matches... ", end="")
        multiples_num = len(tuple(1 for match in blast_matches if len(match.candidates) > 1))
        print(f"Done\nFound {multiples_num} multiple matches")

        if multiples_num > 0:
            print("\nFixing multiple matches... ", end="")
            edited_blast_matches = tuple(
                # If there are multiple possible matches:
                DefiniteMatch(possible_match.main_gene, func.get_best_match(possible_match.candidates))
                if len(possible_match.candidates) > 1

                # If not:
                else
                DefiniteMatch(possible_match.main_gene, tuple(possible_match.candidates)[0])

                for possible_match in blast_matches
                if len(possible_match.candidates) > 0
            )

        else:
            print("No multiples to fix. Proceeding")
            edited_blast_matches = tuple(
                DefiniteMatch(possible_match.main_gene, tuple(possible_match.candidates)[0])

                for possible_match in blast_matches
                if not len(possible_match.candidates) == 0
            )

        print(f"\nblast_matches NONEs (Should be empty if everything's working well): {tuple(m.gene_1 for m in edited_blast_matches if not m)[:20]}")
        print(f"Edited Blast matches: {[x.gene_1 for x in edited_blast_matches][:20]}")
        # endregion

        # region Create gene dictionary
        print("\nCreating a gene dictionary... ", end="")
        gene_lexicon = tuple((d_match.gene_1, d_match.gene_2)
                             for d_match in edited_blast_matches
                             if d_match.gene_1 != d_match.gene_2)
        print("Done")
        print(f"\nGene lexicon: {len(gene_lexicon)} items | {gene_lexicon[:20]}")
        # endregion

        # region Free up memory
        blast_matches = None
        blast_matches_extracted = None
        edited_blast_matches = None
        total_blast_cands = None
        parallel_real_matches = None
        parallels = None
        grouped_parallels = None
        starmap_args = None
        results_txt = None
        fasta_by_gene_1 = None
        fasta_by_gene_2 = None
        # endregion

        # Replace gene names in genelist according to blast
        new_genome_1 = Genome(genome_1.annotations, genome_1.sequence)
        new_genome_2 = Genome(genome_2.annotations, genome_2.sequence)
        new_genome_1.replace(gene_lexicon)
        new_genome_2.replace(gene_lexicon)

        print(f"New Genome 1: {tuple(gene.name for gene in new_genome_1.genelist)}")

        # Make a new genelist
        final_genome_1 = new_genome_1
        final_genome_2 = new_genome_2
        # fixme: final double genelist is missing a few genes
        final_double_genelist = DoubleGenelist(new_genome_1, new_genome_2)
        print("Blasting is DONE!")

    else:
        final_genome_1 = genome_1
        final_genome_2 = genome_2
        final_double_genelist = double_genelist

    # print(f"Final double genelist: {final_double_genelist.list[:4]}")

    # region Convert double genelist -> string tuple
    genome_1_names = tuple(gene.name for gene in final_genome_1.genelist)
    genome_2_names = tuple(gene.name for gene in final_genome_2.genelist)

    # final_double_genelist_names = tuple()

    with open(memory_dir + "Final genelist.txt", 'w') as genlst_file:
        to_write = tuple(gene.name for tup in final_double_genelist.list for gene in tup if gene)
        tuple(genlst_file.write(item + ", ") for item in to_write)

    outp_genelist = tuple(
          f"{tup[0].name} {genome_1_names.index(tup[0].name) + 1} {genome_2_names.index(tup[1].name) + 1}"
          if tup[0] and tup[1]

          else
          f"{tup[1].name} {no_gene_symbol} {genome_2_names.index(tup[1].name) + 1}"
          if not tup[0]

          else
          f"{tup[0].name} {genome_1_names.index(tup[0].name) + 1} {no_gene_symbol}"
          if not tup[1]

          else
          print(f"Both None: {tup}")

          for tup in final_double_genelist.list
    )

    # Check
    print("Checking... ", end="")
    splits_1 = tuple(int(line.split()[1]) for line in outp_genelist if no_gene_symbol not in line.split()[1])
    splits_2 = tuple(int(line.split()[2]) for line in outp_genelist if no_gene_symbol not in line.split()[2])
    final_dgl_names = tuple(gene.name for tup in final_double_genelist.list for gene in tup if gene)
    print("Done")

    for i in range(1, sorted(splits_1)[-1]):
        if i not in splits_1:
            print(f"Missing gene in genome 1: {i} {genome_1_names[i-2:i+2]}; Is in final gnlst? {(t := next((tup for tup in final_double_genelist.list if next(x for x in tup if x).name == genome_1_names[i-1]), None))}, {next(item for item in t if item).name}")
        if i not in splits_2:
            print(f"Missing gene in genome 2: {i} {genome_2_names[i-2:i+2]}; Is in final gnlst? {(t := next((tup for tup in final_double_genelist.list if next(x for x in tup if x).name == genome_1_names[i-1]), None))}, {next(item for item in t if item).name}")

    # print(genome_1_names[353])
    # print(genome_2_names[474])
    # print("Unnumeric locations: " + str(tuple(x.split()[1] for x in outp_genelist if not x.split()[1].isnumeric())))

    # endregion

    return (genome_1.name, genome_2.name), outp_genelist


def order_seqs_in_genelist(filename, target_filename="", input_dir=const_input_dir, output_dir=const_output_dir):
    format = filename.split('.')[-1]

    def sort_by_pos(line_):
        if line_.split()[1] == no_gene_symbol:
            return int(line_.split()[2]) ** 2 + 1
        elif line_.split()[2] == no_gene_symbol:
            return int(line_.split()[1]) ** 2 + 2
        else:
            return int(line_.split()[1]) ** 2 + 3

    if target_filename == "":
        target_filename = filename

    # Extract text
    with open(input_dir + filename, 'r') as input_file:
        input_txt = input_file.readlines()

    output_file = open(output_dir + target_filename, 'w')

    # Adjust top-line description
    if input_txt[0].startswith("1"):
        output_file.write(input_txt[0])
    output_file.write("##Location: In List")

    if format == "genelist1":

        # Collect sequences locations
        locations = list()
        for line in input_txt:
            if not line.startswith("1") and not line.startswith("##"):
                locations.append(int(line.split()[1]))

        # Sort locations by size (who comes before in genome)
        locations.sort()
        line_num = 0
        sorted_op_txt = []

        # Add to list
        for line in input_txt:
            line_num += 1
            print("\rProcessing... (" + str(int((line_num / len(locations)) * 100)) + "%)", end="")
            if not line.startswith("1") and not line.startswith("##"):
                # Add genes to sorted_op_txt in sorted order
                for i in range(len(locations)):
                    if locations[i] == int(line.split()[1]):
                        sorted_op_txt.append("\n" + line.split()[0] + " " + str(i + 1))
                        break

        # Sort and write output file
        # sorted_op_txt.sort(key=sort_by_pos)
        for line in sorted_op_txt:
            output_file.write(line)

    elif format == "genelist2":

        # Collect sequences locations
        locations1 = list()
        locations2 = list()
        for line in input_txt:
            if not line.startswith("1") and not line.startswith("##"):
                if no_gene_symbol not in line.split()[1]:
                    locations1.append(int(line.split()[1]))
                if no_gene_symbol not in line.split()[2]:
                    locations2.append(int(line.split()[2]))

        # Sort locations by size (who comes before in genome)
        locations1.sort()
        locations2.sort()
        line_num = 0
        sorted_op_txt = []

        # Add to list
        for line in [x for x in input_txt
                     if not x.startswith("1") and
                        not x.startswith("##")]:
            line_num += 1
            if len(locations1) > 0:
                print("\rProcessing... (" + str(int((line_num / len(locations1)) * 100)) + "%)", end="")
            # if not line.startswith("1") and not line.startswith("##"):

            # Locations1 (first genome location number)
            if True:
                if no_gene_symbol not in line.split()[1:3]:
                    # Add genes to sorted_op_txt in sorted order
                    for i1, location1 in enumerate(locations1):
                        if location1 == int(line.split()[1]):

                            # Locations2 (second genome location number)
                            for i2, location2 in enumerate(locations2):
                                if location2 == int(line.split()[2]):
                                    # print(f"i2: {i2}")
                                    sorted_op_txt.append(f"\n{line.split()[0]} {i1 + 1} {i2 + 1}")
                                    break

                else:
                    sorted_op_txt.append("\n" + line.removesuffix("\n"))

        # Sort and write output file
        # sorted_op_txt.sort(key=sort_by_pos)

        # Double-check order
        i1 = 1
        i2 = 1

        # todo: Remove this if statement when experiments are done
        if True:
            for num, line in enumerate(sorted_op_txt):
                if line.split()[1] == no_gene_symbol:
                    if not int(line.split()[2]) == i2:
                        sorted_op_txt[num] = f"\n{line.split()[0]} {line.split()[1]} {i2}"
                    i2 += 1
                else:
                    if not int(line.split()[1]) == i1:
                        sorted_op_txt[num] = f"\n{line.split()[0]} {i1} {line.split()[2]}"
                    i1 += 1

        for line in sorted_op_txt:
            output_file.write(line)

    output_file.close()


def close_program():
    exit()


if __name__ == "__main__":
    func.delete_old_logs(log_dir, limit=10)
    main_menu()
