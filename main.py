from libs_and_dirs import *
import func


def main_menu():
    print(
        "Choose an action to apply on the files in the input folder:\n"
        "0. Close program\n"
        "1. Convert\n"
        "2. Analyze\n"
        "3. Text Editor\n"
        "4. Clean Folder"
    )
    command = input()
    if command.isnumeric():
        command = int(command)

        if command == 0:
            close_program()
        elif command == 1:
            convert_menu()
        elif command == 2:
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
        "\n GENELIST1:\n\n"
        "2. -> GENELIST2\n"
        "3. -> GENELIST1 numbered according to location in list"
    )
    command = input()
    if command.isnumeric():
        command = int(command)

        if command == 0:
            main_menu()
        elif command == 1:
            convert_genomes_to_genelist2()
            main_menu()
        elif command == 2:
            print("This command does not exist yet!")
            convert_menu()
        elif command == 3:
            # Get input
            file = func.show_and_get_files(1, "Choose a genelist(1/2) file:\n",
                                           acceptable_formats=["genelist1", "genelist2"])
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
            func.apply_algorithm()
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


def convert_genomes_to_genelist2() -> None:
    """
    Takes GFF & FASTA files of both genomes from input directory, and converts them to a single GENELIST2 file ready
    for SI comparison.\n
    :return: NONE
    """
    start = time.perf_counter()

    do_blast = False
    # Clean Folders
    func.delete_files(memr_inp_dir)
    func.delete_files(memr_op_dir)
    func.delete_files(temp_dir)

    # Validate input files
    chosen_gff3_name_1 = func.find_genome(file_1_dir)[0]
    chosen_fasta_name_1 = func.find_genome(file_1_dir)[1]
    chosen_gff3_name_2 = func.find_genome(file_2_dir)[0]
    chosen_fasta_name_2 = func.find_genome(file_2_dir)[1]

    print(f"File 1: GFF3= {func.find_genome(file_1_dir)[0]}\n"
          f"File 1: Fasta= {func.find_genome(file_1_dir)[1]}\n"
          f"File 2: GFF3= {func.find_genome(file_2_dir)[0]}\n"
          f"File 2: Fasta=  {func.find_genome(file_2_dir)[1]}")
    logging.debug(f"File 1: GFF3= {func.find_genome(file_1_dir)[0]}\n"
                  f"File 1: Fasta= {func.find_genome(file_1_dir)[1]}\n"
                  f"File 2: GFF3= {func.find_genome(file_2_dir)[0]}\n"
                  f"File 2: Fasta=  {func.find_genome(file_2_dir)[1]}")

    # Convert GFF3 -> GENELIST1
    func.create_genelist1(chosen_gff3_name_1, gene_type="Name",
                          input_dir=file_1_dir, output_dir=memr_op_dir, order=func.Order.Location)
    func.create_genelist1(chosen_gff3_name_2, gene_type="Name",
                          input_dir=file_2_dir, output_dir=memr_op_dir, order=func.Order.Location)

    # Convert GENELIST1 -> GENELIST2
    do_blast = func.merge_genelist(
        chosen_gff3_name_1.replace("." + chosen_gff3_name_1.split('.')[-1], func.conversion_1_output_format),
        chosen_gff3_name_2.replace("." + chosen_gff3_name_2.split('.')[-1], func.conversion_1_output_format),
        input_dir=memr_op_dir,
        output_dir=memr_inp_dir
    )

    if do_blast:
        # Convert GFF3 & Fasta -> Unidentified Genes Sequences (Fasta)
        func.separate_by_gene(gff3=chosen_gff3_name_1, fasta=chosen_fasta_name_1, output_filename="Genome1",
                              input_dir=file_1_dir, output_dir=blast_dir, memory_dir=memory_dir,
                              one_outp_file=True, specific_genes_only=False, specific_genes=func.unfamiliar_genes)
        func.separate_by_gene(gff3=chosen_gff3_name_2, fasta=chosen_fasta_name_2, output_filename="Genome2",
                              input_dir=file_2_dir, output_dir=blast_dir, memory_dir=memory_dir,
                              one_outp_file=True, specific_genes_only=False, specific_genes=func.unfamiliar_genes)

        # Blast
        func.blast(blast_dir=blast_dir, is_protein=False)

        # Change GENELIST1s according to blast results
        pool = mp.Pool()
        if __name__ == "__main__":
            # #Creating a bilateral (Can be accessed by either keys or values) dictionary of sequences
            # #seq_dicts = [func.blast_results_to_dict([p1, ])]

            #
            with open(blast_dir + "results.txt", 'r') as results_file:
                results_txt = results_file.readlines()

            args_list = []
            for counter_num, line in enumerate(results_txt):
                args_list.append((line, counter_num, len(results_txt), temp_dir), )

            print(f"args_list includes {len(args_list)} items")

            print("\n\033[4m\033[1m" "Matching unrecognized genes" "\033[0m")
            print(f"Total number of lines: {len(results_txt)}")

            counter_list = list(pool.starmap(func.process_blast_results_line, args_list))
            # counter_list = list(itertools.starmap(func.process_blast_results_line, args_list))

            counter_list = [x for x in counter_list if x]

            print(f"\n\nUnmatched Genes: {[x[0] for x in counter_list[0:30]]}")

            # print(f"Counter list: {counter_list[0:6]}")
            # func.process_blast_results_line(args_list[0], args_list[1], args_list[2], args_list[3])

            # Read counter list from file
            # with open(temp_dir + "Counter List.txt", 'r') as counter_list_f:
            #     counter_list_txt = counter_list_f.read()

            # for line in counter_list_txt.splitlines():
            #     counter_obj = Counter()
            #     counter_obj[line.split()[1]] = line.split()[2].removesuffix(',')
            #     if len(line.split(',')) > 0:
            #         for line_num, match in enumerate(line.split(',')):
            #             if line_num != 0:
            #                 counter_obj[match.split()[0]] = match.split()[1]
            #
            #     counter_list.append([line.split()[0], counter_obj])
            #     print("\rItems in counter_list: " + str(len(counter_list)), end="")
            #
            # print(f"\nCounter list: {counter_list[0:5]}...")
            # print(f"Counter Items: {[x[0] for x in counter_list]}")
            #

            # for item in [x[1] for x in counter_list]:
            #     cl_subjects.extend(item)
            # cl_counters.extend([dict for counter in counter_list for dict in counter[1]])

            # cl_subjects.extend([x for x in [x[1] for x in counter_list]])

            # Adjust multiple matches problems ↓ ↓ ↓ ↓ ↓
            # Create 1 list of possible matches
            cl_counters = [dict for counter in counter_list for dict in counter[1]]  # Counter list subjects
            print(f"cl_counters done ({len(cl_counters)})")
            # parallels = [{x: cl_counters.count(x)} for x in cl_counters if cl_counters.count(x) > 1]
            # parallels = [{x: cl_counters.count(x)} for x in cl_counters if cl_counters.count(x) > 1]
            parallels = list(zip(
                [counter for counter in cl_counters if cl_counters.count(counter) > 1],
                [cl_counters.count(counter) for counter in cl_counters if cl_counters.count(counter) > 1]
            ))
            parallels_temp_1 = [counter for counter in cl_counters if cl_counters.count(counter) > 1]
            parallels_temp_2 = [cl_counters.count(counter) for counter in cl_counters if cl_counters.count(counter) > 1]
            parallels = parallels_temp_1.

            print(f"length: {len(cl_counters)}; {len(parallels)}")

            # print(f"Counter Keys: {cl_counters[0:20]}")

            if len(parallels) > 0:
                print(f"Parallels: ({len(parallels)}) | {parallels[0:100]}")
                parals_more_than_1 = [x for x in parallels if len(x) > 1]
                print(f"Parallels > 1: ({len(parals_more_than_1)}) | {parals_more_than_1[0:100]}")

            # Add all matches to the list | Fix parallels
            parals_to_remove = []
            to_delete = []
            temp_parals_to_compare = []

            parallels_args_list_ = [x for x in enumerate(parallels)]
            parallels_args_list = []
            length = len(parallels)
            # for item in parallels_args_list_:
            #     parallels_args_list.append((item[0], item[1], length))
            parallels_args_list = [(parallel[0], parallel[1], length) for parallel in parallels_args_list_]
            print(f"parallels_args_list: {parallels_args_list[0:3]}; length: {len(parallels_args_list)}")

            starmap_results = pool.starmap(func.fix_parallels, parallels_args_list)
            temp_parals_to_compare_, temp_to_delete = [x[0] for x in starmap_results], [x[1] for x in starmap_results]
            temp_parals_to_compare.extend(temp_parals_to_compare_)
            to_delete.extend(temp_to_delete)

            print(f"\n\ntemp_list: {temp_parals_to_compare[0:10]}; objs: {temp_parals_to_compare_[0:10]}")
            print(f"to_delete: {to_delete[0:10]}; objs: {temp_to_delete[0:10]}")

            to_delete = [trash for list_ in to_delete for trash in list_]
            print(f"To delete: {to_delete[0:8]}")
            for trash in to_delete:
                counter_list[int(trash.split('-')[-1])].pop(trash.split('-')[0])

            parals_to_remove.extend(temp_parals_to_compare)

            # Remove duplicates
            print(f"\n\n\nparallels_to_compare: {parals_to_remove[0:10]}\n")
            print("\nRemoving duplicates...")
            # for num, duplicate_type in enumerate(parals_to_compare):
            #     print(f"\rRemoving duplicates: {int(100 * (num / len(parals_to_compare)))}%", end="")
            #     # print(f"collection: {collection}")
            #     if len(duplicate_type) > 0:
            #         duplicate_type.sort(reverse=True, key=lambda collec: float(list(collec.values())[0]))
            #         # collection.sort(reverse=True, key=lambda e: float(list()[0]))
            #         # collection.sort(reverse=True, key=lambda e: float(e[0]))
            #         duplicate_type.remove(duplicate_type[0])
            #     else:
            #         parals_to_compare.remove(duplicate_type)

            # Sort duplicate collections
            parals_to_remove = [
                sorted(duplicate_type, reverse=True, key=lambda dupl: float(list(dupl.values())[0]))
                for duplicate_type in parals_to_remove
                if len(duplicate_type) > 0
            ]

            # Remove first item in each collection (didn't work without range:()
            for paral in parals_to_remove:
                paral.pop(0)

            # Unpack list
            print(f"{parals_to_remove[:10]=}")
            # parals_to_remove = [dict(seq_type.items()) for seq_type in parals_to_remove]
            parals_to_remove = [parallel for seq_type in parals_to_remove for parallel in seq_type]
            print(f"{parals_to_remove[:10]=}")

            # for i in range(len(parals_to_remove)):
            #     # parals_to_compare[i] = [item for item in parals_to_compare[i] if item != parals_to_compare[i][0]]
            #     parals_to_remove[i].pop(0)

            for counter_num, counter in enumerate([counter[1] for counter in counter_list]):
                if counter_num % 15 == 0:
                    print(f"\rIntersection {int(100 * (counter_num / len([counter[1] for counter in counter_list])))}%",
                          end="")
                # intersection = set(counter.items()) & set([list(x.items()) for x in parals_to_remove])
                # {x.items() for x in parals_to_remove}
                # list(set(counter.items()))
                # intersection = list(set(counter.items()) & {x.items() for x in parals_to_remove})
                intersection = [list(paral.items())[0] for paral in parals_to_remove if list(paral.items())[0][0] in counter.items()]
                for parallel in [paral for paral in intersection if paral in counter]:
                    print(counter_list[counter_num][1].pop(parallel))

            # intersection = [list(paral.items())[0] for paral in parals_to_remove]
            # for i in range(len(counter_list)):
            #     try:
            #         counter_list[i][1].pop(parallel)

            print("Successfully removed duplicates")

                # print(f"dict_counter: {len(dicts_counter)} | List: {list(dicts_counter)}")
                # counter_list[line_num][1].pop(list([seq_match.keys() for seq_match in list(match_dict)
                #                                     if any(
                #         {seq_match: match_dict[seq_match]} in x for x in parals_to_compare)])[0])


                # for dict_item in list(dicts_counter):
                #     dict_ = {dict_item: dicts_counter[dict_item]}
                #
                #     if any(dict_ in x for x in parals_to_compare):
                #         counter_list[line_num][1].pop(list(dict_.keys())[0])
                # print(f"dicts_counter: {dicts_counter}")

            # counter_list.pop(
            #     [dict_item in counter_list for line_num, dicts_counter in counter_list for dict_item in list(dicts_counter)
            #      if any({dict_item: dicts_counter[dict_item]} in x for x in parals_to_compare)
            #
            # ])

            print(f"New Counter List: {counter_list[0: 10]}")

            # Adjust multiple genes problem ↑ ↑ ↑ ↑ ↑

            seq_dicts = [func.make_gene_dict(counter_list)]
            # print(seq_dicts[0].items[0:100])
            # fixme: seq dictionaries aren't equal
            print(f"zipped 1: {list(zip([x[1] for x in seq_dicts[0].items()], seq_dicts[0].keys()))}")
            print(f"zipped 2: {list(zip(seq_dicts[0].keys(), [x[1] for x in seq_dicts[0].items()]))}")

            seq_dicts.append(dict(zip(
                [x[1] for x in seq_dicts[0].items()],
                list(seq_dicts[0].keys())
            )))
            # print(f"\n\nseq_dicts[1]: {[x[1] for x in seq_dicts[0].items()]}")
            print(f"\nValues: {[x[1] for x in seq_dicts[0].items()][:10]}")
            print(f"\nKeys: {seq_dicts[0].keys()[:10]}")
            print(f"\n\nGene Dictionary includes[0]: {seq_dicts[0]}\n")
            print(f"\n\nGene Dictionary includes[1]: {seq_dicts[1]}\n")
            # logging.debug(f"seq_dicts (Gene dictionary): \n{seq_dicts[0][0:100]}")

        # Replace sequence names in GENELIST1 files according to the new data
        # Find original GENELIST1s
        for path, dirs, files in os.walk(memr_op_dir):
            for file in files:
                if file.lower().endswith(".genelist1"):
                    gnlst_f1 = open(memr_op_dir + file, 'r')
                    break

        for path, dirs, files in os.walk(memr_op_dir):
            for file in files:
                if file.lower().endswith(".genelist1") and file != gnlst_f1.name.split('\\')[-1]:
                    gnlst_f2 = open(memr_op_dir + file, 'r')
                    break

        gnlst_txt1 = gnlst_f1.read()
        gnlst_txt2 = gnlst_f2.read()

        # for line in gnlst1_txt2.splitlines():
        #     if line.split()[0] in seq_dicts[0]:
        #         gnlst1_txt2 = gnlst1_txt2.replace(line.split()[0], line.split()[0] + "|" + seq_dicts[0][line.split()[0]])
        # for line in gnlst1_txt1.splitlines():
        #     if seq_dicts[1].__contains__(line.split()[0]):
        #         gnlst1_txt1 = gnlst1_txt1.replace(line.split()[0], seq_dicts[1][line.split()[0]] + "|" + line.split()[0])

        # Change names
        gnlst_txt1 = "\n".join([line.replace(line.split()[0], f"{seq_dicts[1][line.split()[0]]}|{line.split()[0]}")
                                for line in gnlst_txt1.splitlines() if line.split()[0] in seq_dicts[1]])
        gnlst_txt2 = "\n".join([line.replace(line.split()[0], f"{line.split()[0]}|{seq_dicts[0][line.split()[0]]}")
                                for line in gnlst_txt2.splitlines() if line.split()[0] in seq_dicts[0]])

        with open(memr_inp_dir + "REBUILT-" + gnlst_f1.name.split('\\')[-1], 'w') as outp_gnlst_f1:
            outp_gnlst_f1.write(gnlst_txt1)
        with open(memr_inp_dir + "REBUILT-" + gnlst_f2.name.split('\\')[-1], 'w') as outp_gnlst_f2:
            outp_gnlst_f2.write(gnlst_txt2)

        gnlst_f1.close()
        gnlst_f2.close()

        # Convert new GENELIST1 -> GENELIST2
        func.merge_genelist \
                (
                "REBUILT-" + chosen_gff3_name_1.replace("." + chosen_gff3_name_1.split('.')[-1],
                                                        func.conversion_1_output_format),
                "REBUILT-" + chosen_gff3_name_2.replace("." + chosen_gff3_name_2.split('.')[-1],
                                                        func.conversion_1_output_format),
                input_dir=memr_inp_dir,
                output_dir=memr_op_dir,
                asterisks=True
            )

    # Fix location numbers in new GENELIST2 file
    if do_blast:
        order_seqs_in_genelist(filename="Merged-" +
                                        "REBUILT-" + chosen_gff3_name_1.replace("." + chosen_gff3_name_1.split('.')[-1],
                                                                                "") + "-" +
                                        "REBUILT-" + chosen_gff3_name_2.replace("." + chosen_gff3_name_2.split('.')[-1],
                                                                                "") + ".genelist2",
                               target_filename=chosen_gff3_name_1.replace("." + chosen_gff3_name_1.split('.')[-1],
                                                                          "") + "-" +
                                               chosen_gff3_name_2.replace("." + chosen_gff3_name_2.split('.')[-1],
                                                                          "") + ".genelist2",
                               input_dir=memr_op_dir, output_dir=const_output_dir)
    else:
        # fixme: Here lies the problem
        order_seqs_in_genelist(filename="Merged-" +
                                        chosen_gff3_name_1.replace("." + chosen_gff3_name_1.split('.')[-1],
                                                                   "") + "-" +
                                        chosen_gff3_name_2.replace("." + chosen_gff3_name_2.split('.')[-1],
                                                                   "") + ".genelist2",
                               target_filename=chosen_gff3_name_1.replace("." + chosen_gff3_name_1.split('.')[-1],
                                                                          "") + "-" +
                                               chosen_gff3_name_2.replace("." + chosen_gff3_name_2.split('.')[-1],
                                                                          "") + ".genelist2",
                               input_dir=memr_inp_dir, output_dir=const_output_dir)

    # Stop timer
    end = time.perf_counter()
    time_passed = end - start
    if time_passed > 60:
        minutes = int(time_passed / 60)
        seconds = time_passed % 60
        logging.debug(f"\nTook {minutes} minutes and {seconds} seconds to finish operation")
        print(f"\nTook {minutes} minutes and {seconds} seconds to finish operation")
    else:
        logging.debug(f"\nTook {time_passed} seconds to finish operation")
        print(f"\nTook {time_passed} seconds to finish operation")


def order_seqs_in_genelist(filename, target_filename="", input_dir=const_input_dir, output_dir=const_output_dir):
    files_dict = dict()
    format = filename.split('.')[-1]

    def sort_by_pos(line_):
        if line_.split()[1] == "X":
            return int(line_.split()[2]) ** 2 + 1
        elif line_.split()[2] == "X":
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
                if "X" not in line.split()[1]:
                    locations1.append(int(line.split()[1]))
                if "X" not in line.split()[2]:
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
                if "X" not in line.split()[1:3]:
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
                if line.split()[1] == "X":
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


def create_genelist1(filename, input_dir_=memr_inp_dir, output_dir_=memr_op_dir):
    from func import create_genelist1
    create_genelist1(filename, gene_type="Name", input_dir=input_dir_, output_dir=output_dir_)


def close_program():
    exit()


if __name__ == "__main__":
    func.delete_old_logs(log_dir, limit=10)
    main_menu()
