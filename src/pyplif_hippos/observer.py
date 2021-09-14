import sys


subscribers = {}
file_handlers = {}


def subscribe(event_type, task):
    if event_type not in subscribers:
        subscribers[event_type] = []

    subscribers[event_type].append(task)


def do_task(event_type, *args):
    if event_type not in subscribers:
        return

    for task in subscribers[event_type]:
        arg_count = task.__code__.co_argcount
        task(*args[0:arg_count])


def handle_open_file(file_name):
    file_handlers[file_name] = open(file_name, "w")


def handle_initial_log(config, ligand_pose):
    ligand_name = ligand_pose[0].split("_")[0]
    pose_count = len(ligand_pose)
    file_handle = file_handlers["log"]
    file_handle.write("Ligand name is %s with %s poses\n\n" % (ligand_name, pose_count))
    if "similarity" in file_handlers:
        coef = ", ".join(config.similarity_coef)
        file_handle.write("similarity coefficient used are %s\n\n" % coef)
    if "simplified" in file_handlers:
        file_handle.write("RESNAME length startbit endbit\n")


def handle_missing_docking_output(docking_results):
    if len(docking_results["docked_ligands"]) == 0:
        missing_docking_output = (
            "The docking output could not be found. Please check your docking result."
        )

        print(missing_docking_output)
        file_handle = file_handlers["log"]
        file_handle.write(missing_docking_output)

        sys.exit(1)


def handle_bit_variation(resname, bitlength, bit_start, bit_end):
    file_handle = file_handlers["log"]
    file_handle.write("%-10s %-6s %-7s %s\n" % (resname, bitlength, bit_start, bit_end))


def handle_write_simplified(ligand_name, score, all_bits):
    file_handle = file_handlers["simplified"]
    file_handle.write("%s %s %s\n" % (ligand_name, score, all_bits[0]))


def handle_write_full(ligand_name, score, all_bits):
    file_handle = file_handlers["full"]
    file_handle.write("%s %s %s\n" % (ligand_name, score, all_bits[1]))


def handle_write_full_nobb(ligand_name, score, all_bits):
    file_handle = file_handlers["full_nobb"]
    file_handle.write("%s %s %s\n" % (ligand_name, score, all_bits[2]))


def handle_write_similarity(ligand_name, coefficient):
    file_handle = file_handlers["similarity"]
    file_handle.write("%s %s\n" % (ligand_name, " ".join(coefficient)))


def handle_bistring_error(bistring_error):
    file_handle = file_handlers["log"]
    file_handle.write(bistring_error)


def handle_write_time(total_time):
    file_handle = file_handlers["log"]
    file_handle.write("\nTotal time taken %.3f s.\n" % total_time)


def close_simplified():
    file_handlers["simplified"].close()


def close_full():
    file_handlers["full"].close()


def close_full_nobb():
    file_handlers["full_nobb"].close()


def close_similarity():
    file_handlers["similarity"].close()


def close_logfile():
    file_handlers["log"].close()


def setup_simplified(config):
    file_handlers["simplified"] = open(config.simplified_outfile, "w")
    subscribe("write_bitstrings", handle_write_simplified)
    subscribe("close_files", close_simplified)


def setup_full(config):
    file_handlers["full"] = open(config.full_outfile, "w")
    subscribe("write_bitstrings", handle_write_full)
    subscribe("close_files", close_full)


def setup_full_nobb(config):
    file_handlers["full_nobb"] = open(config.full_nobb_outfile, "w")
    subscribe("write_bitstrings", handle_write_full_nobb)
    subscribe("close_files", close_full_nobb)


def setup_similarity(config):
    file_handlers["similarity"] = open(config.sim_outfile, "w")
    subscribe("write_similarity", handle_write_similarity)
    subscribe("close_files", close_similarity)


def setup_logfile(config):
    file_handlers["log"] = open(config.logfile, "w")
    subscribe("initial_information", handle_initial_log)
    subscribe("is_docking_output_missing", handle_missing_docking_output)
    subscribe("simplified_bit_log", handle_bit_variation)
    subscribe("biststring_error", handle_bistring_error)
    subscribe("write_time", handle_write_time)
    subscribe("close_files", close_logfile)


setup_dict = dict(
    simplified=setup_simplified,
    full=setup_full,
    full_nobb=setup_full_nobb,
    similarity=setup_similarity,
    logfile=setup_logfile
)
