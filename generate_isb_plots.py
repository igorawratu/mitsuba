import sys
import statistics
import math

path = sys.argv[1]
filename = sys.argv[2]
scene_name = sys.argv[3]
output = sys.argv[4]

configurations = {
    "Boolean" : ["boolge_1_1000", "boolge_15_1000", "boolge_2_1000", "boolge_1_2000", "boolge_1_3000"],
    #["bool_1_1000", "nbool_1_1000", "boolge_1_1000", "bool_15_1000", "nbool_15_1000", "boolge_15_1000", 
    #"bool_2_1000", "nbool_2_1000", "boolge_2_1000", "bool_1_2000", "nbool_1_2000", "boolge_1_2000", "bool_1_3000", "nbool_1_3000", "boolge_1_3000"],
    "Importance Sampling" : ["import_05", "nimport_05", "import_1", "nimport_1", "import_15", "nimport_15", "import_2", "nimport_2"]
}

bool_dat = {}
import_dat = {}

bool_dat["sample percentage"] = []
bool_dat["time"] = []
bool_dat["error"] = []
import_dat["sample percentage"] = []
import_dat["time"] = []
import_dat["error"] = []

for config_title, config_codes in configurations.items():
    for config_code in config_codes:
        timings_filename = path + "/timings_" + config_code
        error_filename = path + "/" + filename + "_" + config_code + "_errors.txt"
        samplerate_filename = path + "/samplerates_" + config_code

        timings = []
        errors = []
        samplerates = []

        total_timings = 0
        with open(timings_filename, "r") as timings_file:
            total_num_timings = 0

            for y in timings_file.read().split('\n'):
                try:
                    val = float(y)
                    total_timings += val
                    timings.append(val)
                    total_num_timings += 1
                except:
                    pass

            total_timings /= total_num_timings

        total_errors = 0
        with open(error_filename, "r") as errors_file:
            total_num_errors = 0

            for y in errors_file.read().split('\n'):
                try:
                    val = float(y)
                    total_errors += val
                    errors.append(val)
                    total_num_errors += 1
                except:
                    pass

            total_errors /= total_num_errors

        total_samplerates = 0
        with open(samplerate_filename, "r") as samplerates_file:
            total_num_samplerates = 0

            for y in samplerates_file.read().split('\n'):
                try:
                    val = float(y)
                    total_samplerates += val
                    samplerates.append(val)
                    total_num_samplerates += 1
                except:
                    pass

            total_samplerates /= total_num_samplerates

        time_standard_dev = statistics.stdev(timings) if len(timings) > 1 else 0
        error_standard_dev = statistics.stdev(errors) if len(errors) > 1 else 0
        samplerate_standard_dev = statistics.stdev(samplerates) if len(samplerates) > 1 else 0

        if config_title == "Boolean":
            bool_dat["sample percentage"].append((total_samplerates, samplerate_standard_dev))
            bool_dat["time"].append((total_timings, time_standard_dev))
            bool_dat["error"].append((total_errors, error_standard_dev))
        else:
            import_dat["sample percentage"].append((total_samplerates, samplerate_standard_dev))
            import_dat["time"].append((total_timings, time_standard_dev))
            import_dat["error"].append((total_errors, error_standard_dev))

output_str = ""
output_str += "Boolean:\n"
for i in range(0, len(configurations["Boolean"])):
    output_str += "samplerate: " + str(bool_dat["sample percentage"][i][0]) + "(" + str(bool_dat["sample percentage"][i][1]) + ")"
    output_str += "time: " + str(bool_dat["time"][i][0]) + "(" + str(bool_dat["time"][i][1]) + ")"
    output_str += "error: " + str(bool_dat["error"][i][0]) + "(" + str(bool_dat["error"][i][1]) + ")\n"

output_str += "\n"

output_str += "Importance Sampling:\n"
for i in range(0, len(configurations["Importance Sampling"])):
    output_str += "samplerate: " + str(import_dat["sample percentage"][i][0]) + "(" + str(import_dat["sample percentage"][i][1]) + ")"
    output_str += "time: " + str(import_dat["time"][i][0]) + "(" + str(import_dat["time"][i][1]) + ")"
    output_str += "error: " + str(import_dat["error"][i][0]) + "(" + str(import_dat["error"][i][1]) + ")\n"

o = open(output, "w")
o.write(output_str)
o.close()
