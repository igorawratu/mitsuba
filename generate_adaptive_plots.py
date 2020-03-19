import sys
import statistics
import math

path = sys.argv[1]
filename = sys.argv[2]
scene_name = sys.argv[3]
output = sys.argv[4]

configurations = ["adaptive", "static_01", "static_04"]

adaptive_dat = {}
adaptive_dat["sample percentage"] = []
adaptive_dat["time"] = []
adaptive_dat["error"] = []

for config_code in configurations:
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

    adaptive_dat["sample percentage"].append((total_samplerates, samplerate_standard_dev))
    adaptive_dat["time"].append((total_timings, time_standard_dev))
    adaptive_dat["error"].append((total_errors, error_standard_dev))

output_str = ""
output_str += "\\begin{center}\n"
output_str += "\\begin{tabular}{c c c c}\n"
output_str += " & Adaptive & Static 0.1 & Static 0.4\n"

output_str += "sample percentage"
for dat in adaptive_dat["sample percentage"]:
    output_str += " & " + str(adaptive_dat["sample percentage"][0]) + "\\% (" + str(adaptive_dat["sample percentage"]) + "\\%)"
output_str += "\n"

output_str += "time"
for dat in adaptive_dat["time"]:
    output_str += " & " + str(adaptive_dat["time"][0]) + "s (" + str(adaptive_dat["time"]) + "s)"
output_str += "\n"

output_str += "error (rmse)"
for dat in adaptive_dat["error"]:
    output_str += " & " + str(adaptive_dat["error"][0]) + " (" + str(adaptive_dat["error"]) + ")"
output_str += "\n"

output_str += "\\end{tabular}\n"
output_str += "\\end{center}\n"

o = open(output, "w")
o.write(output_str)
o.close()

