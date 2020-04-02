import sys
import statistics
import math

path = sys.argv[1]
filename = sys.argv[2]
scene_name = sys.argv[3]
output = sys.argv[4]
material = sys.argv[5]

configurations = {}
styles = {}

if material == "diffuse":
	configurations = {
		"AMC" : ["amr_500", "amr_1000", "amr_2000", "amr_4000", "amr_8000"],
		"MDLC" : ["mdlc_250", "mdlc_500", "mdlc_1000", "mdlc_2000", "mdlc_4000"],
		"LS" : ["ls_250", "ls_500", "ls_1000", "ls_2000", "ls_4000"],
		"Mat Sep" : ["amr_500", "amr_1000", "amr_2000", "amr_4000", "amr_8000"]
	}

	styles = {
		"AMC" : ("blue", "square"),
		"MDLC" : ("red", "square"),
		"LS" : ("green", "square"),
		"Mat Sep" : ("black", "square")
	}
else:
	configurations = {
		"AMC" : ["amr_500", "amr_1000", "amr_2000", "amr_4000", "amr_8000"],
		"MDLC" : ["mdlc_250", "mdlc_500", "mdlc_1000", "mdlc_2000", "mdlc_4000"],
		"LS" : ["ls_250", "ls_500", "ls_1000", "ls_2000", "ls_4000"]
	}

	styles = {
		"AMC" : ("blue", "square"),
		"MDLC" : ("red", "square"),
		"LS" : ("green", "square")
	}

coords = {}
stats = {}

max_x = -99999
max_y = -99999

for config_title, config_codes in configurations.items():
	coords[config_title] = []
	stats[config_title] = []
	for config_code in config_codes:
		timings_filename = path + "/timings_" + config_code
		error_filename = path + "/" + filename + "_" + config_code + "_errors.txt"
		samplerates_filename = path + "/samplerates_" + config_code

		timings = []
		errors = []
		samplerates = []

		total_timings = 0
		with open(timings_filename, "r") as timings_file:
			total_num_timings = 0

			for y in timings_file.read().split('\n'):
				try:
					val = float(y)
					max_x = max(max_x, val)
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
					max_y = max(max_y, val)
					total_errors += val
					errors.append(val)
					total_num_errors += 1
				except:
					pass

			total_errors /= total_num_errors

		total_samplerates = 0
		with open(samplerates_filename, "r") as samplerates_file:
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
		samplerates_standard_dev = statistics.stdev(samplerates) if len(samplerates) > 1 else 0

		ebar_x = time_standard_dev / math.sqrt(len(timings))
		ebar_y = error_standard_dev / math.sqrt(len(errors))

		coords[config_title].append((total_timings, total_errors, ebar_x, ebar_y))
		stats[config_title].append((config_code, total_timings, time_standard_dev, total_errors, error_standard_dev, total_samplerates, samplerates_standard_dev))

output_str = ""
output_str += "\\begin{figure}\n"
output_str += "\\centering\n"
output_str += "\\begin{tikzpicture}\n"
output_str += "\\begin{axis}[\n"
output_str += "title={" + scene_name + "},\n"
output_str += "xlabel={Time [s]},\n"
output_str += "ylabel={Error [rmse]},\n"
output_str += "xmin=0, xmax=" + str(max_x) + ",\n"
output_str += "ymin=0, ymax=" + str(max_y) + ",\n"
output_str += "legend pos=north east,\n"
output_str += "ymajorgrids=true,\n"
output_str += "error bars/y dir      = both,\n"
output_str += "error bars/y explicit = true,\n"
output_str += "error bars/x dir      = both,\n"
output_str += "error bars/x explicit = true,\n"
output_str += "]\n\n"

for config_title, ccc in coords.items():
	output_str += "\\addplot[\n"
	output_str += "color=" + styles[config_title][0] + ",\n"
	output_str += "mark=" + styles[config_title][1] + ",\n"
	output_str += "]\n"
	output_str += "coordinates {\n"
	for cc in ccc:
		output_str += "(" + str(cc[0]) + ", " + str(cc[1]) + ") +- (" + str(cc[2]) + ", " + str(cc[3]) + ")"
	output_str += "\n"
	output_str += "};\n"
	output_str += "\\addlegendentry{" + config_title + "}\n"

output_str += "\\end{axis}\n"
output_str += "\\end{tikzpicture}\n"
output_str += "\\end{figure}\n"

output_table_str = ""
for config_title, ccc in stats.items():
	for cc in ccc:
		output_table_str += config_title + " " + cc[0] + ":\n"
		output_table_str += "Time: " + str(cc[1]) + " (" + str(cc[2]) + ")\n"
		output_table_str += "Error: " + str(cc[3]) + " (" + str(cc[4]) + ")\n"
		output_table_str += "Samples: " + str(cc[5]) + " (" + str(cc[6]) + ")\n"
		output_table_str += "\n"
	output_table_str += "\n\n"

o = open(output + ".tex", "w")
o.write(output_str)
o.close()

t = open(output + ".stats", "w")
t.write(output_table_str)
t.close()