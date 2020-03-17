import sys
import statistics
import math

path = sys.argv[1]
filename = sys.argv[2]
scene_name = sys.argv[3]
output = sys.argv[4]

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

min_x = sys.float_info.max
max_x = -sys.float_info.max
min_y = sys.float_info.max
max_y = -sys.float_info.max

for config_title, config_codes in configurations.items():
	coords[config_title] = []
	for config_code in config_codes:
		timings_filename = path + "/timings_" + config_code
		error_filename = path + "/" + filename + "_" + config_code + "_errors.txt"

		timings = []
		errors = []

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

		min_x = min(min_x, total_timings)
		max_x = max(max_x, total_timings)
		min_y = min(min_y, total_errors)
		max_y = max(max_y, total_errors)

		time_standard_dev = statistics.stdev(timings) if len(timings) > 1 else 0
		error_standard_dev = statistics.stdev(errors) if len(errors) > 1 else 0

		ebar_x = time_standard_dev / math.sqrt(len(timings))
		ebar_y = error_standard_dev / math.sqrt(len(errors))

		coords[config_title].append((total_timings, total_errors, ebar_x, ebar_y))

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

o = open(output, "w")
o.write(output_str)
o.close()