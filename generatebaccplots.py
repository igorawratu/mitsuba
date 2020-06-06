import sys
import statistics
import math

path = sys.argv[1]
filename = sys.argv[2]
scene_name = sys.argv[3]
output = sys.argv[4]

configurations = {
    "BAMC" : "boolacc",
    "AMC" : "nboolacc",
    "BAMC (GE)" : "boolgeacc"
}

styles = {
    "BAMC" : ("blue", "square"),
    "AMC" : ("red", "square"),
    "BAMC (GE)" : ("green", "square"),
}

coords = {}

max_x = -99999
max_y = -99999

for config_title, config_code in configurations.items():
    error_filename = path + "/" + filename + "_" + config_code + "_errors.txt"
    samplerates_filename = path + "/samplerates_" + config_code

    errors = []
    samplerates = []

    with open(error_filename, "r") as errors_file:
        for y in errors_file.read().split('\n'):
            try:
                val = float(y)
                max_y = max(max_y, val)
                errors.append(val)
            except:
                pass

    with open(samplerates_filename, "r") as samplerates_file:
        for y in samplerates_file.read().split('\n'):
            try:
                val = float(y)
                max_x = max(max_x, val)
                samplerates.append(val)
            except:
                pass

    unsorted = zip(samplerates, errors)
    coords[config_title] = sorted(unsorted, key=lambda tup: tup[0])

output_str = ""
output_str += "\\begin{tikzpicture}\n"
output_str += "\\begin{axis}[\n"
output_str += "title={" + scene_name + "},\n"
output_str += "xlabel={Total samples (\%)},\n"
output_str += "ylabel={Error (rmse)},\n"
output_str += "xmin=0, xmax=" + str(max_x) + ",\n"
output_str += "ymin=0, ymax=" + str(max_y) + ",\n"
output_str += "legend pos=north east,\n"
output_str += "ymajorgrids=true,\n"
output_str += "line width=2pt,\n"
output_str += "]\n\n"

for config_title, ccc in coords.items():
    output_str += "\\addplot[\n"
    output_str += "color=" + styles[config_title][0] + ",\n"
    output_str += "]\n"
    output_str += "coordinates {\n"
    for cc in ccc:
        output_str += "(" + str(cc[0]) + ", " + str(cc[1]) + ")"
    output_str += "\n"
    output_str += "};\n"
    output_str += "\\addlegendentry{" + config_title + "}\n"

output_str += "\\end{axis}\n"
output_str += "\\end{tikzpicture}\n"

o = open(output + ".tex", "w")
o.write(output_str)
o.close()
