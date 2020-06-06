#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

chmod +x "${source_dir}/run_time_error_cbox.sh"
chmod +x "${source_dir}/run_time_error_breakfast-room.sh"
chmod +x "${source_dir}/run_time_error_classroom.sh"
chmod +x "${source_dir}/run_time_error_hairball.sh"
chmod +x "${source_dir}/run_time_error_living_room_3.sh"
chmod +x "${source_dir}/run_time_error_san_miguel.sh"
chmod +x "${source_dir}/run_time_error_sponza.sh"
chmod +x "${source_dir}/run_time_error_staircase.sh"
chmod +x "${source_dir}/run_time_error_kitchen.sh"
chmod +x "${source_dir}/run_time_error_bathroom.sh"
chmod +x "${source_dir}/run_time_error_staircase2.sh"
chmod +x "${source_dir}/run_time_error_living_room.sh"
chmod +x "${source_dir}/run_adaptive.sh"

#${source_dir}/run_time_error_cbox.sh
#${source_dir}/run_time_error_breakfast-room.sh
#${source_dir}/run_time_error_classroom.sh
#${source_dir}/run_time_error_hairball.sh
#${source_dir}/run_time_error_living_room_3.sh
#${source_dir}/run_time_error_san_miguel.sh
#${source_dir}/run_time_error_sponza.sh
#${source_dir}/run_time_error_staircase.sh
#${source_dir}/run_time_error_kitchen.sh
#${source_dir}/run_time_error_bathroom.sh
#${source_dir}/run_time_error_staircase2.sh
#${source_dir}/run_time_error_living_room.sh

python3 "${source_dir}/../generate_plots.py" "${source_dir}/cbox" "cbox" "Cornell box" "${source_dir}/cbox" "diffuse"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/breakfast-room/breakfast-room" "scene" "Breakfast room" "${source_dir}/breakfast-room" "diffuse"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/classroom" "scene" "Classroom" "${source_dir}/classroom" "diffuse"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/hairball" "hairball" "Hairball" "${source_dir}/hairball" "diffuse"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/living_room_3/living-room-3" "scene" "Living room" "${source_dir}/living-room" "diffuse"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/san_miguel" "san-miguel" "San-miguel" "${source_dir}/san-miguel" "diffuse"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/sponza" "sponza" "Sponza" "${source_dir}/sponza" "diffuse"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/staircase/staircase" "scene" "Staircase" "${source_dir}/staircase" "diffuse"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/kitchen" "scene" "Kitchen" "${source_dir}/kitchen" "glossy"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/bathroom/bathroom" "scene" "Bathroom" "${source_dir}/bathroom" "glossy"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/living_room/living-room" "scene" "Grey & White Room" "${source_dir}/living_room/living-room" "glossy"
python3 "${source_dir}/../generate_plots.py" "${source_dir}/staircase2/staircase2" "scene" "Modern Hall" "${source_dir}/staircase2/staircase2" "glossy"

chmod +x "${source_dir}/run_bool_sanmiguel.sh"
chmod +x "${source_dir}/run_is_sanmiguel.sh"
chmod +x "${source_dir}/run_bool_classroom.sh"
chmod +x "${source_dir}/run_is_classroom.sh"
chmod +x "${source_dir}/run_bool_kitchen.sh"
chmod +x "${source_dir}/run_is_kitchen.sh"

#${source_dir}/run_bool_sanmiguel.sh
# ${source_dir}/run_is_sanmiguel.sh
#${source_dir}/run_bool_classroom.sh
# ${source_dir}/run_is_classroom.sh
${source_dir}/run_bool_kitchen.sh
# ${source_dir}/run_is_kitchen.sh

python3 "${source_dir}/../generate_isb_plots.py" "${source_dir}/san_miguel" "san-miguel" "San-miguel" "${source_dir}/san-miguel_boolimport"
python3 "${source_dir}/../generate_isb_plots.py" "${source_dir}/classroom" "classroom" "Classroom" "${source_dir}/classroom_boolimport"
python3 "${source_dir}/../generate_isb_plots.py" "${source_dir}/kitchen" "kitchen" "Kitchen" "${source_dir}/kitchen_boolimport"

python3 "${source_dir}/../generatebaccplots.py" "${source_dir}/san_miguel" "san-miguel" "San Miguel" "${source_dir}/san-miguel_bacc"
python3 "${source_dir}/../generatebaccplots.py" "${source_dir}/classroom" "classroom" "Classroom" "${source_dir}/classroom_bacc"
python3 "${source_dir}/../generatebaccplots.py" "${source_dir}/kitchen" "kitchen" "Kitchen" "${source_dir}/kitchen_bacc"

# ${source_dir}/run_adaptive.sh

python3 "${source_dir}/../generate_adaptive_plots.py" "${source_dir}/classroom" "classroom" "Classroom" "${source_dir}/classroom_adaptive.tex"
