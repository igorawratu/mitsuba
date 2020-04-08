#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/kitchen/timings"
rm "${source_dir}/kitchen/samplerates"

#IMPORTANCE SAMPLING
for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.05 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -iet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_import_05"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_import_05"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_import_05"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -iet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_import_1"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_import_1"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_import_1"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.15 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -iet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_import_15"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_import_15"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_import_15"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -iet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_import_2"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_import_2"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_import_2"


#NO IMPORTANCE SAMPLING
for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.05 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=0 -Dbv=0 -Dvp=0 -iet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nimport_05"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nimport_05"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nimport_05"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=0 -Dbv=0 -Dvp=0 -iet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nimport_1"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nimport_1"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nimport_1"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.15 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=0 -Dbv=0 -Dvp=0 -iet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nimport_15"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nimport_15"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nimport_15"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=0 -Dbv=0 -Dvp=0 -iet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nimport_2"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nimport_2"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nimport_2"
