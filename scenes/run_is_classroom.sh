#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/classroom/timings"
rm "${source_dir}/classroom/samplerates"

#IMPORTANCE SAMPLING
for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.05 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_import_05"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_import_05"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_import_05"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_import_1"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_import_1"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_import_1"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.15 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_import_15"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_import_15"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_import_15"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_import_2"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_import_2"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_import_2"


#NO IMPORTANCE SAMPLING
for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.05 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=0 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nimport_05"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nimport_05"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nimport_05"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=0 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nimport_1"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nimport_1"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nimport_1"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.15 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=0 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nimport_15"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nimport_15"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nimport_15"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=0 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nimport_2"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nimport_2"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nimport_2"