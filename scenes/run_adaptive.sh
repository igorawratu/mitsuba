#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/classroom/timings"
rm "${source_dir}/classroom/samplerates"

#ADAPTIVE
for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0.25
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_adaptive"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_adaptive"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_adaptive"

#0.1 STATIC
for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_static_01"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_static_01"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_static_01"

#0.4 STATIC
for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.4 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_static_04"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_static_04"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_static_04"