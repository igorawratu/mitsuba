#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/classroom/timings"
rm "${source_dir}/classroom/samplerates"

#BOOL
for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.05 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_bool_05"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_bool_05"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_bool_05"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_bool_1"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_bool_1"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_bool_1"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_bool_2"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_bool_2"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_bool_2"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.4 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_bool_4"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_bool_4"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_bool_4"


#NON BOOL
for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.05 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nbool_05"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nbool_05"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nbool_05"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nbool_1"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nbool_1"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nbool_1"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nbool_2"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nbool_2"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nbool_2"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.4 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nbool_4"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nbool_4"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nbool_4"