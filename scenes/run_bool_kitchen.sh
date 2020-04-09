#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/kitchen/timings"
rm "${source_dir}/kitchen/samplerates"

#BOOL
for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_bool_1_1000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_bool_1_1000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_bool_1_1000"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_bool_2_1000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_bool_2_1000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_bool_2_1000"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_bool_1_3000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_bool_1_3000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_bool_1_3000"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_bool_2_3000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_bool_2_3000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_bool_2_3000"


#NON BOOL
for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nbool_1_1000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nbool_1_1000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nbool_1_1000"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nbool_2_1000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nbool_2_1000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nbool_2_1000"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nbool_1_3000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nbool_1_3000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nbool_1_3000"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nbool_2_3000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nbool_2_3000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nbool_2_3000"