#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/kitchen/timings"
rm "${source_dir}/kitchen/samplerates"

#BOOL
for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.05 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_bool_05"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_bool_05"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_bool_05"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_bool_1"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_bool_1"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_bool_1"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_bool_2"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_bool_2"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_bool_2"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.4 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_bool_4"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_bool_4"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_bool_4"


#NON BOOL
for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.05 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nbool_05"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nbool_05"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nbool_05"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nbool_1"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nbool_1"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nbool_1"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nbool_2"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nbool_2"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nbool_2"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.4 -Dcompcstrat=mdlc -Dcps=5000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nbool_4"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nbool_4"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nbool_4"