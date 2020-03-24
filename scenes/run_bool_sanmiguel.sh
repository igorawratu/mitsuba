#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/san_miguel/timings"
rm "${source_dir}/san_miguel/samplerates"

#BOOL
for i in {1..5}
do
    mitsuba "${source_dir}/san_miguel/san-miguel.xml" -Dcstrat=3 -Dsamplerate=0.05 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/san_miguel/san-miguel.exr" "${source_dir}/san_miguel/groundtruth.exr" "${source_dir}/san_miguel/san-miguel_bool_05"
done
mv "${source_dir}/san_miguel/timings" "${source_dir}/san_miguel/timings_bool_05"
mv "${source_dir}/san_miguel/samplerates" "${source_dir}/san_miguel/samplerates_bool_05"

for i in {1..5}
do
    mitsuba "${source_dir}/san_miguel/san-miguel.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/san_miguel/san-miguel.exr" "${source_dir}/san_miguel/groundtruth.exr" "${source_dir}/san_miguel/san-miguel_bool_1"
done
mv "${source_dir}/san_miguel/timings" "${source_dir}/san_miguel/timings_bool_1"
mv "${source_dir}/san_miguel/samplerates" "${source_dir}/san_miguel/samplerates_bool_1"

for i in {1..5}
do
    mitsuba "${source_dir}/san_miguel/san-miguel.xml" -Dcstrat=3 -Dsamplerate=0.15 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/san_miguel/san-miguel.exr" "${source_dir}/san_miguel/groundtruth.exr" "${source_dir}/san_miguel/san-miguel_bool_15"
done
mv "${source_dir}/san_miguel/timings" "${source_dir}/san_miguel/timings_bool_15"
mv "${source_dir}/san_miguel/samplerates" "${source_dir}/san_miguel/samplerates_bool_15"

for i in {1..5}
do
    mitsuba "${source_dir}/san_miguel/san-miguel.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/san_miguel/san-miguel.exr" "${source_dir}/san_miguel/groundtruth.exr" "${source_dir}/san_miguel/san-miguel_bool_2"
done
mv "${source_dir}/san_miguel/timings" "${source_dir}/san_miguel/timings_bool_2"
mv "${source_dir}/san_miguel/samplerates" "${source_dir}/san_miguel/samplerates_bool_2"


#NON BOOL
for i in {1..5}
do
    mitsuba "${source_dir}/san_miguel/san-miguel.xml" -Dcstrat=3 -Dsamplerate=0.05 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/san_miguel/san-miguel.exr" "${source_dir}/san_miguel/groundtruth.exr" "${source_dir}/san_miguel/san-miguel_nbool_05"
done
mv "${source_dir}/san_miguel/timings" "${source_dir}/san_miguel/timings_nbool_05"
mv "${source_dir}/san_miguel/samplerates" "${source_dir}/san_miguel/samplerates_nbool_05"

for i in {1..5}
do
    mitsuba "${source_dir}/san_miguel/san-miguel.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/san_miguel/san-miguel.exr" "${source_dir}/san_miguel/groundtruth.exr" "${source_dir}/san_miguel/san-miguel_nbool_1"
done
mv "${source_dir}/san_miguel/timings" "${source_dir}/san_miguel/timings_nbool_1"
mv "${source_dir}/san_miguel/samplerates" "${source_dir}/san_miguel/samplerates_nbool_1"

for i in {1..5}
do
    mitsuba "${source_dir}/san_miguel/san-miguel.xml" -Dcstrat=3 -Dsamplerate=0.15 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/san_miguel/san-miguel.exr" "${source_dir}/san_miguel/groundtruth.exr" "${source_dir}/san_miguel/san-miguel_nbool_15"
done
mv "${source_dir}/san_miguel/timings" "${source_dir}/san_miguel/timings_nbool_15"
mv "${source_dir}/san_miguel/samplerates" "${source_dir}/san_miguel/samplerates_nbool_15"

for i in {1..5}
do
    mitsuba "${source_dir}/san_miguel/san-miguel.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/san_miguel/san-miguel.exr" "${source_dir}/san_miguel/groundtruth.exr" "${source_dir}/san_miguel/san-miguel_nbool_2"
done
mv "${source_dir}/san_miguel/timings" "${source_dir}/san_miguel/timings_nbool_2"
mv "${source_dir}/san_miguel/samplerates" "${source_dir}/san_miguel/samplerates_nbool_2"
