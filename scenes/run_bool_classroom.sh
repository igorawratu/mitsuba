#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/classroom/timings"
rm "${source_dir}/classroom/samplerates"

# #BOOL
# for i in {1..5}
# do
#     mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_bool_1_1000"
# done
# mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_bool_1_1000"
# mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_bool_1_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_bool_2_1000"
# done
# mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_bool_2_1000"
# mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_bool_2_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_bool_1_3000"
# done
# mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_bool_1_3000"
# mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_bool_1_3000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_bool_2_3000"
# done
# mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_bool_2_3000"
# mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_bool_2_3000"


# #NON BOOL
# for i in {1..5}
# do
#     mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nbool_1_1000"
# done
# mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nbool_1_1000"
# mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nbool_1_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nbool_2_1000"
# done
# mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nbool_2_1000"
# mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nbool_2_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nbool_1_3000"
# done
# mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nbool_1_3000"
# mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nbool_1_3000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nbool_2_3000"
# done
# mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nbool_2_3000"
# mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nbool_2_3000"


#GE
for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_gebool_1_1000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_gebool_1_1000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_gebool_1_1000"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_gebool_2_1000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_gebool_2_1000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_gebool_2_1000"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_gebool_1_3000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_gebool_1_3000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_gebool_1_3000"

for i in {1..5}
do
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_gebool_2_3000"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_gebool_2_3000"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_gebool_2_3000"


#BOOL
for i in {1..50}
do
    sr=$(bc <<< "0.01 * $i")
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=$sr -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_boolacc"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_boolacc"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_boolacc"

#NBOOL
for i in {1..50}
do
    sr=$(bc <<< "0.01 * $i")
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=$sr -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0 -Dge=0
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_nboolacc"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_nboolacc"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_nboolacc"

#GE
for i in {1..50}
do
    sr=$(bc <<< "0.01 * $i")
    mitsuba "${source_dir}/classroom/scene.xml" -Dcstrat=3 -Dsamplerate=$sr -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
    python3 "${source_dir}/../compare_err.py" "${source_dir}/classroom/scene.exr" "${source_dir}/classroom/groundtruth.exr" "${source_dir}/classroom/classroom_boolgeacc"
done
mv "${source_dir}/classroom/timings" "${source_dir}/classroom/timings_boolgeacc"
mv "${source_dir}/classroom/samplerates" "${source_dir}/classroom/samplerates_boolgeacc"