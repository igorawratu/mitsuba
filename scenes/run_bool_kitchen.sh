#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/kitchen/timings"
rm "${source_dir}/kitchen/samplerates"

# #BOOL
# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_1_1000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_1_1000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_1_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.15 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_2_1000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_15_1000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_15_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_2_1000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_2_1000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_2_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=2000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_1_2000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_1_2000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_1_2000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_1_3000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_1_3000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_1_3000"


# #NON BOOL
# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0 -Dge=1
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_1_1000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_1_1000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_1_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.15 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0 -Dge=1
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_2_1000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_15_1000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_15_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0 -Dge=1
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_2_1000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_2_1000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_2_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=2000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0 -Dge=1
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_1_2000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_1_2000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_1_2000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0 -Dge=1
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_1_3000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_1_3000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_1_3000"

#GE
# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_1_1000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_1_1000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_1_1000"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.15 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_2_1000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_15_1000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_15_1000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.2 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_2_1000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_2_1000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_2_1000"

for i in {1..5}
do
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=2000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_1_2000"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_1_2000"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_1_2000"

# for i in {1..5}
# do
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=0.1 -Dcompcstrat=mdlc -Dcps=3000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolge_1_3000"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolge_1_3000"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolge_1_3000"




# #BOOL
# for i in {1..25}
# do
#     sr=$(bc <<< "0.01 * $i")
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=$sr -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolacc"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolacc"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolacc"

# #NBOOL
# for i in {1..25}
# do
#     sr=$(bc <<< "0.01 * $i")
#     mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=$sr -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=0 -Dvp=0 -Diet=0 -Dge=0
#     python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_nboolacc"
# done
# mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_nboolacc"
# mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_nboolacc"

#GE
for i in {1..25}
do
    sr=$(bc <<< "0.01 * $i")
    mitsuba "${source_dir}/kitchen/scene.xml" -Dcstrat=3 -Dsamplerate=$sr -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000 -Dis=1 -Dbv=1 -Dvp=0 -Diet=0 -Dge=1
    python3 "${source_dir}/../compare_err.py" "${source_dir}/kitchen/scene.exr" "${source_dir}/kitchen/groundtruth.exr" "${source_dir}/kitchen/kitchen_boolgeacc"
done
mv "${source_dir}/kitchen/timings" "${source_dir}/kitchen/timings_boolgeacc"
mv "${source_dir}/kitchen/samplerates" "${source_dir}/kitchen/samplerates_boolgeacc"