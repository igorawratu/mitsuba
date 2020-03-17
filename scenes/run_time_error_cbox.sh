#!/bin/bash
set -x
source_path="${BASH_SOURCE[0]}"
source_dir=$(dirname "${source_path}")

rm "${source_dir}/cbox/cbox_"*
rm "${source_dir}/cbox/timings"*
rm "${source_dir}/cbox/samplerates"*


#MDLC with matrix sep
for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=500 -Dsps=1 -Dspp=4 -Dslice=1000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_amr_500"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_amr_500"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_amr_500"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_amr_1000"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_amr_1000"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_amr_1000"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=2000 -Dsps=1 -Dspp=4 -Dslice=1000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_amr_2000"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_amr_2000"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_amr_2000"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=4000 -Dsps=1 -Dspp=4 -Dslice=1000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_amr_4000"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_amr_4000"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_amr_4000"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=0.025 -Dcompcstrat=mdlc -Dcps=8000 -Dsps=1 -Dspp=4 -Dslice=1000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_amr_8000"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_amr_8000"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_amr_8000"


#MDLC
for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=250 -Dsps=1 -Dspp=4 -Dslice=1000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_mdlc_250"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_mdlc_250"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_mdlc_250"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=500 -Dsps=1 -Dspp=4 -Dslice=1000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_mdlc_500"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_mdlc_500"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_mdlc_500"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=1000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_mdlc_1000"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_mdlc_1000"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_mdlc_1000"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=2000 -Dsps=1 -Dspp=4 -Dslice=1000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_mdlc_2000"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_mdlc_2000"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_mdlc_2000"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=mdlc -Dcps=4000 -Dsps=1 -Dspp=4 -Dslice=1000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_mdlc_4000"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_mdlc_4000"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_mdlc_4000"


#LS
for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=250 -Dsps=1 -Dspp=4 -Dslice=50000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_ls_250"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_ls_250"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_ls_250"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=500 -Dsps=1 -Dspp=4 -Dslice=50000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_ls_500"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_ls_500"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_ls_500"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=1000 -Dsps=1 -Dspp=4 -Dslice=50000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_ls_1000"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_ls_1000"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_ls_1000"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=2000 -Dsps=1 -Dspp=4 -Dslice=50000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_ls_2000"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_ls_2000"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_ls_2000"

for i in {1..5}
do
	mitsuba "${source_dir}/cbox/cbox.xml" -Dcstrat=3 -Dsamplerate=1 -Dcompcstrat=ls -Dcps=4000 -Dsps=1 -Dspp=4 -Dslice=50000
	python3 "${source_dir}/../compare_err.py" "${source_dir}/cbox/cbox.exr" "${source_dir}/cbox/groundtruth.exr" "${source_dir}/cbox/cbox_ls_4000"
done
mv "${source_dir}/cbox/timings" "${source_dir}/cbox/timings_ls_4000"
mv "${source_dir}/cbox/samplerates" "${source_dir}/cbox/samplerates_ls_4000"