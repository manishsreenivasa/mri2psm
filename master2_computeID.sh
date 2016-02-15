#!/bin/bash
dataPath='./'
for file in $dataPath/sampleData/c3d/*.c3d
do
	nameis=$(basename $file .c3d)
	echo 'Processing GSM inverse kinematics with c3d file, ' $nameis
	fit_motion $dataPath/model/data_sampleGSM.lua $file
	mv animation.csv ./sampleResults/ik_gsm_$nameis.csv
	mv fitting_log.csv ./sampleResults/ik_gsm_$nameis.fit
	echo 'Processing GSM inverse dynamics with c3d file, ' $nameis
	$dataPath/programs/2_dynamicsComputations/build/inverse_dynamics $dataPath/sampleResults/ik_gsm_$nameis.csv $dataPath/model/data_sampleGSM.lua $file 
	mv id_res.txt ./sampleResults/id_gsm_$nameis.txt
	echo 'Processing PSM inverse kinematics with c3d file, ' $nameis
	fit_motion $dataPath/model/data_samplePSM.lua $file
	mv animation.csv ./sampleResults/ik_psm_$nameis.csv
	mv fitting_log.csv ./sampleResults/ik_psm_$nameis.fit
	echo 'Processing PSM inverse dynamics with c3d file, ' $nameis
	$dataPath/programs/2_dynamicsComputations/build/inverse_dynamics $dataPath/sampleResults/ik_psm_$nameis.csv $dataPath/model/data_samplePSM.lua $file 
	mv id_res.txt ./sampleResults/id_psm_$nameis.txt
done
