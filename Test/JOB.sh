#!/bin/bash
# declare -a coefficients=(0 12 24 32 34 35 36 37 38 40 48)
set -e
# set -x
#functions
ModifyLineTripleAppend() {
    local line=$(sed -n '/'$1'/=' ${BenchmarkFile})
    test -z $line && exit -1
    sed -i 's/^\('$1'[[:blank:]]*=[[:blank:]]*\).*/\1'$2'  '$2'  '$2'/' ${BenchmarkFile}
}
ModifyLineSingleAppend() {
    local line=$(sed -n '/'$1'/=' ${BenchmarkFile})
    test -z $line && exit -1
    sed -i 's/^\('$1'[[:blank:]]*=[[:blank:]]*\).*/\1'$2'/' ${BenchmarkFile}
} 
ModifyBenchmark() {
    local TauM=$(awk "BEGIN {print ${M}*3}")
    local U0=$(awk "BEGIN {printf \"%.16f\",${Pe}*${M}/${ChL}}")
    local EndStep=$(awk "BEGIN {print ${ChL}*${ChL}/${Pe}/${M}/${CFL}}")
#modify time
    ModifyLineSingleAppend 'time.CFL' ${CFL}
    ModifyLineSingleAppend 'time.advance_time_step' ${CFL}	
#modify geometry
    ModifyLineTripleAppend 'amr.n_cell' ${ChL}
    ModifyLineTripleAppend 'geometry.prob_hi' ${ChL}
    ModifyLineSingleAppend 'geometry.ch_num' ${ChL}
    ModifyLineSingleAppend 'geometry.ch_length' ${ChL}
#modify base
    ModifyLineSingleAppend 'base.U0' ${U0}
#modify phase field model
    ModifyLineSingleAppend 'PF.widthI' ${wI}
    ModifyLineSingleAppend 'PF.mobility_phi' ${M}
    ModifyLineSingleAppend 'PF.tau_phi' ${TauM}
#modify control
    ModifyLineSingleAppend 'control.end_step' ${EndStep}
    ModifyLineSingleAppend 'control.conv_interval' ${EndStep}
}
#set parameters
Pe=256
M=0.02
Cn=0.015625
ChL=256
CFL=0.5
FluxScheme="UW"
CaseName="TransDrop"
BenchmarkFile="benchmark.transDrop"
#loop over Peclet number or Mobility coefficients
for Pe in 256 512 640 819.2 1024 1280
# for M in 0.01 0.02 0.04 0.05 0.064 0.08 0.1
# for ChL in 128 192 256 384 512
do
    echo "Peclet : $Pe    Mobility : $M"
    wI=$(awk "BEGIN {print ${ChL}*${Cn}}")
    ModifyBenchmark
    L2File=L2_norm_Pe${Pe}_M${M}_wI${wI}_NL${ChL}_${CaseName}.dat
    TargetFile=L2_Error_Pe${Pe}_M${M}_NL${ChL}_wI${wI}_CFL${CFL}\
_FluxScheme${FluxScheme}_${CaseName}_FreeParameterEffects
    touch ${TargetFile}
    echo "#####################BENCHMARK BEGIN##############################" >> ${TargetFile}
    cat ${BenchmarkFile} >> ${TargetFile}
    echo "#####################BENCHMARK END################################" >> ${TargetFile}
#loop over coefficients
    for var in -24 24 32 34 35 36 37 38 40 48
    #for var in 64 68 72 74 76 78 80
    do
        ModifyLineSingleAppend 'case.scheme_index' ${var}
        make cleanOF > /dev/null
        sh run.sh
        # until [ -f "$L2File" ]
        for iii in {1..50}  
        do
            echo "The $iii-th sleep in $var"
            sleep 30
            [ $iii -eq 50 ] && exit
            [ -f "$L2File" ] || continue
            LAST_MODIFY_TIMESTAMP=$(stat -c %Y $L2File) 
            [ $TIME_STAMP -eq $LAST_MODIFY_TIMESTAMP ] && break
            TIME_STAMP=$LAST_MODIFY_TIMESTAMP
        done
        echo -n "-$var/48    " >> ${TargetFile}
        awk 'END{print}' $L2File >> ${TargetFile}
        echo "write result at -$var/48"
    done
# mv ${TargetFile} ${TargetFile}.dat
done
