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
    local MinL=$(awk "BEGIN {printf \"%.16f\",${ChL}/${ChNum}}")
    local U0=$(awk "BEGIN {printf \"%.16f\",${Pe}*${M}/${ChL}}")
    local Dt=$(awk "BEGIN {printf \"%.16f\",${CFL}*${ChL}/${ChNum}}")
    local EndStep=$(awk "BEGIN {print ${Period}*${ChL}*${ChNum}/${Pe}/${M}/${CFL}}")
#modify time
    ModifyLineSingleAppend 'time.CFL' ${CFL}
    ModifyLineSingleAppend 'time.advance_time_step' ${Dt}
#modify geometry
    ModifyLineTripleAppend 'amr.n_cell' ${ChNum}
    ModifyLineTripleAppend 'geometry.prob_hi' ${ChL}
    ModifyLineSingleAppend 'geometry.ch_num' ${ChNum}
    ModifyLineSingleAppend 'geometry.ch_length' ${ChL}
    ModifyLineSingleAppend 'geometry.min_length' ${MinL}
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
ModifyBenchmarkPseudoPotential() {
    local Mu=$(awk "BEGIN {printf \"%.16f\",${TauF}/3}")

    ModifyLineSingleAppend 'time.CFL' ${CFL}
    ModifyLineSingleAppend 'time.advance_time_step' ${CFL}
    ModifyLineSingleAppend 'base.T0' ${T}
    ModifyLineSingleAppend 'base.Mu0' ${Mu}
    ModifyLineSingleAppend 'base.Nu0' ${Mu}
    ModifyLineSingleAppend 'base.TauMomentum' ${TauF}
}
#phase field interface capture
# #set parameters
# Period=1
# Pe=512
# M=0.02
# Cn=0.0078125
# ChL=512
# ChNum=512
# CFL=0.8
# FluxScheme="DVM"
# CaseName="SmoothedDeformation"
# BenchmarkFile="./Benchmark/benchmark.${CaseName}"
# #loop over Peclet number or Mobility coefficients
# for Pe in 2048
# #for Pe in 512 1024 1280 1638.4 2048 2560
# # for Pe in 256 512 640 819.2 1024 1280
# # for M in 0.01 0.02 0.04 0.05 0.064 0.08 0.1
# # for ChNum in 128 192 256 384 512
# do
#     echo "Peclet : $Pe    Mobility : $M"
#     ChL=${ChNum}
#     wI=$(awk "BEGIN {print ${ChL}*${Cn}}")
#     ModifyBenchmark
#     L2File=L2_norm_Pe${Pe}_M${M}_wI${wI}_NL${ChNum}_Flux${FluxScheme}_${CaseName}.dat
#     TimeFile=elapsed_time_Pe${Pe}_M${M}_wI${wI}_NL${ChNum}_Flux${FluxScheme}_${CaseName}.dat
#     echo "L2-File : $L2File"
#     TargetFile=L2_Error_Pe${Pe}_M${M}_NL${ChNum}_wI${wI}_CFL${CFL}_Flux${FluxScheme}.${CaseName}
#     touch ${TargetFile}
#     echo "#####################BENCHMARK BEGIN##############################" >> ${TargetFile}
#     cat ${BenchmarkFile} >> ${TargetFile}
#     echo "#####################BENCHMARK END################################" >> ${TargetFile}
# #loop over coefficients
#     # for var in 0 16 24 32 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100 104 108 112 116 120
#     for var in 0
#     do
#         ModifyLineSingleAppend 'case.scheme_index' ${var}
#         make cleanOF > /dev/null
#         sh run.sh
#         # until [ -f "$L2File" ]
#         for iii in {1..400}
#         do
#             echo "The $iii-th sleep in $var"
#             sleep 30
#             [ $iii -eq 400 ] && exit
#             [ -f "$L2File" ] || continue
#             EndLine=$(awk '{print NR}' $L2File)
#             [ $EndLine -eq 1 ] && break
#             # LAST_MODIFY_TIMESTAMP=$(stat -c %Y $L2File)
#             # [ $TIME_STAMP -eq $LAST_MODIFY_TIMESTAMP ] && break
#             # TIME_STAMP=$LAST_MODIFY_TIMESTAMP
#         done
#         echo -n "-$var/64    " >> ${TargetFile}
#         awk 'END{print}' $L2File | tr -d '\n' >> ${TargetFile}
#         echo -n -e "\t" >> ${TargetFile}
#         awk 'END {print $NF}' $TimeFile >> ${TargetFile}
#         echo "write result at -$var/64"
#     done
# # mv ${TargetFile} ${TargetFile}.dat
# done

#pseudo-potential
#set parameters
ChL=1
ChNum=16
TauF=0.0003
CFL=0.0048
T=0.9
FluxScheme="CD"
CaseName="UnsteadyTGV"
BenchmarkFile="./Benchmark/${CaseName}.benchmark"
#loop over Peclet number or Mobility coefficients
for TauF in 0.3
do
    echo "TauF : $TauF"
    ModifyBenchmarkPseudoPotential
    VelFile=Level0_velocity_norm_U00_TauF${TauF}_T0${T}_NL${ChNum}_Flux${FluxScheme}_${CaseName}.dat
    L2File=L2_norm_U00_TauF${TauF}_T0${T}_NL${ChNum}_Flux${FluxScheme}_${CaseName}.dat
    echo "Vel-File : $VelFile"
    TargetFile=Velocity_Magnitude_TauF${TauF}_NL${ChNum}_T0${T}_CFL${CFL}_Flux${FluxScheme}.${CaseName}
    touch ${TargetFile}
    echo "#####################BENCHMARK BEGIN##############################" >> ${TargetFile}
    cat ${BenchmarkFile} >> ${TargetFile}
    echo "" >> ${TargetFile}
    echo "#####################BENCHMARK END################################" >> ${TargetFile}
#loop over coefficients
    # for var in 0 16 24 32 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100 104 108 112 116 120
    for var in 1.0
    do
        ModifyLineSingleAppend 'base.Alpha' ${var}
        make cleanOF > /dev/null
        sh run.sh
        # until [ -f "$L2File" ]
        for iii in {1..400}
        do
            echo "The $iii-th sleep in $var"
            sleep 30
            [ $iii -eq 400 ] && exit
            [ -f "$VelFile" ] || continue
            Endword=$(awk 'END {print $NF}' nohup.out)
            [[ $Endword == "finalized" ]] && break
            # LAST_MODIFY_TIMESTAMP=$(stat -c %Y $L2File)
            # [ $TIME_STAMP -eq $LAST_MODIFY_TIMESTAMP ] && break
            # TIME_STAMP=$LAST_MODIFY_TIMESTAMP
        done
        echo -n "$var    " >> ${TargetFile}
        awk 'END {printf $NF}' $VelFile >> ${TargetFile}
        echo -n "    " >> ${TargetFile}
        awk 'END {print $NF}' $L2File >> ${TargetFile}
        echo "write result at Alpha=$var"
    done
# mv ${TargetFile} ${TargetFile}.dat
done
