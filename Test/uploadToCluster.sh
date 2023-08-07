#!bash
set -e
workingdir=$(pwd)
parentdir=${workingdir%/*}
localdir=${parentdir##*/}
echo "folder name : $localdir"
#
if [ "$1" = "hn" ];
then
	hostIP=10.168.36.131
	rootdir=home
elif [ "$1" = "amd" ];
then
	hostIP=10.168.47.197
	hostName=amd
	rootdir=home
else
	echo "unknown host name"
	exit 1
fi
#
if [ "$2" = "HO" -o "$2" = "ho" ];
then
	remotedir="HODUGKSAMReX"
elif [ "$2" = "DVM" -o "$2" = "dvm" ];
then
	remotedir="DVMDUGKSAMReX"
elif [ "$2" = "" ];
then
	remotedir="DUGKSAMReX"
else
	echo "unknown host name"
	exit 1
fi	
echo "working directory : ${workingdir}"
echo "remote host name : $1"
echo "remote directory : ${remotedir}"
scp GNUmakefile yangzr@${hostName}:/${rootdir}/yangzr/code/${remotedir}/Test
scp JOB.sh yangzr@${hostName}:/${rootdir}/yangzr/code/${remotedir}/Test
scp run.sh yangzr@${hostName}:/${rootdir}/yangzr/code/${remotedir}/Test
echo "uploading Test folder to $1 ${remotedir}"
scp -r ./Benchmark/* yangzr@${hostName}:/${rootdir}/yangzr/code/${remotedir}/Test/Benchmark >/dev/null
echo "uploading Ark folder to $1 ${remotedir}"
scp -r ../Ark/* yangzr@${hostName}:/${rootdir}/yangzr/code/${remotedir}/Ark >/dev/null
echo "upload done!"
