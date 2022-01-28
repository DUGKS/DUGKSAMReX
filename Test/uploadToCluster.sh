#!bash
set -e
workingdir=$(pwd)
parentdir=${workingdir%/*}
localdir=${parentdir##*/}
echo "folder name : $localdir"
#
if [ "$1" = "hn" ];
then
	hostIP=10.168.17.85
	rootdir=home
elif [ "$1" = "mh" ];
then
	hostIP=10.168.6.48
	rootdir=data
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
scp JOB.sh yangzr@${hostIP}:/${rootdir}/yangzr/Code/${remotedir}/Test >/dev/null
echo "uploading Test folder to $1 ${remotedir}"
scp -r ./Benchmark/* yangzr@${hostIP}:/${rootdir}/yangzr/Code/${remotedir}/Test/Benchmark >/dev/null
echo "uploading Ark folder to $1 ${remotedir}"
scp -r ../Ark/* yangzr@${hostIP}:/${rootdir}/yangzr/Code/${remotedir}/Ark >/dev/null
echo "upload done!"
