#!bin/bash
set -e
DUGKSDIR=/mnt/d/home/Code/DUGKSAMReX
rsync -avz --exclude '.git/*' ~/DUGKSAMReX/Test/dgksrsync.sh ${DUGKSDIR}/Test/
rsync -avz --exclude '.*' ${DUGKSDIR} ~/
find ~/DUGKSAMReX -type d -exec chmod 755 {} +
find ~/DUGKSAMReX -type f ! -name "*.sh" ! -name "*.ex" -exec chmod 644 {} +
