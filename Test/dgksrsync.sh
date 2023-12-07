#!bin/bash
set -e
DUGKSDIR=/mnt/d/HOME/Code/DUGKSAMReX
rsync -avz --exclude '.git/*' ~/code/DUGKSAMReX/Test/dgksrsync.sh ${DUGKSDIR}/Test/
rsync -avz --exclude '.*' ${DUGKSDIR} ~/code/
find ~/code/DUGKSAMReX -type d -exec chmod 755 {} +
find ~/code/DUGKSAMReX -type f ! -name "*.sh" ! -name "*.ex" -exec chmod 644 {} +
