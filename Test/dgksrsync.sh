#!bin/bash
set -e
if [ "$USER" = "Ark" ]; then
    DUGKSDIR=/mnt/d/HOME/Code/DUGKSAMReX
    echo "User is $USER"
elif [ "$USER" = "cc" ]; then
    DUGKSDIR=/mnt/c/Home/Code/DUGKSAMReX
    echo "User is $USER"
fi
rsync -avz --exclude '.*' ${DUGKSDIR} ~/
find ~/DUGKSAMReX -type d -exec chmod 755 {} +
find ~/DUGKSAMReX -type f ! -name "*.sh" ! -name "*.ex" -exec chmod 644 {} +
