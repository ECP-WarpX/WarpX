#!/bin/sh
# This must be run from the local directory containing the Dockerfile!!!
# Requires an argument of the tag to provide, usually the git branch, and a
# suffix, usually the image type.
set -e

if [ -z "$1" ]
then
    echo "No argument supplied for the docker tag (usually current git branch)"
    exit 101
fi

DTAG=$1

# Per https://serverfault.com/a/382740 this accepts the empty string, but not
# no argument.
if [ -z "${2+set}" ]
then
    echo "No argument supplied for the docker suffix (image type)"
    exit 101
fi

SUFFIX=$2

FULLTAG="$1:$2"

if [ ! -d "../../../../WarpX" ]; then
    echo "No WarpX install found in ../../.."
    exit 111
elif [ ! -d "../../../../warpx-data" ]; then
    echo "No warpx-data install found in ../../.."
    exit 111
fi

# Do the docker build
cp .dockerignore ../../../../
$(aws ecr get-login --no-include-email --region us-west-2)
docker build --build-arg BUILDPLATFORM=$(uname -m) -f Dockerfile -t simteam/mewarpx-"$FULLTAG" ../../../../
docker tag simteam/mewarpx-"$FULLTAG" 167833485543.dkr.ecr.us-west-2.amazonaws.com/simteam/mewarpx-"$FULLTAG"
docker push 167833485543.dkr.ecr.us-west-2.amazonaws.com/simteam/mewarpx-"$FULLTAG"
