#!/bin/bash

# This script handles running MPI/WarpX in docker containers on AWS, and the
# associated copy and monitoring operations that accompany it. It detects runs
# that don't output to stdout for $TIMEOUT seconds and will terminate at that
# point. It will move files either as the run is ongoing, or just at the end,
# back to S3. Dump files are never moved or copied to S3; they will be copied
# from S3 if they're already there.

# S3 can sometimes make directories with no name ("") if two slashes are
# present in sequence
# https://stackoverflow.com/questions/1848415/remove-slash-from-the-end-of-a-variable
DIRNAME=${1%/}
BUCKET=$2
SCRIPTNAME='run_script.sh'

# MVFILES can be an environment variable set to True to move files as
# generated.
echo "Move diagnostic files as generated:" $MVFILES

# TIMEOUT can be an environment variable. If no output to stdout has been made
# in the previous TIMEOUT seconds, this script terminates.
# TIMEOUT defaults to 3h = 10800 s
# https://stackoverflow.com/questions/3601515/how-to-check-if-a-variable-is-set-in-bash
if [ -z ${TIMEOUT:x} ] || ! (( "$TIMEOUT" > 0 )); then
    echo TIMEOUT set to default of 10800 sec
    TIMEOUT=10800
else
    echo TIMEOUT set to "$TIMEOUT"
fi

# Use EFS directory as working directory.
echo Make directory ${DIRNAME} in EFS and cd there
mkdir -p /efs/WarpX_simulation_runs/${DIRNAME}
cd /efs/WarpX_simulation_runs/${DIRNAME}

# Copy files from S3
aws s3 cp --recursive s3://${BUCKET}/${DIRNAME} ./ --exclude "*" \
    --include "run*" --include "*.dump" --include "std*" \
    --include "*runinfo.dpkl"
ls

# Move diagnostic files that are done
function checkptmove {
    # Check file is not open. Since they shouldn't be re-opened, this should be
    # fine.
    if ! lsof $1 > /dev/null
    then
        # http://www.linuxjournal.com/content/normalizing-path-names-bash
        # Remove all /./ and // sequences.
        local path=${DIRNAME}/"$1"
        path=${path//\/.\//\/}
        path=${path//\/\//\/}
        path=${path/diags_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]/diags}

        # Remove dir/.. sequences.
        while [[ $path =~ ([^/][^/]*/\.\./) ]]
        do
            path=${path/${BASH_REMATCH[0]}/}
        done
        echo $path
        aws s3 mv $1 s3://${BUCKET}/${path}
    fi
}

# Every 5 minutes, move (to keep space free, if $MVFILES is True) or sync
# per-timestep diagnostic files as they are done being written. Always sync
# other files (which may still be opened!)
function checkpoint {
while true
do
    sleep 300
    if [[ "$MVFILES" == "True" ]]
    then
        find . -name "*Header" | while read file; do checkptmove "$file"; done
        find . -name "Cell_D_*[0-9][0-9][0-9][0-9][0-9]" | while read file; do checkptmove "$file"; done
        find . -name "warpx_job_info" | while read file; do checkptmove "$file"; done
        find . -name "*.h5" | while read file; do checkptmove "$file"; done
        find . -name "*.npy" | while read file; do checkptmove "$file"; done
        find . -name "*.dpkl" | while read file; do checkptmove "$file"; done
        find . -name "*.png" | while read file; do checkptmove "$file"; done
        find . -name "*.txt" | while read file; do checkptmove "$file"; done
        find . -name "*.json" | while read file; do checkptmove "$file"; done
        aws s3 sync ./ s3://${BUCKET}/${DIRNAME} --exclude "*" --include "std*" --include "*.json"
    else
        for d in ./diags_* ; do
            aws s3 sync $d s3://${BUCKET}/${DIRNAME}/diags/
        done
        aws s3 sync ./diags/ s3://${BUCKET}/${DIRNAME}/diags/
        aws s3 sync ./ s3://${BUCKET}/${DIRNAME} --exclude "*" --include "std*" --include "*.json"
    fi

    # In addition, check if no output to stdout has been made for last TIMEOUT
    # seconds, if so terminate this run.
    # https://stackoverflow.com/questions/28337961/find-out-if-file-has-been-modified-within-the-last-2-minutes
    CURTIME=$(date +%s)
    FILETIME=$(stat stdout.out -c %Y)
    TIMEDIFF=$(expr $CURTIME - $FILETIME)
    echo stdout was written to $TIMEDIFF seconds ago

    if [ $TIMEDIFF -gt $TIMEOUT ]; then
        echo Timeout - no output to stdout after ${TIMEOUT} seconds
        exit 33
    fi

    # Finally, 'touch' each file in the run directory so files for runs that
    # last many days won't be deleted from EFS before the run is completed
    find . -type f -exec touch -a -c {} +
done
}
# https://stackoverflow.com/questions/21465297/tee-stdout-and-stderr-to-separate-files-while-retaining-them-on-their-respective
{ { bash $SCRIPTNAME; } > >( tee -a stdout.out ); } \
                       2> >( tee -a stderr.err >&2 ) &

MAINID=$!
checkpoint &
CHECKID=$!
wait -n $MAINID $CHECKID
EXITFLAG=$?
# kill, but silently, as one will error no matter what. It looks like killing
# $MAINID actually isn't enough -- child processes still run -- but it
# shouldn't matter as the docker container should exit when this script exits,
# I think.
kill $CHECKID &> /dev/null
kill $MAINID &> /dev/null

for d in ./diags_* ; do
    aws s3 mv --recursive $d s3://${BUCKET}/${DIRNAME}/diags/
done
aws s3 mv --recursive ./diags/ s3://${BUCKET}/${DIRNAME}/diags/
aws s3 mv --recursive ./ s3://${BUCKET}/${DIRNAME} --exclude "*" --include "std*" --include "*.json"

echo Exit with status $EXITFLAG
exit $EXITFLAG
