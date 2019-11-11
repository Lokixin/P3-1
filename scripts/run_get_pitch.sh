#!/bin/bash

# Put here the program (maybe with path)
GETF0="get_pitch"

for fwav in pitch_db/train/*.wav; do
    echo "$fwav ----"
    ff0=${fwav/.wav/.f0}
    $GETF0 $fwav $ff0 > /dev/null || (echo "Error in $GETF0 $fwav $ff0"; exit 1)
done

exit 0
