#!/bin/sh

BW=$1
REGIONS=$2

cat $REGIONS | while read chr start end n; do \
    echo $chr $start $end $n mean `bigWigSummary -type=mean $BW $chr $start $end $n`;
    echo $chr $start $end $n min `bigWigSummary -type=min $BW $chr $start $end $n`;
    echo $chr $start $end $n max `bigWigSummary -type=max $BW $chr $start $end $n`;
    echo $chr $start $end $n std `bigWigSummary -type=std $BW $chr $start $end $n`;
done;
