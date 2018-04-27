#!/bin/bash

SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
CHR_REGION=$SHELL_FOLDER/../src/chr_region.txt

echo '''
****** Generate Chromatin Regions ******'''
echo '
>>>> The following 25Mb chromatin regions are generated as `./src/chr_region.txt` as default.'
echo '   > chromosome_id	start_position(Mb)	end_position(Mb)'

# Start positions (Mb) of chromatin segment (chromosomes 1 to 22)
# This array can be modified based on location of chromatin segments that are of interested
gSta_ARR=(20 20 20 20 20 20 20 60 80 88 60 45 25 25 35 50 30 20 34 35 20 20)

echo -n > $CHR_REGION
for ((i=0;i<22;i++));do
	let chrId="$i+1"
	gSta=${gSta_ARR[$i]}
	let gEnd="$gSta+25"
	echo "     $chrId			$gSta			$gEnd"
	echo "     $chrId			$gSta			$gEnd" >> $CHR_REGION
done
echo