# Run multiqc after performing fastqc for the trimmed or adapted data.

#!/bin/bash

echo -e "Moving output messages"
#Create a folder for the output messages.
mkdir -p outputMessages && mv errfile* outfile* $_

echo -e "Launching multiqc"
module add UHTS/Analysis/MultiQC/1.8
multiqc .
echo -e "Done"
