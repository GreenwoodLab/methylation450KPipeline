#! /bin/bash
# Description: pipeline 450K methylation
# Author : Voisin Greg Lady Davis institute
# Date: June 2015
#
# Adaptation for CBRAIN: Natacha Beck MCIN
# Date: July 2015


# NOTE REMOVED OUTPUT ! PUT IN RUBY PART
# input     = $1
# output    = $2
# threshold = $3

# Removed sslVerify for git
# git config --global http.sslVerify false
# git clone https://github.com/GreenwoodLab/methylation450KPipeline.git $PIPELINE_450K

# In this case it's the second pass
if [ -d "$1/IDAT" ]
then
  rsync -a -v $1 $2
else
  # Create output directory
  mkdir -p $2
  # Set the path to analyze
  Rscript $PIPELINE_450K/initialize_analysis.R $2
  # Create the output and copy *.idat and samplesheet data in it
  cp $1/*.* $2/IDAT
fi

# Run the script for the processing of *.idat file (before QC).
Rscript $PIPELINE_450K/processingIdatData_before_QC.R $2 distant
# Rscript $PIPELINE_450K/processingQCData.R $2

# Run the script for the QC step
cp $PIPELINE_450K/QCReport.Rmd $2/QC_REPORT/
Rscript $PIPELINE_450K/processingQCData.R $2 distant

# Run the processingIdatData_after_QC step
Rscript $PIPELINE_450K/processingIdatData_after_QC.R $2 distant $3
