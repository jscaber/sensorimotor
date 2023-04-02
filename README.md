# Sensorimotor 
This repository includes a number of analysis scripts used for the Manuscript "Cellular and axonal transport phenotypes due to the C9orf72 HRE in iPSC-motor and sensory neurons"
The code is shared under the GPL v3 licence, see LICENSE for full condition.
Written by Jakub Scaber

- live - folder containing live imaging analysis
  - WrapDifferenceTracker.py
    - Jython script for automating "Difference Tracker", an ImageJ plugin
    - Difference tracker consists of two plugins for ImageJ and is designed to track a large collection of faint objects 
    - Andrews S, Gilley J, Coleman MP. Difference Tracker: ImageJ plugins for fully automated analysis of multiple axonal transport parameters. J Neurosci Methods. 2010 Nov 30;193(2):281-7. doi: 10.1016/j.jneumeth.2010.09.007. Epub 2010 Sep 24. PMID: 20869987.
    - The script requires the user to specify a directory and will process all files with an extension specified by the user within the directory and its subdirectories
  - Figure_6_S6.R
    - R script for analysis of "Difference Tracker" output
    - Takes difference tracker output from above and cumulates data from subdirectories into a single file
    - Adds tracking of particles that don't run horizontally in grooves (not used)
    - Adds additional features such as directionality to output, and divides moving particles into 'dynamic stationary' and 'moving'
- rnaseq - folder containing rnaseq analysis
- Randomisation script (R Script)
- Pipeline for filtering 10X bam
