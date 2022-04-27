# Sensorimotor 
This repository includes a number of analysis scripts used for the Manuscript "iPSC-derived motor and sensory neurons have distinct molecular profiles but share axonal vulnerability due to the C9orf72 hexonucleotide repeat expansion"
The code is shared under the GPL v3 licence, see LICENSE for full condition.
Written by Jakub Scaber

1. Jython script for automating "Difference Tracker", an ImageJ plugin
   - WrapDifferenceTracker.py
   - Difference tracker consists of two plugins for ImageJ and is designed to track a large collection of faint objects 
   - Andrews S, Gilley J, Coleman MP. Difference Tracker: ImageJ plugins for fully automated analysis of multiple axonal transport parameters. J Neurosci Methods. 2010 Nov 30;193(2):281-7. doi: 10.1016/j.jneumeth.2010.09.007. Epub 2010 Sep 24. PMID: 20869987.
   - The script requires the user to specify a directory and will process all files with an extension specified by the user within the directory and its subdirectories
2. R script for analysis of "Difference Tracker" output
   - NameOfScript.R
   - Takes difference tracker output from above and cumulates data from subdirectories into a single file
   - Adds tracking of particles that don't run horizontally in grooves
   - Adds additional features such as directionality to output, and divides moving particles into 'dynamic stationary' and 'moving'
3. Randomisation script (R Script)
   - here
