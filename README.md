# nuskybgd:
# An IDL module for producing simulated background for NuSTAR.


## Overview

nuskybgd is code for simulating the NuSTAR background the Cosmix X-ray
Background (CXB). It simulates the X-rays that are focued throught the
optic as well the "stray light" from the CXB that leaks in around the
optics bench and through the aperture stop.

In /docs there is a detailed documentation and walk through. Please
read this first before attempting to use this software.

## Reference:

This code was originally written by Dan Wik, and was originally described in the
appendix of the paper here:

http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1403.2722

If you use nuskbygd, please reference this paper.


## Warnings:

This code is offered "as is" for the user. Users should be warned that
using nuskbygd occasionally requires some tweaking of the scripts to
work for specific cases.

Please also note that nuskbygd relies (heavily) on the AstroLib, which
can be checked out here:

git clone https://github.com/wlandsman/IDLAstro.git astrolib-idl


Installation
------------

1. Clone our project from github or download a release tarball

    $ git clone https://github.com/NuSTAR/nuskybgd.git

2. Go to the nuskybgd directory and initialize the environment variaibles

    $ source initialize_nuskybgd.sh

3. Add the the nuskybgd directory to your default IDL path.

	If you have an IDL_STARTUP file already, add the following lines to it:
	
	nuskybgd_code = getenv('NUSKYBGD')                                                                                                                                                                
	!path = expand_path('+'+nuskybgd_code)+':'+ $
                    !path

## UPDATE 8/28/2015

With this 'release', the use of the 'nuabs' XSpec model has been phased out
of nuskybgd routines, so it no longer needs to be installed as a local model.
Routines now assume that the absorption has been included directly in the
response matrices, WHICH IS NOT DONE in nuproducts/numkrmf routines.  The
new routine 'getspecrmf.pro' can extract spectral and RMF files appropriate 
for use with 'nuskybgd_fitab.pro', or the default routines (which the
'getspecnoarf.py' script calls) can be used to produce initial RMFs.  The
'addabs2rmf.pro' routine can then add the detector absorption to the RMFs --
just be sure to update the RESPFILE keyword in the spectra if the resulting
RMF filenames change.
