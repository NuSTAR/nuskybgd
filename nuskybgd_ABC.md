# nuskybgd ABD Guide

---

# Table of Contents
1. [Disclaimer and Scope](#disclaimer)
2. [About](#about)
3. [Setup](#setup)
4. [Walkthrough](#walkthrough)

---

## Disclaimer

Questions, comments, or bugs should (please!) be reported via the [GitHub issues page](https://github.com/NuSTAR/nuskybgd).

## Scope

This is an ABC-like example for using the `nuskbygd` tasks.

It's assumed that you have a working installation of the following:

1. [HEASoft](https://heasarc.gsfc.nasa.gov/lheasoft/) including `XSPEC`.
2. The [NuSTAR CALDB](http://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_intro.html). 
3. IDL (for now).
4. ds9 (for defining source and background regions).
4. `nuskybgd` 
	* See the [README](README.md) file for how to set up the `nuskybgd` code in your local IDL installation.

> **NOTE**: We recommend maintaining a record of your work in a separate IDL script file that can be copy and pasted to the IDL command line. `nuskbygd` is interactive enough that we do *not* generally recommend automating the analysis. 

---

## About

`nuskybgd` tasks are primarily meant to determine the background from blank sky
regions within an observation and use that to produce background spectra for
source regions and to produce background images.

It is also possible to  create 'default' background spectra and images if there are no pure background regions in an observation.  However, one must take care because the various background components vary; for example, the unfocused, stray light 'aperture bgd' component that dominates below ~20 keV can vary +/-10% due to cosmic variance.  That said, for a quick background subtraction, especially of images for presentation purposes, the default parameters do a good job.

 The CXB components assume a fiducial normalization from HEAO-1 by default.
 For the Boldt '87 CXB normalization from HEAO-1 in the 10-15 keV band,
 the expected count rate is 3.2e-4 cts/s/cm^2/deg^2 (pulled from simulated
 aperture spectra created below).  The image created by projinitbgds is in
 the 10-15 keV band assuming a 100 ks exposure; it's just a counts image,
 but the function that produces it has a normalization in units of
 cts/deg^2/mm^2, so the normalization factor should be 0.32 cts/deg^2/mm^2.
 This is the value given in the data file nomapbgdparams.dat.

---

## Setup

Two environmental variables need to be correctly set:

* `$CALDB` points to your *NuSTAR* CALDB, which should already be set. This should already be the case if you are analyzing *NuSTAR* data.
   
* `$NUSKYBGD` and `$NUSKYBGD_AUXIL` point to this repository and the 'auxil' directory in this repository. If you have initialized your installation using the `initialize_nuskybgd.sh` script then this should already be set.
   
You can always check to see if these variables are set by doing the following:

```bash
echo $NUSKBYGD
/Users/yourname/git/nuskybgd
```





---

# Walkthrough

## Set the path that you're going to use.

Here you tell `nuskybgd` where the data resides. For simplicity, we currently assume that all of the data and spectral products live underneath the SEQID directory that you downloaded from the HEASARC.

`dir` is the directory with respect to wherever you invoked IDL to run `nuskybgd`. Normally this is just above your data directory, in which case you get the following:

```IDL
dir='./'
obsid='50002031004'
```

## Make an instrument map

This grabs the instrument map from the caldb and makes a 'pixel map' that
assigns each RAW pixel to a unique DET1 pixel.  This allows bad pixels to 
be excluded in the same way as is done in the pipeline.  It automatically 
excludes all identified bad pixels -- including any supplied by the
user -- during pipeline processing.  You can exclude additional pixels,
but in practice should not -- instead, rerun the pipeline with the
updated bad pixel list and run this again.

```IDL
nuskybgd_instrmap,dir,obsid,'A','bgd'
nuskybgd_instrmap,dir,obsid,'B','bgd'
```

The last input, `bgd`, specifies the directory all the ancilliary observation-specific files are placed.  This directory is placed under the observation's event_cl/ and will need to be fed to later routines so it puts all the intermediate files in the same place and knows where to grab them from.

## Create an image and background regions with DS9

The image is necessary because its header is pasted on the reference
and background images so it has WCS coordinates. We have provided a simple `python` wrapper script (`mkimgs.py`) to do this that calls `nuproducts` to produce the image:

```bash
./mkimgs.py ./ 50002031004 3 20
```

The first argument is the path (in this example we assume that you have copied the `mkimgs.py` script to the same location as where you're running `nuskybgd`. The second argument is the SEQID, the third and fourth are the low and high energy bounds for the image. These can generally be set to 3 and 20 keV.

Although not strictly necessary, it helps to have multiple bgd regions
that span the gradient of the aperture component, so it can be isolated
from the other two low energy components.  Annuli are good choices, with
perhaps two or three concentric annuli around the source.



## Define background regions and extract background spectra

**Note that the following information needs to be updated on how to use the new scripts to properly scale the RMF for the photon path absorption!!! What follows is slightly out of date, but will work except if you really, really care about what happens below ~5 keV.**

It's up to the user to define regions that are source-free to determine the background. In the Word document in the [`docs`](./docs) directory there is an example of background region for a source. For the remainder we assume that you have specified two background regions for each instrument.

We have provided a `python` wrapper script to `nuproducts` that generates background spectra without producing any ancillary response functions (ARFs) since ARFs are not necessary for the `nuskybgd` analysis.

```bash
./getspecnoarf.py ./ 50002031004 bgd1A bgdspec A
./getspecnoarf.py ./ 50002031004 bgd1B bgdspec B
./getspecnoarf.py ./ 50002031004 bgd2A bgdspec A
./getspecnoarf.py ./ 50002031004 bgd2B bgdspec B
```


> **NOTE**: The chi^2 statistic is used, so be sure to group the bgd spectra.
In the above script, the output spectra are bgd1A_sr_g30.pha, etc.,
with bins grouped to at least 30 counts.  Also, background regions
do not need to avoid gaps, dead pixels, or areas where there is no
data -- these areas are weighted appropriately.  You are limited to
ds9 region types 'circle', 'ellipse', 'box', and 'polygon', and
multiple include and exclude regions can be used (they are 'ORed' 
together).

## Fit the background spectra

In this step, all the background spectra are fit simultaneously to an
empirical model in `XSPEC`, assuming the spatial shape of the aperture bgd
and the relative instrumental line strengths between detectors.

> **NOTE**: The current assumptions are based on ~0.5 Msec of data, some of which are from early in the mission, and will likely be revised.

First, adjust the following lines based on how many background regions you have:

```IDL
core='bgd'+['1','2','0']
bgdreg=[core+'A'+'.reg', core+'B'+'.reg']
spec=[core+'A'+'_sr.pha', core+'B'+'_sr.pha']
```

Now we call the main `nuskbygd_fitab` script below. 



```IDL
nuskybgd_fitab,dir,obsid,bgdreg,'bgdspec',spec,'AB','bgd'
```


The last input, 'bgd', specifies the directory all the ancilliary 
observation-specific files are placed.  The files take some time to create,
but only need to be created once as long as the PA used doesn't change.

The default behavior for `nuskbygd_fitab` is to run a fit in `XSPEC` and thenstop so that the user can evaluate the fit to the background. 

If there is non-bgd emission in your regions, such as from the wings of a bright source or extended source emission, you are free to add a model that describes the spectrum of the additional component.

When finished, run the script it tells you to (before you started fiddling), which will save the parameters and exit. This looks like:

```

XSPEC12> @50002031004/event_cl/bgdspec/bgdparams.xcm
```

> **NOTE**: The output parameter file is used to create spectra and images and by default is called bgdfitparams[A/B].dat and is placed in the same directory as your bgd spectra.  It is assumed to be there unless the file is directly pointed to using the paramfile keyword.
    
## Produce simulated background spectra for analysis 

> **NOTE**: This assumes that you have already produced spectral files (PHA, RMF, and ARF) for your source using `nuproducts`.

In principle we now know what the background is doing everywhere, so we can generate (using `fakeit` in `XSPEC`) a spectrum for any arbitrary region.

The scripts below assume that there is a source region called 'src.reg' in your `event_cl` directory and the actual spectra in `event_cl/src/`. Alternatively you can specify the name of the spectrum and the location of the spectra by using various inputs to `nuskbygd_spec`.

Example:

```IDL
nuskybgd_spec,dir,obsid,'srcA.reg','src','bgdspec','A','bgd', specname='nu50002031004A_sr.pha'
nuskybgd_spec,dir,obsid,'srcA.reg','src','bgdspec','B','bgd', specname='nu50002031004B_sr.pha'
```

This will produced files called `bgdnu50002031004A_sr.pha` and `bgdnu50002031004B_sr.pha` in the `src` directory. This automatically updates the relevant keywords in the source PHA files so that `XSPEC` will automatically load the simulated backgrounds during spectral analysis.

> **NOTE**: `nuskbgyd_spec` by default produces a "smooth" background by simulating 100x more background exposure than source exposure. This is controlled by the `expfactor` keyword. If you want to simulate more/less background (i.e. if you want to get some estimate of Poisson effects on the background), then you may want to set `expfactor=1` in the above calls and generate a large number of background spectra. See the help text in `nuskbygd_spec.pro` for other options.





## Produce simulated background images

You can also produce simulated background images using the following:

```IDL
nuskybgd_image,dir,obsid,'imA3to20keV.fits',3,20,'A','bgd','bgdspec'
nuskybgd_image,dir,obsid,'imB3to20keV.fits',3,20,'B','bgd','bgdspec'
```

Here `imA3to20keV.fits` is an image either produced using nuproducts or using the example `python` wrapper script as described above.

See the help text in `nuskbygd_image.pro` for other available options.

> **NOTE**: The background images produced in this way are scaled to the exposure in the input file but the counts are *not* quantized, so the units on the resulting background image are in fractional counts per pixel.
