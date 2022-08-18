# CLEAN codes

Here are the python codes of several types of CLEAN. The inputs (.mat files) are generated based on the Measurement Set [HIP-274](https://jira.skatelescope.org/browse/HIP-274).

To run the codes, simply use:

    python HogCLEAN.py
  
or

    python MsCLEAN.py

## Hogbom CLEAN

INPUT:

    psf.mat   : PSF (dirty beam)
    dirty.mat : Dirty image
    cbeam.mat : CLEAN beam
    
OUTPUT:

    maximum  : Intensity and position of the peaks
    residual : Residual image
    modelnew : Model image
    skymodel : Restored image (the final output we want!)

## Multi-Scale CLEAN

INPUT:

    res_scalestack.mat      : Scaled residual images
    psf_scalescalestack.mat : Scaled and crossed PSFs
    cbeam.mat               : CLEAN beam
    
OUTPUT:

    maximum  : Intensity and position of the peaks
    residual : Residual images
    modelnew : Model image
    skymodel : Restored image (the final output we want!)
