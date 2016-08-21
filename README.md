# Quickclump

Identify clumps within a 3D FITS datacube.

## Obtaining Quickclump

GitHub: <https://github.com/vojtech-sidorin/quickclump>  
Direct download link (GitHub):
[quickclump-master.zip](https://github.com/vojtech-sidorin/quickclump/archive/master.zip)  
Web: <http://galaxy.asu.cas.cz/~vosidorin/quickclump.html>

## Quick start guide

To find clumps in your FITS datacube, type

    $ python qc.py my_datacube.fits

To show the usage help, type

    $ python qc.py -h

To run Quickclump in interactive mode, type

     $ python
     >>> import qc
     >>> qc.main(["my_datacube.fits"])

## Description

Quickclump is an automated, fast and accurate clump-finding code.  It was
inspired by another clump-finding code called
[DENDROFIND](http://galaxy.asu.cas.cz/~richard/dendrofind/) ([Wunsch et al.,
2012](http://adsabs.harvard.edu/abs/2012A%26A...539A.116W)).  Quickclump
reaches the same accuracy as DENDROFIND run with the parameter `Nlevels` set to
infinity.  It also runs much faster and its running time and memory usage scale
linearly with the input datacube size.

Quickclump was originally written for analysing astronomical data -- radio
observations of interstellar clouds -- but its use is not limited to astronomy.
Quickclump doesn't presume anything about the origin of the input data.  It
will find clumps in any 3-dimensional rectangular array served in a FITS file.

## Requirements

- Python 2.7/3.4 or above.
- Python packages:
    - numpy
    - astropy or pyfits

## Performance tests

Running Quickclump on the whole LAB HI Survey (721x361x485 pixels)
([Kalberla et al., 2005](http://adsabs.harvard.edu//abs/2005A%26A...440..775K))
on CPU Intel i7-4770S took 10 minutes and the memory utilisation peaked at
870 MiB.  It found 4178 clumps.  Quickclump was run with option
`--otext None`, i.e. producing no text output.

Tests with real CO data from the
[Galactic Ring Survey](http://www.bu.edu/galacticring/)
([Jackson et al., 2006](http://adsabs.harvard.edu//abs/2006ApJS..163..145J))
show that Quickclump consumes up to 10 times the size of the input datacube.
Numpy's `std()` method is especially eager for memory and takes about 6 times
the size of the array (input datacube).  However, if you provide the
parameters `--dTleaf` and `--Tcutoff` at the command-line, the memory-hungry
numpy routines won't be called and the memory usage should stay below 5 times
the size of your input datacube.

## Notes

-   Tested in Python 2.7 and 3.4.
