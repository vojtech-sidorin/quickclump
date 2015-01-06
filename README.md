# clumpit

Identify clumps within a 3D FITS datacube.

## Obtaining clumpit

GitHub: <https://github.com/vojtech-sidorin/clumpit>  
Direct download link (GitHub):
[clumpit-master.zip](https://github.com/vojtech-sidorin/clumpit/archive/master.zip)  
Web: <http://galaxy.asu.cas.cz/~vosidorin/clumpit.html>

## Quick start guide

To find clumps in your data cube, type

    $ python clumpit.py my_datacube.fits

To show usage help, type

    $ python clumpit.py -h

To run in interactive mode, type

     $ python
     >>> import clumpit
     >>> clumpit.main(["my_datacube.fits"])

## Description

clumpit is an improved implementation of
[DENDROFIND](http://galaxy.asu.cas.cz/~richard/dendrofind/)
([Wunsch et al., 2012](http://adsabs.harvard.edu/abs/2012A%26A...539A.116W))
-- a clump-finding algorithm inspired by
[Clumpfind](http://www.ifa.hawaii.edu/users/jpw/clumpfind.shtml)
([Williams et al., 1994](http://adsabs.harvard.edu/abs/1994ApJ...428..693W)).
DENDROFIND was originally conceived by Richard Wunsch, who also published its
first implementation in Python, later rewritten in C.  Compared to DENDROFIND,
clumpit uses different data structures and doesn't need parameter `Nlevels`.
clumpit is also faster (about 50 000 times) and scales linearly with the data
cube volume (number of pixels).

clumpit was originally written for analysing astronomical data -- radio
observations of molecular clouds -- but its use is not limited to astronomy.
clumpit doesn't presume anything about the origin of the input data.  It will
find clumps in any 3-dimensional rectangular array served in a FITS file.

## Notes

-   Following my tests with real CO data, this program consumes up to 10 times
    the size of the input data cube.  Numpy's `std()` method is especially
    eager for memory and takes about 6 times the size of the array (input data
    cube).  However, if you provide the parameters `--dTleaf` and `--Tcutoff`
    at the command-line, the memory-hungry numpy routines won't be called and
    the memory usage should stay below 5 times the size of your input data
    cube.

-   Tested in Python 2.7 and 3.4.
