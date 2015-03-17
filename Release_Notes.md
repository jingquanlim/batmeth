#Release\_Notes

### Note ###
BatMeth scans up to 5 mismatches per read. If the rate of mismatch in the read is high, we recommend you to try the -F option which trims a read from the 3' end prior to mapping.

### Next Version ###
Fixed Total Read count output to display.

### Version 1.04b ###
Added a per-base methylation information string to each hit (BatMeth-Illumina version). Please read wiki page for BatMeth-Illumina for more details. http://code.google.com/p/batmeth/wiki/README_ILLUMINA

Fixed some issues that prevent smooth compiling in various Linux platforms. (Maths functions: pow, sqrt, etc were causing some problems under some versions of GCC)

Fixed index-loading glitch.

Fixed index-building glitch.

Added the ability to handle genomes with an arbitrary number of chromosomes.

Reduced the size of index for BatMeth\_SOLiD (temporary files are removed)

Changed BatMeth\_SOLiD default mapping read-length to 36bp. Can be changed by `-F <integer>` option.