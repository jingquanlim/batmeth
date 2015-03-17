

# Installation #


To install the package, do the following.
  1. untar the package:
```
tar -zxvf batmethxxxx_SOLiD
```
  1. change to installation directory batmetxxxx\_SOLiD
```
cd batmethxxxx_SOLiD
```
  1. Configure and Install
```
./configure
make
make copy
```

# Building Indexes. #

All executables of this package can be found in batmethxxxx\_SOLiD/bin folder after `make copy`.

Indexes need to be built by calling a script as follows.
  1. build\_meth GENOME

# Mapping #


To map the data set METH to the genome with base name GENOME do,
```
./batmeth -g GENOME -i INPUT_METH -n <#mismatch on fully-converted genome> -N <#mismatch on non-CpG-converted genome> -o <Output_Prefix>
./split <FINAL_RESULT> GENOME <number of mismatches> n <all...output...files...from...batmeth...>
```

# Output #

`<Final_result>` will contain all uniquely mapped hits in the following format:
> Line 1: `Header`

> Line 2: `<read_ID> <Color Read>`

> Line 3: `<read-origin> <chrNum> <strand> <chrPos> <mismatch> <read_length> <0> <Color read in base space> <genomic reference region>`

> `<read-origin>` is a number from [1..4]. {1,3} means the read is mapped to the + strand and {2,4} means otherwise. Hits labelled with {1,4} are from fully-converted genomes while those labelled with {2,3} are from non-CpG converted genomes.

> `<strand>` is always "+" to indicate coordinates are left-most position of an alignment.

Please ignore the other fields which are appended to the end of the lines.

# Sample Run #
`[linux]$ .batmeth -g genome.fa -i query.fa -n 4 -N 0 -p 4 -o tmp`<br><code>[linux]$ .split output.out genome.fa 3 n tmp.0 tmp.1 tmp.2 tmp.3</code>

<h1>Note on Input format</h1>
SOLID reads usually have their reads and qualities stored in two files. BatMeth requires input to be in csfasta or csfastq (similar to fastq).