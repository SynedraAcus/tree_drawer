### Tree drawer
Tree visualisation code for a paper on intragenic duplications in diatoms.

#### Design

This script takes a Newick tree of domains and tries to highlight different
parts of the same multidomain sequence in the same color, thus facilitating the
manual search for matching clades. The sequences that cannot be matched, or
domains that form a pair of sister leaves (*ie* are not part of a larger clade)
are highlighted in dark gray. The sequences that are the only domain within
their respective sequences are not highlighted. Example tree is shown below:

![example.png](example tree) 
#### Usage
Clone or download this repo, then call `./draw.py`. See `./draw.py -h` for
command line interface.

#### Dependencies
Python3.6+, ete3

#### Licensing
This code is available for free use under the terms of Creative Commons
[CC-By 3.0](https://creativecommons.org/licenses/by/3.0/legalcode) license. This
means that you're welcome to do whatever you want with the contents of this repo
, assuming you do mention that it was originally created by Alexey Morozov aka
synedraacus.