# metaSLR
phase stLFR reads based on reference ( or assembly scaffold ) of each species.

## INSTALL

```
git clone https://github.com/BGI-Qingdao/metaSLR.git
cd metaSLR
make
```

## USAGE

```
Usage    :
    ./metaSLR.sh [OPTION]

phase meta stLFR reads based on reference of each species.

Options  :
        --haplotype   haplotype reference file in fasta format.
                      ( note : this script use haplotype file name as species name.)
                      (        please make sure the basename of each haplotype is not exactly same !!! )
        --meta        meta stLFR reads file in fastq format.
                      file in gzip format is accepted, but filename must end by ".gz".
        --thread      threads num.
                      [ optional , default 8 thread. ]
        --memory      x (GB) of memory to initial hash table by jellyfish.
                      (noted: real memory used maybe greater than this. )
                      [ optional , default 20GB. ]
        --jellyfish   jellyfish path.
                      [ optional , default jellyfish. ]
        --mer         mer-size
                      [ optional , default 21. ]
        --lower       ignore mer with count < lower.
                      [ optional , default 1. ]
        --upper       ignore mer with count > upper.
                      [ optional , default 33. ]
        --help        print this usage message.

Examples :
    ./metaSLR.sh --haplotype  s1.fa --haplotype s2.fa --haplotype s3.fa --meta read.fq.gz

    ./metaSLR.sh --haplotype  s1.fa --haplotype s2.fa --haplotype s3.fa \
                     --meta read1.fq --meta read2.fq

    ./metaSLR.sh --haplotype  s1.fa --haplotype s2.fa --haplotype s3.fa \
                     --meta read1.fq --meta read2.fq \
                     --memory 20 --thread 20 \
                     --mer 21 --lower=1 --upper=33 \
                     --jellyfish /home/software/jellyfish/jellyfish-linux
```

Enjoy !
