#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
    echo "Usage    :"
    echo "    ./metaSLR.sh [OPTION]" 
    echo ""
    echo "phase meta stLFR reads based on reference of each species."
    echo ""
    echo "Options  :"
    echo "        --haplotype   haplotype reference file in fasta format."
    echo "                      ( note : gzip format IS NOT supported. )"
    echo "                      ( note : this script use haplotype file name as species name.)"
    echo "                      (        please make sure the basename of each haplotype is not exactly same !!! ) "
    echo "        --meta        meta stLFR reads file in fastq format."
    echo "                      file in gzip format is accepted, but filename must end by \".gz\"."
    echo "        --thread      threads num."
    echo "                      [ optional , default 8 thread. ]"
    echo "        --memory      x (GB) of memory to initial hash table by jellyfish."
    echo "                      (noted: real memory used maybe greater than this. )"
    echo "                      [ optional , default 20GB. ]"
    echo "        --jellyfish   jellyfish path."
    echo "                      [ optional , default jellyfish. ]"
    echo "        --mer         mer-size"
    echo "                      [ optional , default 21. ]"
    echo "        --lower       ignore mer with count < lower."
    echo "                      [ optional , default 1. ]"
    echo "        --upper       ignore mer with count > upper."
    echo "                      [ optional , default 33. ]"
    echo "        --help        print this usage message."
    echo "        "
    echo "Examples :"
    echo "    ./metaSLR.sh --haplotype  s1.fa --haplotype s2.fa --haplotype s3.fa --meta read.fq.gz"
    echo ""
    echo "    ./metaSLR.sh --haplotype  s1.fa --haplotype s2.fa --haplotype s3.fa \\"
    echo "                     --meta read1.fq --meta read2.fq  "
    echo ""
    echo "    ./metaSLR.sh --haplotype  s1.fa --haplotype s2.fa --haplotype s3.fa \\"
    echo "                     --meta read1.fq --meta read2.fq \\ "
    echo "                     --memory 20 --thread 20 \\"
    echo "                     --mer 21 --lower=1 --upper=33 \\"
    echo "                     --jellyfish /home/software/jellyfish/jellyfish-linux"
}

###############################################################################
# basic variables 
###############################################################################
MER=21
JELLY=jellyfish
CPU=8
MEMORY=10
LOWER=9
UPPER=33
HAPS=""
META=""
SPATH=`dirname $0`
###############################################################################
# parse arguments
###############################################################################
if [[ $# == 0 ]] ; then 
    usage
    exit 0
fi
echo "CMD :$0 $*"
while [[ $# > 0 ]] 
do
    case $1 in
        "-h")
            usage
            exit 0
            ;;
        "--help")
            usage
            exit 0
            ;;
        "--jellyfish")
            JELLY=$2
            shift
            ;;
        "--memory")
            MEMORY=$2
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--lower")
            LOWER=$2
            shift
            ;;
        "--upper")
            UPPER=$2
            shift
            ;;
        "--mer")
            MER=$2
            shift
            ;;
        "--haplotype")
            HAPS=$HAPS" "$2
            shift
            ;;
        "--meta")
            META=$META" "$2
            shift 
            ;;
        *)
            echo "invalid params : \"$1\" . exit ... "
            exit
        ;;
    esac
    shift
done
# print arguments
echo "metaSLR starting with : "
echo "    haplotype input : $HAPS"
echo "    meta input      : $META"
echo "    jellyfish       : $JELLY"
echo "    memory          : $MEMORY GB"
echo "    thread          : $CPU "
echo "    mer             : $MER "
echo "    lower           : $LOWER"
echo "    upper           : $UPPER"
echo "metaSLR.sh in dir   : $SPATH"

CLASSIFY=$SPATH"/classify"
FILTER_FQ_BY_BARCODES_AWK=$SPATH"/filter_fq_by_barcodes.awk"
# sanity check
if [[ $MEMORY -lt 1  || $CPU -lt 1 || \
    -z $HAPS || -z $META || \
    -z $JELLY  || $MER -lt 11 || \
    $LOWER -lt 1 || $UPPER -gt 100000000 ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi
if [[ ! -e $CLASSIFY ]] ; then 
    echo "ERROR : please run \"make\" command in $SPATH before using this script! exit..."
    exit 1
fi
if [[ ! -e $FILTER_FQ_BY_BARCODES_AWK ]] ; then
    echo "ERROR : \"$FILTER_FQ_BY_BARCODES_AWK\"  is missing. please download it from github. exit..."
    exit 1
fi
for x in $HAPS $META
do
   if [[ ! -e $x ]] ; then 
       echo "ERROR : input file \"$x\" is not exist ! exit ..."
       exit 1
   fi
done
date
echo "__START__"
###############################################################################
# extract paternal.unique.filter.mer & maternal.unique.filter.mer
###############################################################################
echo "extract unique mers by jellyfish ..."
echo " stage 0 extract filter mers by jellyfish ..."
for x in $HAPS
do
    species=`basename $x`
    # mer-count
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C    -o $species".t0.js"  $x
    # dump filter mers
    $JELLY dump -L $LOWER -U $UPPER $species".t0.js" -o $species".t0.mer.filter.fa"
    # dump all mers
    $JELLY dump $species".t0.js"                     -o $species".t0.mer.fa"
done
date
echo " stage 1 extract unique mers by jellyfish ..."
for x in $HAPS
do
    species=`basename $x`
    others=""
    for ot in $HAPS
    do
        tn=`basename $ot`
        if [[ $tn != $species ]] ; then
            name=$tn".t0.mer.fa"
            others=$others" "$name" "$name
        fi
    done
    #echo "    Mix : "$species".t0.mer.fa  "$others
    # mix 1 copy of species mers and 2 copy of all other species mers
    cat $species".t0.mer.fa"  $others >$species".t1.mixed.fa"
    # count mixed mers
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o $species".t1.js" $species".t1.mixed.fa"
    # count==1 refer to species unique mers
    $JELLY dump -L 1 -U 1 $species".t1.js"          >$species".t1.mer.unique.fa"
done
rm -rf *.t1.mixed.fa
date
echo " stage 2 extract unique & filter mers by jellyfish ..."
for x in $HAPS
do
    species=`basename $x`
    # mix unique mers and filter mers
    cat  $species".t1.mer.unique.fa" $species".t0.mer.filter.fa" >$species".t2.mixed.fa"
    #echo "Mix "$species".t1.mer.unique.fa "$species".t0.mer.filter.fa"
    # count unique and filer mers
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o $species".t2.js" $species".t2.mixed.fa"
    # extrat both unique and filter mers
    $JELLY dump -t -c -L 2 -U 2    $species".t2.js" | awk '{print $1}' >$species".t2.unique.filter.mer"
done
rm -rf *.t2.mixed.fa
echo "extract unique mers done..."
date
###############################################################################
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
echo "extract unique barcode by classify ..."
for x in $HAPS
do 
    species=`basename $x`
    HAPINPUT=$HAPINPUT" --hap "$species".t2.unique.filter.mer"
done

for x in $META
do 
    READ="$READ"" --read ""$x"
done

$CLASSIFY $HAPINPUT $READ  --thread $CPU >phased.barcodes 2>phased.log
date
index=0
echo "parase phased.barcodes now ..."
for x in $HAPS
do
    species=`basename $x`
    awk '{if($2 == i) print $1;}' i=$index  phased.barcodes >$species".unique.barcodes"
    ((index++))
done
awk '{if($2 == "-1") print $1;}' phased.barcodes >homozygous.unique.barcodes
echo "extract unique barcode done"

###############################################################################
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
date
echo "phase reads ..."
for x in $META
do
    name=`basename $x`
    if [[ ${name: -3} == ".gz" ]] ; then 
        for file in $HAPS
        do
            species=`basename $file`
            gzip -dc $x |  awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK $species".unique.barcodes" - > $species"."$name
        done
        gzip -dc $x | awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK homozygous.unique.barcodes - >"homozygous."$name
    else
        for file in $HAPS
        do
            species=`basename $file`
            awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK $species".unique.barcodes" $x > $species"."$name
        done
        awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK homozygous.unique.barcodes $x >"homozygous."$name
    fi
done
echo "phase reads done"
date
echo "__END__"
