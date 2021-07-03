# Parse FindER file into chipByRegion file
if [ "$1" == "-h" ] ; then
    echo -e "Creates vcfByRegion file, which parses genome into same sized chunks (\"size\"), and reports the chunks where one or more variant calls are present, in the format: \"chr:positionStart:positionEnd\",  \"adjustedResidual/max(adjustedResidual)\", and \"chiSquare\". When Requires enrichment BED file: FindER.bed.gz.3col.sorted
Usage: `basename $0` -i <dirIn> -o <dirOut> -f <file> -s <size> -n <name>
    <dirIn>: input directory
    <dirOut>: output directory
    <file>: input file
    <size>: the size of chunks into which the genome will be parsed
    <name>: new name for file. Current file name will be used if this is empty. The postfix is .vcfByRegion"
    exit 0
fi

while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f) file="$2"; name="$2"; shift;;
        -s) size="$2"; shift;;
        -n) name="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut

#generating a file of variant calls by region. the presence of a variant call (VC) is indicated by the "info" column text, separated by a "%" symbol when more than one VC is present. The format is "chr:posStart:posEnd POS|REF|ALT|GT|EFF", where EFF=effect(impact|functional class|codon change or distance|amino acid change|amino acid length|gene name|biotype|gene coding|Transcript ID|exon intron rank|genotype number|[warnings]).
#Region size is determined at "-v size". Currently 400bp on all analyses since it gives ~40% genome coverage with the DMR data chi square constraints, while still being granular compared to the 2kbp region of influence of each CpG (according to Tony's data in 2018). 
#"n" is used to indicate the index of the genome region being addressed, where for chromosome 1 and region size=400, aaddress(n=0) is "chr1:1:400" and address(n=4) is "chr1:1601:2000"
bcftools query -H -f '%CHROM\t%POS\t%REF|%ALT[|%GT|%EFF]\n' "$dirIn/$file" | awk -v size=$size 'BEGIN{print "chr:posStart:posEnd\tVCinfo_"name"_pos|ref|alt|gt|eff"}
(NR==2){n=int($2/size); address=$1":"(n*size+1)":"((n+1)*size); info=$2"|"$3; prevAddress=address}
(NR>2){n=int($2/size); address=$1":"(n*size+1)":"((n+1)*size); if(address==prevAddress){info=info"%"$2"|"$3; prevAddress=address} else{print prevAddress"\t"info; info=$2"|"$3; prevAddress=address}}
END{print prevAddress"\t"info; info=$2"|"$3; prevAddress=address}' > "$dirOut/$name.vcfByRegion.txt"
