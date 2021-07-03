# Parse FindER file into chipByRegion file
if [ "$1" == "-h" ] ; then
    echo -e "Creates chipByRegion file, which parses genome into same sized chunks (\"size\"), and reports the chunks where most bps are considered enriched by this histone modification, in the format \"chr:positionStart:positionEnd\" and \"1\". Requires enrichment BED file: \"FindER.bed.gz.3col.sorted\" The postfix will be \".chipByRegion\".
Usage: `basename $0` -i <dirIn> -o <dirOut> -f <file> -c <ChIPtype> -s <size> -n <name>
    <dirIn>: input directory
    <dirOut>: output directory
    <file>: input file
    <ChIPtype>: enrichment type to be added as a postfix. examples: \"H3K9me3\", \"H3K27ac\"
    <size>: the size in bps of chunks into which the genome will be parsed. Common size is 400
    <name>: new name for file. Current file name will be used if this is empty. The postfix is .ChIPtype.chipByRegion"
    exit 0
fi

while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f) file="$2"; name="$2"; shift;;
        -c) ChIPtype="$2"; shift;;
        -s) size="$2"; shift;;
        -n) name="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut

#"n" is used to indicate the index of the genome region being addressed. For example, for chromosome 1 and region size=400, address(n=0) is "chr1:1:400" and address(n=4) is "chr1:1601:2000"
awk -v size=$size -v name=$name 'BEGIN{print "chr:posStart:posEnd\tChIP_"ChIPtype"_"name} {nStart=int($2/size)+int($2*2/size)%2 ; nStop=int($3/size)+int($3*2/size)%2 ; for(n=nStart; n<nStop; n++) print $1":"(n*size+1)":"((n+1)*size)"\t1" ; done}' "$dirIn/$file" > "$dirOut/$name.$ChIPtype.chipByRegion.txt"
