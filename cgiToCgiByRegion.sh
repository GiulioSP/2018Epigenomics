# Parse CpG BED file into dmrByRegion file
if [ "$1" == "-h" ] ; then
    echo -e "Creates .cgiByRegion file, which parses genome into same sized chunks (\"size\"), and reports regions where most bps are part of a CGI. Each region is reported as: \"chr:positionStart:positionEnd\" and \"1\". Requires CGI position BED file, in this case \"hg38_CGIs\".
Usage: `basename $0` -i <dirIn> -o <dirOut> -f <file> -s <size>
    <dirIn>: input directory
    <dirOut>: output directory
    <file>: input file
    <size>: the size of chunks into which the genome will be parsed"
    exit 0
fi

while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f) file="$2"; shift;;
        -s) size="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut

#generating a file of CGI by region. 
#"n" is used to indicate the index of the genome region being addressed, where for chromosome 1 and region size=400, aaddress(n=0) is "chr1:1:400" and address(n=4) is "chr1:1601:2000"
awk -v size=$size 'BEGIN{print "chr:posStart:posEnd\tis_CGI"} 
{nStart=int($2/size)+int($2*2/size)%2 ; nStop=int($3/size)+int($3*2/size)%2 ; for(n=nStart; n<nStop; n++) print $1":"(n*size+1)":"((n+1)*size)"\t1" ; done}' "$dirIn/$file" > "$dirOut/$file.cgiByRegion.txt"
