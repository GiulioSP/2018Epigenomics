# Contrast VCF files to isolate variant calls unique to a primary sample, meaning it is not present in the other samples provided. Also outputs amount of entries in original file and new file. 
if [ "$1" == "-h" ] ; then
    echo -e "Contrasts VCF files to isolate variants unique to the main sample (the first in fileList). Can apply filters for higher specificity, at the cost of sensitivity, when -specif is added.
Usage: `basename $0` -i <dirIn> -o <dirOut> -fl <fileList> -fn <fileNickname> -s <size> -specif
    <dirIn>: input directory
    <dirOut>: output directory
    <fileList>: txt file with list of input files. One entry per line
    <fileNickname>: txt file with list of input file nicknames in the same order as fileList. One entry per line. If not specified, full names in fileList will be used.
    <size>: the size of chunks into which the genome will be parsed
    <specif>: "
    exit 0
fi

while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
	-fl) fileList="$2"; fileNickname="$2"; shift;;
	-fn) fileNickname="$2"; shift;;
        -s) size="$2"; shift;;
	-f) specif="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut
declare -a fileAddress
readarray -t fileList < "$dirIn/$fileList"
readarray -t fileNickname < "$dirIn/$fileNickname"

for fileInd in "${!fileList[@]}" 
do
bgzip "$dirIn/${fileList[$fileInd]}"
tabix "$dirIn/${fileList[$fileInd]}.gz" 
fileAddress[$fileInd]="$dirIn/${fileList[$fileInd]}.gz"
done

fileOutName="${fileNickname[*]}" 

bcftools isec -w1 -C ${fileAddress[*]} > "$dirOut/${fileOutName//" "/-}.vcfcontrast.vcf"
echo "amount of VCF entries in contrast for file ${fileOutName//" "/-}.vcfcontrast.vcf:"
grep -c -v ^"#" "$dirOut/${fileOutName//" "/-}.vcfcontrast.vcf"
