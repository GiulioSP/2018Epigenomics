# Parse CpG BED file into dmrByRegion file
if [ "$1" == "-h" ] ; then
    echo -e "Creates .dmrByRegion file, differentialy methylated region, which parses genome into same sized chunks (\"size\"), and reports regions where the methylated CpG count of the files inputed is sufficiently different to pass a chi squared test. Each differentially methylated region is reported as: \"chr:positionStart:positionEnd\", for each file \"adjustedResidual/max(adjustedResidual)\", and \"chiSquare\". Where adjustedResidual=(obs–exp)*total/√[exp*(total-rowTotal)*(total-columnTotal)]. Requires methylation counts BED file: \".5mC.CpG.gz\". Will create \".methyByRegion.txt\" files for each input file, as well as a \".methyByRegion.join.txt\". This script can compare at most 28 files at once; to increase this number, edit the code to add more chi squared values to the \"chi\" vector. 
Usage: `basename $0` -i <dirIn> -o <dirOut> -f1 <file1> -fn1 <fileNickname1> -f2 <file2> -fn2 <fileNickname2> -f3 <file3> -fn3 <fileNickname3> -f4 <file4> -fn4 <fileNickname4> -f5 <file5> -fn5 <fileNickname5> -f6 <file6> -fn6 <fileNickname6> -jfn <joinFileName> -s <size>
    <dirIn>: input directory
    <dirOut>: output directory
    <f1>: name of file 1
    <fn1>: nickname of file 1, used in column names, prefferably short. If absent, full file name will be used
    <f2>: name of file 2
    <fn2>: nickname of file 2, used in column names, prefferably short. If absent, full file name will be used
    <f3>: name of file 3
    <fn3>: nickname of file 3, used in column names, prefferably short. If absent, full file name will be used
    <f4>: name of file 4
    <fn4>: nickname of file 4, used in column names, prefferably short. If absent, full file name will be used
    <f5>: name of file 5
    <fn5>: nickname of file 5, used in column names, prefferably short. If absent, full file name will be used
    <f6>: name of file 6
    <fn6>: nickname of file 6, used in column names, prefferably short. If absent, full file name will be used
    <joinFileName>: name of the \".methyByRegion.join\" file being created for this grouping 
    <size>: the size of chunks into which the genome will be parsed"
    exit 0
fi

while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f1) file1="$2"; fileNickname1="$2"; shift;;
        -fn1) fileNickname1="$2"; shift;;
        -f2) file2="$2"; fileNickname2="$2"; shift;;
        -fn2) fileNickname2="$2"; shift;;
        -f3) file3="$2"; fileNickname3="$2"; shift;;
        -fn3) fileNickname3="$2"; shift;;
        -f4) file4="$2"; fileNickname4="$2"; shift;;
        -fn4) fileNickname4="$2"; shift;;
        -f5) file5="$2"; fileNickname5="$2"; shift;;
        -fn5) fileNickname5="$2"; shift;;
        -f6) file6="$2"; fileNickname6="$2"; shift;;
        -fn6) fileNickname6="$2"; shift;;
        -jfn) joinFileName="$2"; shift;;
        -s) size="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut
fileList=("$file1" "$file2" "$file3" "$file4" "$file5" "$file6")
fileNickname=("$fileNickname1" "$fileNickname2" "$fileNickname3" "$fileNickname4" "$fileNickname5" "$fileNickname6")

initFlag=1
for fileInd in "${!fileList[@]}" 
do
#generating a join file of methylation by region
awk -v size=$size -v fileNickname=${fileNickname[$fileInd]} 'BEGIN{print "chr:posStart:posEnd\tM_"fileNickname"\tUM_"fileNickname; chrom=$1; countM=0; countUM=0; pos=0}
{if(chrom!=$1){if((countM!=0)||(countUM!=0)){print chrom":"(pos+1)":"(pos+size)"\t"countM"\t"countUM; countM=0; countUM=0} chrom=$1; pos=0}
while($2>=pos+size){if((countM!=0)||(countUM!=0)){print chrom":"(pos+1)":"(pos+size)"\t"countM"\t"countUM; countM=0; countUM=0} pos+=size}
if(($2>=pos)&&($2<pos+size)&&(chrom==$1)){countM+=$4; countUM+=$5}}
END{if((countM!=0)||(countUM!=0)){print chrom":"(pos+1)":"(pos+size)"\t"countM"\t"countUM}}' <(gzip -dc "$dirIn/${fileList[$fileInd]}") > "$dirOut/${fileNickname[$fileInd]}.methyByRegion"
echo "Generated .methyByRegion file for sample ${fileNickname[$fileInd]}"
#joining methylationByRegion files into a single one for analysis
#initiation of "join" with the first file called
if (($initFlag))
    then
        cat "$dirOut/${fileNickname[$fileInd]}.methyByRegion" > "$dirOut/$joinFileName.methyByRegion.join"
        initFlag=0
        echo "Sample ${fileNickname[$fileInd]} incorporated in $joinFileName.methyByRegion.join"
    else
        join -j1 "$dirOut/$joinFileName.methyByRegion.join" "$dirOut/${fileNickname[$fileInd]}.methyByRegion" > "$dirOut/methyByRegion.jointemporary"
        mv "$dirOut/methyByRegion.jointemporary" "$dirOut/$joinFileName.methyByRegion.join"
        echo "Sample ${fileNickname[$fileInd]} incorporated in $joinFileName.methyByRegion.join"
    fi
done

#collecting chi square test statistic from df value and 95% confidence, such that ${chi[$df]}=chiTest:
chi=(0 3.841 5.991 7.815 9.488 11.07 12.592 14.067 15.507 16.919 18.307 19.675 21.026 22.362 23.685 24.996 26.296 27.587 28.869 30.144 31.41 32.671 33.924 35.172 36.415 37.652 38.885 40.113)
df=5

#Chi squared analysis to identify differentially methylated regions
awk -v chiTest=${chi[$df]} -v df=$df 'NR==1{print $1"\tStdRes_"$2"\tStdRes_"$4"\tStdRes_"$6"\tStdRes_"$8"\tStdRes_"$10"\tStdRes_"$12"\tChi_Total>"chiTest}
NR>1{ cellsLessThan5=0
totColMeth= $2+$4+$6+$8+$10+$12
totColUnmeth= $3+$5+$7+$9+$11+$13
totRow1= $2+$3
totRow2= $4+$5
totRow3= $6+$7
totRow4= $8+$9
totRow5= $10+$11
totRow6= $12+$13
total= totColMeth+ totColUnmeth
if (totColMeth==0||totColUnmeth==0||totRow1==0||totRow2==0||totRow3==0||totRow4==0||totRow5==0||totRow6==0) {next}

e_1M= totColMeth*totRow1/total
e_2M= totColMeth*totRow2/total
e_3M= totColMeth*totRow3/total
e_4M= totColMeth*totRow4/total
e_5M= totColMeth*totRow5/total
e_6M= totColMeth*totRow6/total

e_1UM= totColUnmeth*totRow1/total
e_2UM= totColUnmeth*totRow2/total
e_3UM= totColUnmeth*totRow3/total
e_4UM= totColUnmeth*totRow4/total
e_5UM= totColUnmeth*totRow5/total
e_6UM= totColUnmeth*totRow6/total

if (e_1M<1||e_1UM<1||e_2M<1||e_2UM<1||e_3M<1||e_3UM<1||e_4M<1||e_4UM<1||e_5M<1||e_5UM<1||e_6M<1||e_6UM<1){next}

if (e_1M<5){cellsLessThan5+=1}
if (e_1UM<5){cellsLessThan5+=1}
if (e_2M<5){cellsLessThan5+=1}
if (e_2UM<5){cellsLessThan5+=1}
if (e_3M<5){cellsLessThan5+=1}
if (e_3UM<5){cellsLessThan5+=1}
if (e_4M<5){cellsLessThan5+=1}
if (e_4UM<5){cellsLessThan5+=1}
if (e_5M<5){cellsLessThan5+=1}
if (e_5UM<5){cellsLessThan5+=1}
if (e_6M<5){cellsLessThan5+=1}
if (e_6UM<5){cellsLessThan5+=1}
if (cellsLessThan5>((df+1)*2/5)){next}

oe_1M= $2-e_1M
oe_1UM= $3-e_1UM
oe_2M= $4-e_2M
oe_2UM= $5-e_2UM
oe_3M= $6-e_3M
oe_3UM= $7-e_3UM
oe_4M= $8-e_4M
oe_4UM= $9-e_4UM
oe_5M= $10-e_5M
oe_5UM= $11-e_5UM
oe_6M= $12-e_6M
oe_6UM= $13-e_6UM

chi_1M= (oe_1M)^2/e_1M
chi_1UM= (oe_1UM)^2/e_1UM
chi_2M= (oe_2M)^2/e_2M
chi_2UM= (oe_2UM)^2/e_2UM
chi_3M= (oe_3M)^2/e_3M
chi_3UM= (oe_3UM)^2/e_3UM
chi_4M= (oe_4M)^2/e_4M
chi_4UM= (oe_4UM)^2/e_4UM
chi_5M= (oe_5M)^2/e_5M
chi_5UM= (oe_5UM)^2/e_5UM
chi_6M= (oe_6M)^2/e_6M
chi_6UM= (oe_6UM)^2/e_6UM

chi_Total= chi_1M+ chi_1UM+ chi_2M+ chi_2UM+ chi_3M+ chi_3UM+ chi_4M+ chi_4UM+ chi_5M+ chi_5UM+ chi_6M+ chi_6UM

ar_1M= oe_1M/sqrt(e_1M*(total-totRow1)*(total-totColMeth)/total^2)
ar_2M= oe_2M/sqrt(e_2M*(total-totRow2)*(total-totColMeth)/total^2)
ar_3M= oe_3M/sqrt(e_3M*(total-totRow3)*(total-totColMeth)/total^2)
ar_4M= oe_4M/sqrt(e_4M*(total-totRow4)*(total-totColMeth)/total^2)
ar_5M= oe_5M/sqrt(e_5M*(total-totRow5)*(total-totColMeth)/total^2)
ar_6M= oe_6M/sqrt(e_6M*(total-totRow6)*(total-totColMeth)/total^2)

if(chi_Total>=chiTest){print $1"\t"ar_1M"\t"ar_2M"\t"ar_3M"\t"ar_4M"\t"ar_5M"\t"ar_6M"\t"chi_Total}}'  "$dirOut/$joinFileName.methyByRegion.join" > "$dirOut/$joinFileName.chiSqMethyByRegion.txt"

awk '(NR==1){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}
(NR>1){max=sqrt($2^2) ; for (i=3; i <=(NF-1); i++) if(sqrt($i^2)>max) max=sqrt($i^2); print $1"\t"$2/max"\t"$3/max"\t"$4/max"\t"$5/max"\t"$6/max"\t"$7/max"\t"$8}' "$dirOut/$joinFileName.chiSqMethyByRegion.txt" > "$dirOut/$joinFileName.DMR.txt"
