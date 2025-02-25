#!/bin/bash
##########


########### ************************ Auburn Easley Specific notes /  May 2nd 2023 **************** ###############

# usage: nohup /home/jgo0012/juicedir/scripts/juicer.sh -g Pan_dem_hyd_era_chrom -z /home/jgo0012/juicedir/references/Pan_dem_hyd_era_chrom.fasta -p /home/jgo0012/juicedir/references/Pan_dem_hyd_era_chrom.fasta.sizes -D /home/jgo0012/juicedir -d /home/jgo0012/juicedir -y /home/jgo0012/juicedir/restriction_sites/pan_erato_DpnII.txt -s DpnII -q gpu2 -l gpu2

# nohup de-attaches the job from login mode so you can logout and the job will continue. 
#-g genome ID prefix
#-z actual genome that samples will be mapped to 
#-p chromozome sizes file 
#-D TOP directory name 
#-d Juicer directory which should have scripts/ references/ and restriction_sites/ underneath it 
#-y restriction site file .txt if rest. enzyme not used write -y none 
#-s name of restriction enzyme (ie. DpnII), if rest. enzyme not used write -y none, if data comes from the Arima Hi-C kit and endonucleases were not used, you may enter Arima here (which contains a combo of restriction enzymes)
#-q queue. Enter partition here (ie gpu2, bac0071_amd. I also manually set this in juicer.sh script to avoid potential conflicts)
#-l long queue. Enter partition here (ie gpu2, bac0071_amd. I also manually set this in juicer.sh script to avoid potential conflicts)
#-S stage: if first run, no stage is required. -S merge, if mapping (.sam) is done but merged_sort.txt file hasnt been created -S dedup after that if merged_nodups.txt hasnt been created, -S final after that to create .hic files

# in setting queues, though it would make sense to write-in different partitions, i've had best results when writing-in the same partition for all (i.e "gpu2")

# many changes done to this script dont make a lot of sense but are part of the reason certain stages work on easley. I.e hashing out all 'module load awk' commands (in Easley, awk is pre-loaded, so it shouldnt matter), or options within this script that should be overridden (but aren't) by our usage statement: ie specifying -q gpu2)
# the debug folder is your best friend, know where error files come from by cmd + f the fixed title of bug file and searching it on this script 

#After mapping completes and all splits are merged into merged_sort.txt in your aligned (results) folder, you may need to run the following script separetly in order to bypass a permissions issue that blocks juicer.sh from creating a 
# job that produces merged_nodups.txt, which is the end file of the standard juicer.sh run. After you run this file and have merged_nodups.txt, run the same juicer.sh command from the start with the same usage except for specifying the stage -final (-S final). 
# This last run on -S final will create the final Hi-C files that you can take to juicebox for visualisation




######### start of additional script for bypassing permissions issue after obtaining merged_sort.txt ######

#module load java

#juicer is being a plonker and made merged_sort.txt file on aligned folder but stopped short of making merged_nodups.txt because of a permission issue where jobs arent allowed to create jobs (is Easley being a plonker?)
#For this we will bypass the permissions and run dups.awk directly on merged_sort.txt to create the file we need.


#cd /home/jgo0012/juicedir


#outputdir=/home/jgo0012/juicedir/aligned_HW2
#juiceDir=/home/jgo0012/juicedir


#touch ${outputdir}/dups.txt
#touch ${outputdir}/optdups.txt
#touch ${outputdir}/merged_nodups.txt
#awk -f ${juiceDir}/scripts/dups.awk -v name=${outputdir}/ ${outputdir}/merged_sort.txt
# for consistency with cluster naming in split_rmdups
#mv ${outputdir}/optdups.txt ${outputdir}/opt_dups.txt


# once the merged_nodups file is complete, relaunch Juicer in final mode by sending in the flag “-S final”





#You may still get an error in producing a contact domains and loops files as "Hi-C map may be too sparse" we will re-run this script with --ignore_sparsity flag, flag --cpu is essential, since the gpu version wont work with our cluster

#/usr/bin/java -Xmx40g -Djava.io.tmpdir=/scratch/Counterman_Lab_Shared/Hi-C/tmp -jar /home/jgo0012/juicedir/scripts/juicer_tools.jar arrowhead --cpu ${outputdir}/inter.hic ${outputdir}/inter_contact_domains/ --ignore_sparsity

#/usr/bin/java -Xmx40g -Djava.io.tmpdir=/scratch/Counterman_Lab_Shared/Hi-C/tmp -jar /home/jgo0012/juicedir/scripts/juicer_tools.jar hiccups --cpu ${outputdir}/inter.hic ${outputdir}/inter_loops/ --ignore_sparsity







##### end of additional script #####


########### ************************  ^^^ Auburn Easley Specific notes ^^^  **************** ###############

#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# Alignment script. Sets the reference genome and genome ID based on the input
# arguments (default human, none). Optional arguments are the queue for the 
# alignment, description for stats file, 
# stage to relaunch at, paths to various files if needed,
# chunk size, path to scripts directory, and the top-level directory (default 
# current directory). In lieu of setting the genome ID, you can instead set the
# reference sequence and the chrom.sizes file path, but the directory 
# containing the reference sequence must also contain the BWA index files.
#
# Splits the fastq files, creates jobs to align them, creates merge jobs that
# wait for the alignment to finish, and creates a final merge job.
#
# Also creates "cleanup" jobs that at each stage, deletes jobs off the cluster
# if any one of them fails.
#
# If all is successful, takes the final merged file, removes name duplicates,
# removes PCR duplicates, and creates the hic job and stats job.  Final
# product will be hic file and stats file in the aligned directory.
#                                                                       
# [topDir]/fastq  - Should contain the fastq files. This code assumes that
#                   there is an "R" in the appropriate files, i.e. *R*.fastq
# From the top-level directory, the following two directories are created:
#                                                                              
# [topDir]/splits  - Where to write the scratch split files (fastq files and
#                    intermediate SAM files). This can be deleted after 
#                    execution.
# [topDir]/aligned - Where to write the final output files.
#
# The following globals should be set correctly before proceeding:
#
# splitsize - The number of lines that each split fastq should contain. Larger
#             means fewer files and longer overall, but too small means there
#             are so many jobs that the cluster won't run them. This can be
#             set with the -C command as well
# read1str  - portion of fastq filename that indicates this is the "read 1"
#             file; used to loop over only the read 1 and within that loop,
#             also align read 2 and merge.  If this is not set correctly,
#             script will not work. The error will often manifest itself
#             through a "*" in the name because the wildcard was not able to
#             match any files with the read1str.   
# Juicer version 1.6
shopt -s extglob
juicer_version="1.6"
## Set the following variables to work with your system

# Aiden Lab specific check
isRice=$(host $(hostname) | gawk '{if ($1~/rice/){print 1}else {print 0}}') #'
isBCM=$(host $(hostname) | gawk '{if ($1~/bcm/){print 1}else {print 0}}') #'
isVoltron=0
## path additionals, make sure paths are correct for your system
## use cluster load commands
if [ $isRice -eq 1 ] 
then
    export PATH=$HOME/bin:$PATH
    isNots=$(host $(hostname) | gawk '{if ($1~/nots/){print 1}else {print 0}}') #'
    if [ $isNots -eq 1 ]
    then
    load_bwa="module load bwa/2.0"
    load_java="module load java/1.8" 
    load_gpu="module load cuda11.0/toolkit" 
    else
    load_bwa="export PATH=/tools/bwa-0.7.17:$PATH"
    load_java="module load java" 
    load_gpu="module load gcccuda/2016a;module load CUDA/8.0.54;module load coreutils;module load cuda11.0/toolkit;" 
    fi
    
    # Juicer directory, contains scripts/, references/, and restriction_sites/
    # can also be set in options via -D
    juiceDir="/home/jgo0012/juicedir" ### RICE
    # default queue, can also be set in options via -q
    queue="gpu2"
    queue_time="24:00:00"
    # default long queue, can also be set in options via -l
    long_queue="gpu2"
    long_queue_time="24:00:00"
elif [ $isBCM -eq 1 ]
then    
    # Juicer directory, contains scripts/, references/, and restriction_sites/
    # can also be set in options via -D
    juiceDir="/home/jgo0012/juicedir"
    # default queue, can also be set in options via -q
    queue="gpu2"
    queue_time="1200"
    # default long queue, can also be set in options via -l
    long_queue="gpu2"
    long_queue_time="3600"
else
    isVoltron=1
    #export PATH=/gpfs0/biobuild/biobuilds-2016.11/bin:$PATH 
    # unset MALLOC_ARENA_MAX # on IBM platform this parameter has significant speed efect but may result in memory leaks
    load_bwa="module load bwa"
    #load_awk="module load awk"
    #load_gpu="spack load cuda@8.0.61 arch=\`spack arch\` && CUDA_VISIBLE_DEVICES=0,1,2,3"
    load_gpu="module load cuda11.0/toolkit"
    #call_gem="/gpfs0/work/neva/gem3-mapper/bin/gem-mapper --3c"
    # Juicer directory, contains scripts/, references/, and restriction_sites/
    # can also be set in options via -D
    juiceDir="/home/jgo0012/juicedir"
    # default queue, can also be set in options
    queue="gpu2"
    queue_time="2880"
    # default long queue, can also be set in options
    long_queue="gpu2"
    long_queue_time="7200"
fi

# size to split fastqs. adjust to match your needs. 4000000=1M reads per split
# can also be changed via the -C flag
splitsize=90000000
#Jamie O: originally was 90,000,000: i doubled this to 180M because appearently this fixes a memory issue that prevents the merged_nodups.txt file from being created.

# fastq files should look like filename_R1.fastq and filename_R2.fastq 
# if your fastq files look different, change this value
read1str="_R1" 
read2str="_R2" 

# unique name for jobs in this run
groupname="a$(date +%s)"

## Default options, overridden by command line arguments

# top level directory, can also be set in options
topDir=$(pwd)
# restriction enzyme, can also be set in options
# default to used in heliconius Hi-C experiment: DpnII
site="DpnII"
# genome ID, default to Heliconius, can also be set in options
genomeID="Pan_dem_hyd_era_chrom"
# description, default empty
about=""
# do not include fragment delimited maps by default
nofrag=1
# use wobble for dedupping by default (not just exact matches)
justexact=0

## Read arguments                                                     
usageHelp="Usage: ${0##*/} [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]\n                 [-a about] [-S stage] [-p chrom.sizes path]\n                 [-y restriction site file] [-z reference genome file]\n                 [-C chunk size] [-D Juicer scripts directory]\n                 [-Q queue time limit] [-L long queue time limit] [-b ligation] [-t threads]\n                 [-A account name] [-e] [-h] [-f] [-j]"
genomeHelp="* [genomeID] must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \n  \"$genomeID\"); alternatively, it can be defined using the -z command"
dirHelp="* [topDir] is the top level directory (default\n  \"$topDir\")\n     [topDir]/fastq must contain the fastq files\n     [topDir]/splits will be created to contain the temporary split files\n     [topDir]/aligned will be created for the final alignment"
queueHelp="* [queue] is the queue for running alignments (default \"$queue\")"
longQueueHelp="* [long queue] is the queue for running longer jobs such as the hic file\n  creation (default \"$long_queue\")"
siteHelp="* [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" \n  (default \"$site\")"
aboutHelp="* [about]: enter description of experiment, enclosed in single quotes"
stageHelp="* [stage]: must be one of \"chimeric\", \"merge\", \"dedup\", \"final\", \"postproc\", or \"early\".\n    -Use \"merge\" when alignment has finished but the merged_sort file has not\n     yet been created.\n    -Use \"dedup\" when the files have been merged into merged_sort but\n     merged_nodups has not yet been created.\n    -Use \"final\" when the reads have been deduped into merged_nodups but the\n     final stats and hic files have not yet been created.\n    -Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use \"early\" for an early exit, before the final creation of the hic files\n    Can also use -e flag to exit early"
pathHelp="* [chrom.sizes path]: enter path for chrom.sizes file"
siteFileHelp="* [restriction site file]: enter path for restriction site file (locations of\n  restriction sites in genome; can be generated with the script\n  misc/generate_site_positions.py)"
chunkHelp="* [chunk size]: number of lines in split files, must be multiple of 4\n  (default ${splitsize}, which equals $(awk -v ss=${splitsize} 'BEGIN{print ss/4000000}') million reads)"
scriptDirHelp="* [Juicer scripts directory]: set the Juicer directory,\n  which should have scripts/ references/ and restriction_sites/ underneath it\n  (default ${juiceDir})"
refSeqHelp="* [reference genome file]: enter path for reference sequence file, BWA index\n  files must be in same directory"
queueTimeHelp="* [queue time limit]: time limit for queue, i.e. -W 12:00 is 12 hours\n  (default ${queue_time})"
longQueueTimeHelp="* [long queue time limit]: time limit for long queue, i.e. -W 168:00 is one week\n  (default ${long_queue_time})"
ligationHelp="* [ligation junction]: use this string when counting ligation junctions"
threadsHelp="* [threads]: number of threads when running BWA alignment"
userHelp="* [account name]: user account name on cluster"
excludeHelp="* -f: include fragment-delimited maps in hic file creation"
justHelp="* -j: just exact duplicates excluded at dedupping step"
earlyexitHelp="* -e: Use for an early exit, before the final creation of the hic files"
gemHelp="* -c: use GEM3 as aligner"
helpHelp="* -h: print this help and exit"


printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$genomeHelp"
    echo -e "$dirHelp"
    echo -e "$queueHelp"
    echo -e "$longQueueHelp"
    echo -e "$siteHelp"
    echo -e "$aboutHelp"
    echo -e "$stageHelp"
    echo -e "$pathHelp"
    echo -e "$siteFileHelp"
    echo -e "$refSeqHelp"
    echo -e "$chunkHelp"
    echo -e "$scriptDirHelp"
    echo -e "$queueTimeHelp"
    echo -e "$longQueueTimeHelp"
    echo -e "$ligationHelp"
    echo -e "$threadsHelp"
    echo -e "$userHelp"
    echo -e "$earlyexitHelp"
    echo -e "$gemHelp"
    echo "$excludeHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:g:a:hq:s:p:l:y:z:S:C:D:Q:L:b:A:t:jfec" opt; do
    case $opt in
    g) genomeID=$OPTARG ;;
    h) printHelpAndExit 0;;
    d) topDir=$OPTARG ;;
    l) long_queue=$OPTARG ;;
    q) queue=$OPTARG ;;
    s) site=$OPTARG ;;
    a) about=$OPTARG ;;
    p) genomePath=$OPTARG ;;  
    y) site_file=$OPTARG ;;
    z) refSeq=$OPTARG ;;
    S) stage=$OPTARG ;;
    C) splitsize=$OPTARG; splitme=1 ;;
    D) juiceDir=$OPTARG ;;
    Q) queue_time=$OPTARG ;;
    L) long_queue_time=$OPTARG ;;
    f) nofrag=0 ;;
    b) ligation=$OPTARG ;;
    t) threads=$OPTARG ;;
    A) user=$OPTARG ;;
    j) justexact=1 ;;
    e) earlyexit=1 ;;
    c) gemmapper=1 ;;
    [?]) printHelpAndExit 1;;
    esac
done

if [ ! -z "$stage" ]
then
    case $stage in
    chimeric) chimeric=1 ;;
        merge) merge=1 ;;
        dedup) dedup=1 ;;
        early) earlyexit=1 ;;
        final) final=1 ;;
    postproc) postproc=1 ;; 
        *)  echo "$usageHelp"
        echo "$stageHelp"
        exit 1
    esac
fi

## Set reference sequence based on genome ID
if [ -z "$refSeq" ]
then 
    case $genomeID in
    mm9)    refSeq="${juiceDir}/references/Mus_musculus_assembly9_norandom.fasta";;
    mm10)   refSeq="${juiceDir}/references/Mus_musculus_assembly10/v0/Mus_musculus_assembly10.fasta";;
    hg38)   refSeq="${juiceDir}/references/hg38/hg38.fa";;
    hg19)   refSeq="${juiceDir}/references/Homo_sapiens_assembly19.fasta";;
    hg18)   refSeq="${juiceDir}/references/hg18.fasta";;
    *)  echo "$usageHelp"
        echo "$genomeHelp"
        exit 1
    esac
else
    ## Reference sequence passed in, so genomePath must be set for the .hic 
    ## file to be properly created
    if [[ -z "$genomePath" ]] && [[ -z $earlyexit ]]
    then
        echo "***! You must define a chrom.sizes file or a standard genome ID via the \"-p\" flag that delineates the lengths of the chromosomes in the genome at $refSeq; you may use \"-p hg19\" or other standard genomes";
        exit 1;
    fi
fi

## Check that refSeq exists 
if [ ! -e "$refSeq" ]; then
    echo "***! Reference sequence $refSeq does not exist";
    exit 1;
fi

## Check that index for refSeq exists
if [[ ! -e "${refSeq}.bwt" ]] && [[ -z $gemmapper ]]
then
    echo "***! Reference sequence $refSeq does not appear to have been indexed. Please run bwa index on this file before running juicer.";
    exit 1;
elif [[ -n $gemmapper ]] && [[ ! -e "${refSeq%.*}.gem" ]]
then
    echo "***! Reference sequence $refSeq does not appear to have been indexed. Please run gem index on this file before running juicer.";
    exit 1;
fi

## Set ligation junction based on restriction enzyme
if [ -z "$ligation" ]; then
    case $site in
    HindIII) ligation="AAGCTAGCTT";;
    MseI)  ligation="TTATAA";;
    DpnII) ligation="GATCGATC";;
    MboI) ligation="GATCGATC";;
        NcoI) ligation="CCATGCATGG";;
    Arima) ligation="'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'" ;;
    none) ligation="XXXX";;
    *)  ligation="XXXX"
        echo "$site not listed as recognized enzyme."
        echo "Ligation junction is undefined"
    esac
fi

if [[ -n $gemmapper ]] 
then
    if [[ "$site" == "none" ]]
    then
    re=""
    else
    re=$(echo $site | tr "+" "\n" |gawk '{str=str" --restriction-enzyme "$1}END{print str}')
    fi
fi

## If DNAse-type experiment, no fragment maps; or way to get around site file
if [[ "$site" == "none" ]] 
then
    nofrag=1;
fi

if [ -z "$site_file" ]
then
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [[ ! -e "$site_file" ]] && [[ "$site" != "none" ]] &&  [[ ! "$site_file" =~ "none" ]]
then
    echo "***! $site_file does not exist. It must be created before running this script."
    exit 1
elif [[ "$site" != "none" ]] && [[ ! "$site_file" =~ "none" ]]
then
    echo  "Using $site_file as site file"
fi

## Set threads for sending appropriate parameters to cluster and string for BWA call
if [ -z "$threads" ]
then
    # default is 8 threads; may need to adjust
    if [ $isRice -eq 1 ]
    then
    threads=8
    threadstring="-t $threads"
    elif [ $isBCM -eq 1 ]
    then
    threads=1
    threadstring="-t $threads"
    else
    threads=8 # VOLTRON; may need to make this separate and specific for Voltron
    ## On voltron with 8 thread per core Power8 CPU bwa can use more threads
    threadstring="-t \$SLURM_JOB_CPUS_PER_NODE"
    fi
else
    if [ $isVoltron -eq 1 ]
    then
    threadstring="-t \$SLURM_JOB_CPUS_PER_NODE"
    else
    threadstring="-t $threads"
    fi
fi

alloc_mem=$(($threads * 8000))

if [ $alloc_mem -gt 80000 ]
then
    alloc_mem=80000
fi

if [ $isBCM -eq 1 ] || [ $isRice -eq 1 ]
then
    alloc_mem=50000
fi

## Directories to be created and regex strings for listing files
splitdir=${topDir}"/splits"
donesplitdir=$topDir"/done_splits"
fastqdir=${topDir}"/fastq/*_R*.fastq*"
outputdir=${topDir}"/aligned"
tmpdir=${topDir}"/scratch/Counterman_Lab_Shared/Hi-C/tmp"
debugdir=${topDir}"/debug"

## Check that fastq directory exists and has proper fastq files
if [ ! -d "$topDir/fastq" ]; then
    echo "Directory \"$topDir/fastq\" does not exist."
    echo "Create \"$topDir/fastq\" and put fastq files to be aligned there."
    echo "Type \"juicer.sh -h\" for help"
    exit 1
else 
    if stat -t ${fastqdir} >/dev/null 2>&1
    then
    echo "(-: Looking for fastq files...fastq files exist"
    else
    if [ ! -d "$splitdir" ]; then 
        echo "***! Failed to find any files matching ${fastqdir}"
        echo "***! Type \"juicer.sh -h \" for help"
        exit 1      
    fi
    fi
fi

## Create output directory, only if not in postproc, dedup or final stages
if [[ -d "$outputdir" && -z "$final" && -z "$dedup" && -z "$postproc" ]] 
then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    echo "***! Type \"juicer.sh -h \" for help"
    exit 1          
else
    if [[ -z "$final" && -z "$dedup" && -z "$postproc" ]]; then
        mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 1; } 
    fi
fi

## Create split directory
if [ -d "$splitdir" ]; then
    splitdirexists=1
else
    mkdir "$splitdir" || { echo "***! Unable to create ${splitdir}, check permissions." ; exit 1; }
fi

## Create temporary directory, used for sort later
if [ ! -d "$tmpdir" ] && [ -z "$final" ] && [ -z "$dedup" ] && [ -z "$postproc" ]; then
    mkdir "$tmpdir"
    chmod 777 "$tmpdir"
fi

## Create output directory, used for reporting commands output
if [ ! -d "$debugdir" ]; then
    mkdir "$debugdir"
    chmod 777 "$debugdir"
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline
# If chunk size sent in, split. Otherwise check size before splitting
if [ $isVoltron -ne 1 ]
then
    if [ -z $splitme ]
    then
    fastqsize=$(ls -lL  ${fastqdir} | gawk '{sum+=$5}END{print sum}')
    if [ "$fastqsize" -gt "2592410750" ]
    then
        splitme=1
    fi
    fi
fi

testname=$(ls -l ${fastqdir} | gawk 'NR==1{print $9}')
if [ "${testname: -3}" == ".gz" ]
then
    read1=${splitdir}"/*${read1str}*.fastq.gz"
    gzipped=1
else
    read1=${splitdir}"/*${read1str}*.fastq"
fi

if [ -z "$user" ]
then
    userstring=""
else
    userstring="#SBATCH -A $user"
fi

# Add header containing command executed and timestamp:
jid=`sbatch <<- HEADER | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l 
        $userstring
    #SBATCH -p $queue
    #SBATCH -t 2
    #SBATCH -c 1
    #SBATCH -o $debugdir/head-%j.out
    #SBATCH -e $debugdir/head-%j.err
    #SBATCH -J "${groupname}_cmd"
    date
    ${load_bwa}
    ${load_java}
    ${load_awk}

    # Experiment description
    if [ -n "${about}" ]
    then
        echo -ne 'Experiment description: ${about}; '
    else
        echo -ne 'Experiment description: '
    fi

    # Get version numbers of all software
    echo -ne "Juicer version $juicer_version;" 
    bwa 2>&1 | gawk '\\\$1=="Version:"{printf(" BWA %s; ", \\\$2)}'
    echo -ne "$threads threads; "
    if [ -n "$splitme" ]
    then
        echo -ne "splitsize $splitsize; "
    fi  
    java -version 2>&1 | gawk 'NR==1{printf("%s; ", \\\$0);}'
    ${juiceDir}/scripts/juicer_tools -V 2>&1 | gawk '\\\$1=="Juicer" && \\\$2=="Tools"{printf("%s; ", \\\$0);}'
    
    echo "$0 $@"
HEADER`
headfile="${debugdir}/head-${jid}.out"

## Record if we failed while aligning, so we don't waste time on other jobs
## Remove file if we're relaunching Juicer 
errorfile=${debugdir}/${groupname}_alignfail
if [ -f $errorfile ]
then
    rm $errorfile
fi

# Not in merge, dedup,  or final stage, i.e. need to split and align files.
if [ -z $merge ] && [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    if [ "$nofrag" -eq 0 ]
    then
    echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $genomeID with site file $site_file"
    else
        echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $genomeID with no fragment delimited maps."
    fi
    
    ## Split fastq files into smaller portions for parallelizing alignment 
    ## Do this by creating a text script file for the job on STDIN and then 
    ## sending it to the cluster    
    dependsplit="afterok"
    if [ ! $splitdirexists ]
    then
    echo "(-: Created $splitdir and $outputdir."
    if [ -n "$splitme" ]
        then
            for i in ${fastqdir}
            do
        filename=$(basename $i)
        filename=${filename%.*}      
                if [ -z "$gzipped" ]
                then    
            jid=`sbatch <<- SPLITEND | egrep -o -e "\b[0-9]+$"
            #!/bin/bash -l
                        #SBATCH -p $queue
            #SBATCH -t $queue_time
            #SBATCH -c 1
            #SBATCH --mem=5G
            #SBATCH -o $debugdir/split-%j.out
            #SBATCH -e $debugdir/split-%j.err
            #SBATCH -J "${groupname}_split_${i}"
                        $userstring         
            date
            echo "Split file: $filename"
            split -a 3 -l $splitsize -d --additional-suffix=.fastq $i $splitdir/$filename
            date
SPLITEND`
        else
            jid=`sbatch <<- SPLITEND | egrep -o -e "\b[0-9]+$"
            #!/bin/bash -l
            #SBATCH -p $queue
            #SBATCH -t $queue_time
            #SBATCH -c 1
            #SBATCH --mem=5G
            #SBATCH -o $debugdir/split-%j.out
            #SBATCH -e $debugdir/split-%j.err
            #SBATCH -J "${groupname}_split_${i}"
                        $userstring         
            date
            echo "Split file: $filename"
            zcat $i | split -a 3 -l $splitsize -d --additional-suffix=.fastq - $splitdir/$filename
            date
SPLITEND`
        fi
        dependsplit="$dependsplit:$jid"
                # if we split files, the splits are named .fastq
                read1=${splitdir}"/*${read1str}*.fastq"
        done
        
        srun -c 1 -p "$queue" -t 1 -o $debugdir/wait-%j.out -e $debugdir/wait-%j.err -d $dependsplit -J "${groupname}_wait" sleep 1
        else
            cp -rs ${fastqdir} ${splitdir}
            wait
        fi
    else
        ## No need to re-split fastqs if they already exist
        echo -e "---  Using already created files in $splitdir\n"
    # unzipped files will have .fastq extension, softlinked gz 
        testname=$(ls -l ${splitdir} | gawk '$9~/fastq$/||$9~/gz$/{print $9; exit}')

        if [[ ${testname: -3} == ".gz" ]]
        then
            read1=${splitdir}"/*${read1str}*.fastq.gz"
        else
        read1=${splitdir}"/*${read1str}*.fastq"
        fi
    fi
    
    ## Launch job. Once split/move is done, set the parameters for the launch. 
    echo "(-: Starting job to launch other jobs once splitting is complete"
    
    ## Loop over all read1/read2 fastq files and create jobs for aligning.
    ## Then call chimeric script on aligned, sort individual
    ## Wait for splits to be individually sorted, then do a big merge sort.
    ## ARRAY holds the names of the jobs as they are submitted
    countjobs=0
    declare -a ARRAY
    declare -a JIDS
    declare -a TOUCH

    dependmerge="afterok"

    for i in ${read1}
    do
    ext=${i#*$read1str}
    name=${i%$read1str*} 
    # these names have to be right or it'll break
    name1=${name}${read1str}
    name2=${name}${read2str}    
    jname=$(basename "$name")${ext}
        usegzip=0
        if [ "${ext: -3}" == ".gz" ]
        then
            usegzip=1
    fi
    touchfile=${tmpdir}/${jname}

    # count ligations
    jid=`sbatch <<- CNTLIG |  egrep -o -e "\b[0-9]+$"
        #!/bin/bash -l
        #SBATCH -p $queue
        #SBATCH -t $queue_time
        #SBATCH -c 1
        #SBATCH -o $debugdir/count_ligation-%j.out
        #SBATCH -e $debugdir/count_ligation-%j.err
        #SBATCH -J "${groupname}_${jname}_Count_Ligation"
        #SBATCH --mem=5G
                $userstring         

        date
        export usegzip=${usegzip}; export name=${name}; export name1=${name1}; export name2=${name2}; export ext=${ext}; export ligation="${ligation}"; ${juiceDir}/scripts/countligations.sh
        date
CNTLIG`
    dependcount="$jid"

    if [ -z "$chimeric" ]
    then
        # align fastqs
        jid=`sbatch <<- ALGNR1 | egrep -o -e "\b[0-9]+$"
        #!/bin/bash -l
        #SBATCH -p $queue
        #SBATCH -o $debugdir/align1-%j.out
        #SBATCH -e $debugdir/align1-%j.err
        #SBATCH -t $queue_time
        #SBATCH -n 1
        #SBATCH -c $threads
        #SBATCH --ntasks=1
        #SBATCH --mem=$alloc_mem
        #SBATCH -J "${groupname}_align1_${jname}"
        #SBATCH --threads-per-core=1        
                $userstring         

        ${load_bwa}

        # Align reads
        date
                if [ \$gemmapper ]; then
            echo "Running command $call_gem $re -I ${refSeq%.*}.gem $threadstring -1 $name1$ext -2 $name2$ext -o $name$ext.sam"
            srun --ntasks=1 $call_gem $re -I ${refSeq%.*}.gem $threadstring -1 $name1$ext -2 $name2$ext -o $name$ext.sam
            if [ \$? -ne 0 ]
            then  
                touch $errorfile
                exit 1
            else
                echo "(-: Gem align of $name$ext.sam done successfully"
            fi
        else
            echo "Running command bwa mem -SP5M $threadstring $refSeq $name1$ext $name2$ext > $name$ext.sam" 
            srun --ntasks=1 bwa mem -SP5M $threadstring $refSeq $name1$ext $name2$ext > $name$ext.sam
            if [ \$? -ne 0 ]
            then  
                touch $errorfile
                exit 1
            else
                echo "(-: Mem align of $name$ext.sam done successfully"
            fi
        fi
        date
ALGNR1`

        dependalign="afterok:$jid:$dependcount"
    else
        dependalign="afterok:$dependcount"
    fi

    if [ $isVoltron -eq 1 ]
    then
        sortthreadstring="--parallel=\$SLURM_JOB_CPUS_PER_NODE"
    else
        sortthreadstring="--parallel=$threads"
    fi

    # wait for alignment, chimeric read handling
    jid=`sbatch <<- MRGALL | egrep -o -e "\b[0-9]+$"
        #!/bin/bash -l
        #SBATCH -p $long_queue
        #SBATCH -o $debugdir/merge-%j.out
        #SBATCH -e $debugdir/merge-%j.err
        #SBATCH --mem=40G
        #SBATCH -t $long_queue_time
        #SBATCH -c 12
        #SBATCH --ntasks=1
        #SBATCH -d $dependalign
        #SBATCH -J "${groupname}_merge_${jname}"
                #SBATCH --threads-per-core=1
                $userstring
        ${load_awk}
        date
        # call chimeric_blacklist.awk to deal with chimeric reads; sorted file is sorted by read name at this point
        touch ${name}${ext}_abnorm.sam ${name}${ext}_unmapped.sam ${name}${ext}_norm.txt
        awk -v "fname1"=${name}${ext}_norm.txt -v "fname2"=${name}${ext}_abnorm.sam -v "fname3"=${name}${ext}_unmapped.sam -f $juiceDir/scripts/chimeric_blacklist.awk ${name}${ext}.sam

        if [ \$? -ne 0 ] 
        then    
            echo "***! Failure during chimera handling of $name${ext}"
            touch $errorfile
            exit 1   
        fi  
        # if any normal reads were written, find what fragment they 
        # correspond to and store that
        # check if site file exists and if so write the fragment number
        # even if nofrag set
        # one is not obligated to provide a site file if nofrag set; 
        # but if one does, frag numbers will be calculated correctly
        if [ -e "$name${ext}_norm.txt" ] && [ "$site" != "none" ] && [ -e "$site_file" ]
        then
            perl ${juiceDir}/scripts/fragment.pl ${name}${ext}_norm.txt ${name}${ext}.frag.txt $site_file
        elif [ "$site" == "none" ] || [ "$nofrag" -eq 1 ]
        then
            awk '{printf("%s %s %s %d %s %s %s %d", \\\$1, \\\$2, \\\$3, 0, \\\$4, \\\$5, \\\$6, 1); for (i=7; i<=NF; i++) {printf(" %s",\\\$i);}printf("\n");}' $name${ext}_norm.txt > $name${ext}.frag.txt
        else
            echo "***! No $name${ext}_norm.txt file created"
            touch $errorfile
            exit 1
        fi
        if [ \$? -ne 0 ]
        then
            echo "***! Failure during fragment assignment of $name${ext}"
            touch $errorfile
            exit 1 
        fi
        # sort by chromosome, fragment, strand, and position
        sort $sortthreadstring -S 35G -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $name${ext}.frag.txt > $name${ext}.sort.txt
        if [ \$? -ne 0 ]   
        then
            echo "***! Failure during sort of $name${ext}"
            touch $errorfile
            exit 1
        else
            rm $name${ext}_norm.txt $name${ext}.frag.txt
        fi
        touch $touchfile
        date
MRGALL`

    dependmerge="${dependmerge}:${jid}"
    ARRAY[countjobs]="${groupname}_merge_${jname}"
    JIDS[countjobs]="${jid}"
    TOUCH[countjobs]="$touchfile"
        countjobs=$(( $countjobs + 1 ))
    done # done looping over all fastq split files
    
    # list of all jobs. print errors if failed    
    for (( i=0; i < $countjobs; i++ ))
    do
    f=${TOUCH[$i]}
    msg="***! Error in job ${ARRAY[$i]}  Type squeue -j ${JIDS[$i]} to see what happened"
    
    # check that alignment finished successfully
    jid=`sbatch <<- EOF
        #!/bin/bash -l
        #SBATCH -o $debugdir/aligncheck-%j.out
        #SBATCH -e $debugdir/aligncheck-%j.err
        #SBATCH -t $queue_time
        #SBATCH -p $queue
        #SBATCH -J "${groupname}_check"
        #SBATCH -d $dependmerge
                $userstring         

        date
        echo "Checking $f"
        if [ ! -e $f ]
        then
            echo $msg
            touch $errorfile
        fi
        date
EOF`
    jid=$(echo $jid | egrep -o -e "\b[0-9]+$")
    dependmergecheck="${dependmerge}:${jid}"
    done
fi  # Not in merge, dedup,  or final stage, i.e. need to split and align files.

# Not in final, dedup, or postproc
if [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    if [ -z $merge ]
    then
    sbatch_wait="#SBATCH -d $dependmergecheck"
    else
        sbatch_wait=""
    fi
    
    # merge the sorted files into one giant file that is also sorted. jid=`sbatch <<- MRGSRT | egrep -o -e "\b[0-9]+$"
    
    if [ $isVoltron -eq 1 ]
    then  
    sbatch_time="#SBATCH -t 10080"
    else
    sbatch_time="#SBATCH -t 1440"
    fi
    if [ $isBCM -eq 1 ]
    then
    sbatch_cpu_alloc="#SBATCH -c 1"
    sbatch_mem_alloc="#SBATCH --mem=80G"
    else
    sbatch_cpu_alloc="#SBATCH -c 8"
    sbatch_mem_alloc="#SBATCH --mem=64G"
    fi


    jid=`sbatch <<- EOF
        #!/bin/bash -l
        #SBATCH -o $debugdir/fragmerge-%j.out
        #SBATCH -e $debugdir/fragmerge-%j.err
        ${sbatch_mem_alloc}
        ${sbatch_time}
        #SBATCH -p $long_queue
        ${sbatch_cpu_alloc}
        #SBATCH -J "${groupname}_fragmerge"
        ${sbatch_wait}
                $userstring         

        date
        if [ -f "${errorfile}" ]
        then
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi
        export LC_COLLATE=C
        if [ -d $donesplitdir ]
        then
            mv $donesplitdir/* $splitdir/.
        fi
        if [ $isRice -eq 1 ]
        then
            if ! ${juiceDir}/scripts/sort --parallel=48 -S 32G -T ${tmpdir} -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/*.sort.txt > $outputdir/merged_sort.txt
            then
                echo "***! Some problems occurred somewhere in creating sorted align files."
                touch $errorfile
                exit 1
            else
                echo "(-: Finished sorting all sorted files into a single merge."
            fi
        else
            if ! sort --parallel=48 -S 32G -T ${tmpdir} -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/*.sort.txt > $outputdir/merged_sort.txt
            then
                echo "***! Some problems occurred somewhere in creating sorted align files."
                touch $errorfile
                exit 1
            else
                echo "(-: Finished sorting all sorted files into a single merge."
            fi
        fi
        date
EOF`

    jid=$(echo $jid | egrep -o -e "\b[0-9]+$")
    dependmrgsrt="afterok:$jid"
fi

# Remove the duplicates from the big sorted file
if [ -z $final ] && [ -z $postproc ]
then
    if [ -z $dedup ]
    then
        sbatch_wait="#SBATCH -d $dependmrgsrt"
    else
        sbatch_wait=""
    fi
    # Guard job for dedup. this job is a placeholder to hold any job submitted after dedup.
    # We keep the ID of this guard, so we can later alter dependencies of inner dedupping phase.
    # After dedup is done, this job will be released. 
    guardjid=`sbatch <<- DEDUPGUARD | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l
    #SBATCH -p $queue
    #SBATCH -o $debugdir/dedupguard-%j.out
    #SBATCH -e $debugdir/dedupguard-%j.err
    #SBATCH -t 10
    #SBATCH -c 1
    #SBATCH -H
    #SBATCH --ntasks=1
    #SBATCH -J "${groupname}_dedup_guard"
    ${sbatch_wait}
        $userstring         

    date
DEDUPGUARD`

    dependguard="afterok:$guardjid"

    # if jobs succeeded, kill the cleanup job, remove the duplicates from the big sorted file
    jid=`sbatch <<- DEDUP | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l
    #SBATCH -p $queue
    #SBATCH --mem-per-cpu=2G
    #SBATCH -o $debugdir/dedup-%j.out
    #SBATCH -e $debugdir/dedup-%j.err
    #SBATCH -t $queue_time
    #SBATCH -c 1
    #SBATCH --ntasks=1
    #SBATCH -J "${groupname}_dedup"
    ${sbatch_wait}
        $userstring
    
    ${load_awk}
    date
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
    squeue -u $USER -o "%A %T %j %E %R" | column -t
    awk -v queue=$long_queue -v groupname=$groupname -v debugdir=$debugdir -v dir=$outputdir -v topDir=$topDir -v juicedir=$juiceDir -v site=$site -v genomeID=$genomeID -v genomePath=$genomePath -v user=$USER -v guardjid=$guardjid -v justexact=$justexact -f $juiceDir/scripts/split_rmdups.awk $outputdir/merged_sort.txt
    ##Schedule new job to run after last dedup part:
    ##Push guard to run after last dedup is completed:
    ##srun --ntasks=1 -c 1 -p "$queue" -t 1 -o ${debugdir}/dedup_requeue-%j.out -e ${debugdir}/dedup-requeue-%j.err -J "$groupname_msplit0" -d singleton echo ID: $ echo "\${!SLURM_JOB_ID}"; scontrol update JobID=$guardjid dependency=afterok:\$SLURM_JOB_ID
    squeue -u $USER -o "%A %T %j %E %R" | column -t
    date
    
    scontrol release $guardjid
DEDUP`

    dependosplit="afterok:$jid"

    #Push dedup guard to run only after dedup is complete:
    scontrol update JobID=$guardjid dependency=afterok:$jid

    #Wait for all parts of split_rmdups to complete:
    jid=`sbatch <<- MSPLITWAIT | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l
    #SBATCH -p $queue
    #SBATCH -o $debugdir/post_dedup-%j.out
    #SBATCH -e $debugdir/post_dedup-%j.err
    #SBATCH -t 100
    #SBATCH -c 1
    #SBATCH --ntasks=1
    #SBATCH -J "${groupname}_post_dedup"
    #SBATCH -d ${dependguard}
        $userstring         

    date
    rm -Rf $tmpdir;
    find $debugdir -type f -size 0 | xargs rm
    squeue -u $USER -o "%A %T %j %E %R" | column -t
    date
MSPLITWAIT`

    dependmsplit="afterok:$jid"
    sbatch_wait="#SBATCH -d $dependmsplit"
else
    sbatch_wait=""
fi

if [ -z "$genomePath" ]
then
    #If no path to genome is give, use genome ID as default.
    genomePath=$genomeID
fi

#Skip if post-processing only is required
if [ -z $postproc ]
    then
    # Check that dedupping worked properly
    # in ideal world, we would check this in split_rmdups and not remove before we know they are correct
    awkscript='BEGIN{sscriptname = sprintf("%s/.%s_rmsplit.slurm", debugdir, groupname);}NR==1{if (NF == 2 && $1 == $2 ){print "Sorted and dups/no dups files add up"; printf("#!/bin/bash -l\n#SBATCH -o %s/dup-rm.out\n#SBATCH -e %s/dup-rm.err\n#SBATCH -p %s\n#SBATCH -J %s_msplit0\n#SBATCH -d singleton\n#SBATCH -t 1440\n#SBATCH -c 1\n#SBATCH --ntasks=1\ndate;\nrm %s/*_msplit*_optdups.txt; rm %s/*_msplit*_dups.txt; rm %s/*_msplit*_merged_nodups.txt;rm %s/split*;\ndate\n", debugdir, debugdir, queue, groupname, dir, dir, dir, dir) > sscriptname; sysstring = sprintf("sbatch %s", sscriptname); system(sysstring);close(sscriptname); }else{print "Problem"; print "***! Error! The sorted file and dups/no dups files do not add up, or were empty."}}'
    jid=`sbatch <<- DUPCHECK | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l
    #SBATCH -p $queue
    #SBATCH -o $debugdir/dupcheck-%j.out
    #SBATCH -e $debugdir/dupcheck-%j.err
    #SBATCH -t $queue_time
    #SBATCH -c 1
    #SBATCH --ntasks=1
    #SBATCH --mem-per-cpu=1G
    #SBATCH -J "${groupname}_dupcheck"
    ${sbatch_wait}
        $userstring         
    #${load_awk}
    date      
    ls -l ${outputdir}/merged_sort.txt | gawk '{printf("%s ", \\\$5)}' > $debugdir/dupcheck-${groupname}
    ls -l ${outputdir}/merged_nodups.txt ${outputdir}/dups.txt ${outputdir}/opt_dups.txt | gawk '{sum = sum + \\\$5}END{print sum}' >> $debugdir/dupcheck-${groupname}
    gawk -v debugdir=$debugdir -v queue=$queue -v groupname=$groupname -v dir=$outputdir '$awkscript' $debugdir/dupcheck-${groupname}
        date                                                                                                           
DUPCHECK`
    sbatch_wait="#SBATCH -d afterok:$jid"

    jid=`sbatch <<- PRESTATS | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l
    #SBATCH -p $queue
    #SBATCH -o $debugdir/prestats-%j.out
    #SBATCH -e $debugdir/prestats-%j.err
    #SBATCH -t $queue_time
    #SBATCH -c 1
    #SBATCH --ntasks=1
    #SBATCH --mem-per-cpu=1G
    #SBATCH -J "${groupname}_prestats"
    ${sbatch_wait}
    $userstring
    #${load_awk}
        date
        ${load_java}
        export IBM_JAVA_OPTIONS="-Xmx1024m -Xgcthreads1"                                                                                                         
        export _JAVA_OPTIONS="-Xmx1024m -Xms1024m"                                                                                                              
        
        tail -n1 $headfile | gawk '{printf"%-1000s\n", \\\$0}' > $outputdir/inter.txt                                                                              
        cat $splitdir/*.res.txt | gawk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter.txt                                                                
        ${juiceDir}/scripts/juicer_tools LibraryComplexity $outputdir inter.txt >> $outputdir/inter.txt                                                           
        cp $outputdir/inter.txt $outputdir/inter_30.txt                                                       
        date
PRESTATS`

    sbatch_wait0="#SBATCH -d afterok:$jid"
    jid=`sbatch <<- STATS | egrep -o -e "\b[0-9]+$"
        #!/bin/bash -l
        #SBATCH -p $long_queue
        #SBATCH -o $debugdir/stats-%j.out
        #SBATCH -e $debugdir/stats-%j.err
        #SBATCH -t $long_queue_time
        #SBATCH -c 1
        #SBATCH --ntasks=1
        #SBATCH --mem=25G
        #SBATCH -J "${groupname}_stats"
        ${sbatch_wait0}
                $userstring         

        date
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
        perl ${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt

        date
STATS`
    sbatch_wait1="#SBATCH -d afterok:$jid"

    dependstats="afterok:$jid"
    jid=`sbatch <<- STATS30 | egrep -o -e "\b[0-9]+$"
        #!/bin/bash -l
        #SBATCH -p $long_queue
        #SBATCH -o $debugdir/stats30-%j.out
        #SBATCH -e $debugdir/stats30-%j.err
        #SBATCH -t $long_queue_time
        #SBATCH -c 1
        #SBATCH --ntasks=1
        #SBATCH --mem=25G
        #SBATCH -J "${groupname}_stats"
        ${sbatch_wait0}
                $userstring         

        perl ${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter_30.txt -q 30 $outputdir/merged_nodups.txt
        date
STATS30`

    dependstats30="afterok:$jid"
    sbatch_wait1="${sbatch_wait1}:$jid"
    # This job is waiting on deduping, thus sbatch_wait (vs sbatch_wait0 or 1) 
    jid=`sbatch <<- CONCATFILES | egrep -o -e "\b[0-9]+$"
        #!/bin/bash -l
        #SBATCH -p $long_queue
        #SBATCH -o $debugdir/stats30-%j.out
        #SBATCH -e $debugdir/stats30-%j.err
        #SBATCH -t $long_queue_time
        #SBATCH -c 1
        #SBATCH --ntasks=1
        #SBATCH --mem=25G
        #SBATCH -J "${groupname}_stats"
        ${sbatch_wait}
                $userstring         
        # ${load_awk}
        cat $splitdir/*_abnorm.sam > $outputdir/abnormal.sam
        cat $splitdir/*_unmapped.sam > $outputdir/unmapped.sam
        # collect collisions and dedup them
        gawk -f ${juiceDir}/scripts/collisions.awk $outputdir/abnormal.sam > $outputdir/collisions.txt
        # dedup: two pass algorithm, ideally would make one pass
        gawk -v fname=$outputdir/collisions.txt -f ${juiceDir}/scripts/collisions_dedup_rearrange_cols.awk $outputdir/collisions.txt | sort -k3,3n -k4,4n -k10,10n -k11,11n -k17,17n -k18,18n -k24,24n -k25,25n -k31,31n -k32,32n | gawk -v name=$outputdir/ -f ${juiceDir}/scripts/collisions_dups.awk
        date
CONCATFILES`

    # if early exit, we stop here, once the stats are calculated
    if [ ! -z "$earlyexit" ]
    then
    jid=`sbatch <<- FINCLN1 | egrep -o -e "\b[0-9]+$" 
    #!/bin/bash -l
    #SBATCH -p $queue
    #SBATCH --mem=2G
    #SBATCH -o $debugdir/fincln1-%j.out
    #SBATCH -e $debugdir/fincln1-%j.err
    #SBATCH -t 1200
    #SBATCH -c 1
    #SBATCH --ntasks=1
    #SBATCH -J "${groupname}_prep_done"     
        #SBATCH --mail-type=END,FAIL
    ${sbatch_wait1}
        $userstring    


    date
    export splitdir=${splitdir}; export outputdir=${outputdir}; export early=1; ${juiceDir}/scripts/check.sh
    date
FINCLN1`
    echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id $jid"
    exit 0
    fi
    
    jid=`sbatch <<- HIC | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l
    #SBATCH -p $long_queue
    #SBATCH -o $debugdir/hic-%j.out
    #SBATCH -e $debugdir/hic-%j.err 
    #SBATCH -t $long_queue_time
    #SBATCH -c 1
    #SBATCH --ntasks=1
    #SBATCH --mem=49G
    #SBATCH -J "${groupname}_hic"
    #SBATCH -d $dependstats
        $userstring         

    ${load_java}
    export IBM_JAVA_OPTIONS="-Xmx49152m -Xgcthreads1"
    export _JAVA_OPTIONS="-Xmx49152m -Xms49152m"
    date
    if [ -f "${errorfile}" ]
    then 
        echo "***! Found errorfile. Exiting." 
        exit 1 
    fi 
    
    if [ "$nofrag" -eq 1 ]
    then 
        ${juiceDir}/scripts/juicer_tools pre -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath
    else
        ${juiceDir}/scripts/juicer_tools pre -f $site_file -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath
    fi
    date
HIC`

    dependhic="afterok:$jid"

    jid=`sbatch <<- HIC30 | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l
    #SBATCH -p $long_queue
    #SBATCH -o $debugdir/hic30-%j.out
    #SBATCH -e $debugdir/hic30-%j.err
    #SBATCH -t $long_queue_time
    #SBATCH -c 1
    #SBATCH --ntasks=1
    #SBATCH --mem=49G
    #SBATCH -J "${groupname}_hic30"
    #SBATCH -d ${dependstats30}
        $userstring         

    ${load_java}
    export IBM_JAVA_OPTIONS="-Xmx49152m -Xgcthreads1"
    export _JAVA_OPTIONS="-Xmx49152m -Xms49152m"
    date
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
        if [ "$nofrag" -eq 1 ]
        then 
        ${juiceDir}/scripts/juicer_tools pre -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath
    else
        ${juiceDir}/scripts/juicer_tools pre -f $site_file -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath
    fi
    date
HIC30`

    dependhic30="${dependhic}:$jid"
    sbatch_wait="#SBATCH -d $dependhic30"
else
    sbatch_wait=""
fi

if [[ "$isNots" -eq 1 ]] || [[ "$isVoltron" -eq 1 ]]
then
    if [[  "$isNots" -eq 1 ]] 
    then
    sbatch_req="#SBATCH --gres=gpu:kepler:1"
    fi
    jid=`sbatch <<- HICCUPS | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l
    #SBATCH -p $queue
    #SBATCH --mem-per-cpu=4G
    ${sbatch_req}
    #SBATCH -o $debugdir/hiccups_wrap-%j.out
    #SBATCH -e $debugdir/hiccups_wrap-%j.err
    #SBATCH -t $queue_time
    #SBATCH --ntasks=1
    #SBATCH -J "${groupname}_hiccups_wrap"
    ${sbatch_wait}
        $userstring

    ${load_gpu}
    echo "load: $load_gpu"
    ${load_java}
    date
    nvcc -V
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
    /usr/bin/java -Xmx40g -Djava.io.tmpdir=/scratch/Counterman_Lab_Shared/Hi-C/tmp -jar ${juiceDir}/scripts/juicer_tools.jar hiccups --cpu ${outputdir}/inter.hic ${outputdir}/inter_loops/ --ignore_sparsity
    date
HICCUPS`
    dependhiccups="afterok:$jid"
else
    dependhiccups="afterok"
fi

jid=`sbatch <<- ARROWS | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l
    #SBATCH -p $queue
    #SBATCH --mem-per-cpu=8G
    #SBATCH -o $debugdir/arrowhead_wrap-%j.out
    #SBATCH -e $debugdir/arrowhead_wrap-%j.err
    #SBATCH -t $queue_time
    #SBATCH --ntasks=1
    #SBATCH -J "${groupname}_arrowhead_wrap"
    ${sbatch_wait}
        $userstring         

    ${load_java}
    date
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
    /usr/bin/java -Xmx40g -Djava.io.tmpdir=/scratch/Counterman_Lab_Shared/Hi-C/tmp -jar ${juiceDir}/scripts/juicer_tools.jar arrowhead --cpu ${outputdir}/inter.hic ${outputdir}/inter_contact_domains/ --ignore_sparsity
    date;
ARROWS`
dependarrows="${dependhiccups}:$jid"

jid=`sbatch <<- FINCLN1 | egrep -o -e "\b[0-9]+$"
    #!/bin/bash -l
    #SBATCH -p $queue
    #SBATCH --mem-per-cpu=2G
    #SBATCH -o $debugdir/fincln-%j.out
    #SBATCH -e $debugdir/fincln-%j.err
    #SBATCH -t 1200
    #SBATCH -c 1
    #SBATCH --ntasks=1
    #SBATCH -J "${groupname}_prep_done"
    #SBATCH -d $dependarrows
        $userstring         

    date
    export splitdir=${splitdir}; export outputdir=${outputdir}; ${juiceDir}/scripts/check.sh
    date
FINCLN1`

echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id $jid"



