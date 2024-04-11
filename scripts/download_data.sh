# Author: Marco Kreuzer <marco.kreuzer@bioinformatics.unibe.ch>
#/
#/ This script downloads data from the NGSP server to a project folder and checks
#/ the md5sum of the files to ensure that the data transfer was done correctly.
#/
help="
#/ Usage: download_data.sh -d <run_folder> -o <output_folder> -l <log_folder>;
#/ ;
#/ OPTIONS;
#/;
#/  -d	Absolute path to /opt/illumina/output/<run_folder>/Unaligned_./<project> folder on ;
#/  -o  Relative path to where the data files should be stored;
#/  -l  path to log folder. The name of the file cannot be changed and is called;
#/      download.log;
#/;
#/ EXAMPLES;
#/ ;
#/ This is the recommended way for the structure with IBU projects:;
#/ scripts/download_data.sh \;
#/		-d /opt/illumina/output/220726_FS10001441_0071_BPL20311-1822/Unaligned_1/Seehausen_Pooja_OS \;
#/		-o raw_data;
#/      -l log;
#/;
#/ OUTPUT: ;
#/   raw_data/*multiqc.txt: Contains the multiqc report;
#/   raw_data/*fastq.gz: copied fastq files;
#/   raw_data/md5sums_{project}.txt: the original md5sums;
#/   raw_data/md5sums_ms01.txt: the newly calculated md5sums of the copied files.;
#/   log/download.log: The log specifies if the md5sums are equal or not.;
"

#{{{ CL arguments


function print_help { echo $help | sed 's/;/\n/g';} 

while getopts ":d:o:l:" opt; do
  case ${opt} in
    d )
      directory=$OPTARG
      ;;
    o )
      outdir=$OPTARG
      ;;
    l )
      log=$OPTARG
      ;;
    ? )
      print_help; exit 2;
      ;;
	\? | h | *)
      print_help; exit 2;
      ;;
  esac
done
shift $((OPTIND -1))
#}}}

# Bash settings

set -o errexit # abort on nonzero exitstatus
set -o nounset # abort on unbound variable     
set -o pipefail # dont hide errors within pipes


#{{{ Variales
readonly script_name=$(basename "${0}")
readonly script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
IFS=$'\t\n'   # Split on newlines and tabs (but not on spaces)
#}}}

main() {
	 if [ $(hostname) = "login8.hpc.binf.unibe.ch" ]; then
	 	printf '%s\n' "You are on the head node $(hostname)."
		printf '%s\n' "Please log onto a node for this operation."
		exit 1;
	 fi
	 run_folder=$(dirname $(dirname ${directory}));
	 project=$(basename ${directory});
	 remote='ibuillumina.binf.unibe.ch'
	 log="${log}/download.log";
	 echo ${directory} 2>&1 | tee -a ${log}; 
	 ssh $remote "ls ${directory}" > ${outdir}/files.txt;
	 cat ${outdir}/files.txt | xargs -n 1 -P 10 \
			 -I {} rsync \
			 -rav ${remote}:${directory}/{} ${outdir}/. ;
	 echo -e "\nCopying of fastq files done. " 2>&1 | tee -a ${log};
	 #
	 rsync -rav  ${remote}:${run_folder}/${project}/analysis/illumina-qc/md5sums/md5sums.txt  ${outdir}/.  ;
	 echo -e "\nCopying of md5sums file done. " 2>&1 | tee -a ${log};
	 rsync -rav  ${remote}:${run_folder}/${project}/analysis/illumina-qc/multiqc/multiqc.html  ${outdir}/.  ;
	 echo -e "\nCopying of multiQC report done. " 2>&1 | tee -a ${log};
	 
	 echo -e "\nCalculating md5sums of copied files. " 2>&1 | tee -a ${log};
	 (cd ${outdir}; 
	 	md5sum -c md5sums.txt 2>&1 | tee -a ../${log} );
}

main "${@}"
