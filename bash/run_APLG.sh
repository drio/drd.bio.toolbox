#!/bin/bash
#
# Author: David Rio Deiros
#
usage()
{
cat << EOF
Usage: `basename $0` options

This scripts generates code to run Broad's AllPathLG assembler against
a set of nextgen data.

OPTIONS:
   -h      Show this message
   -i      Ploidy of the genome [$PLOIDY]
   -d      Input data quality is in Phread_64 [$PHRED_64]
   -m      Max memory to use in assembly [$ASSEMBLY_MAX_MEM]
   -c      Minimun size to consider a contig [$MIN_CONTIG]
   -t      Max number of threads in assembly [$ASSEMBLY_N_THREADS]
   -s      Specie name (Mandatory)
EOF
}

error()
{
  local func=$1
  local msg=$2
  echo "ERROR: ${func}() ${msg}" >&2
  exit 1
}

check_requirements()
{
  local func=check_requirements
  # Let's make sure we have the necessary input files
  [ ! -f ${CSV_GROUPS} ] && error "${func}" "I cannot find ${CSV_GROUPS}"
  [ ! -f ${CSV_LIBS}   ] && error "${func}" "I cannot find ${CSV_LIBS}"
  # Let's make sure we have the software we need
  PATH=$PATH:${AP_PATH}
  for b in RunAllPathsLG PrepareAllPathsInputs.pl
  do
    hash ${b} || error "${func}" "I cannot find: -${b}-"
  done
}

check_inputs()
{
  if [[ -z $SPECIE ]]
  then
    usage
    exit 1
  fi
}

set_defaults()
{
  DEF_PLOIDY=2             ; PLOIDY=$DEF_PLOIDY
  DEF_PHRED_64=0           ; PHRED_64=$DEF_PHRED_64
  DEF_MIN_CONTIG=300       ; MIN_CONTIG=$DEF_MIN_CONTIG
  DEF_ASSEMBLY_MAX_MEM=50  ; ASSEMBLY_MAX_MEM=$DEF_ASSEMBLY_MAX_MEM
  DEF_ASSEMBLY_N_THREADS=2 ; ASSEMBLY_N_THREADS=$DEF_ASSEMBLY_N_THREADS
  DEF_SPECIE=""            ; SPECIE=$DEF_SPECIE

  #AP_PATH="/stornext/snfs6/rogers/drio_scratch/tmp/allpaths/local/bin"
  AP_PATH="/users/yl131317/work_dir/allpathslg-39784/bin"
  #PICARD_TOOLS_DIR="/stornext/snfs6/rogers/drio_scratch/bb/local/picard-tools-1.60"
  PICARD_TOOLS_DIR="/stornext/snfs6/assembly/local-bin/picard-tools-1.38"
  COMMON_CMD="source $HOME/.bashrc; ulimit -s 100000"
  CSV_GROUPS=in_groups.csv
  CSV_LIBS=in_libs.csv
}

process_arguments()
{
  while getopts “hi:d:m:c:s:t:” OPTION
  do
    case $OPTION in
      h) usage; exit 1              ;;
      i) PLOIDY=$OPTARG             ;;
      d) PHRED_64=$OPTARG           ;;
      m) ASSEMBLY_MAX_MEM=$OPTARG   ;;
      c) MIN_CONTIG=$OPTARG         ;;
      s) SPECIE=$OPTARG             ;;
      t) ASSEMBLY_N_THREADS=$OPTARG ;;
      ?) usage; exit                ;;
    esac
  done
}

set_cmds()
{
  # Rationale of the APLG directory structure:
  # ${PWD}/REFERENCE_NAME/DATA/RUN/ASSEMBLIES/SUBDIR
  # REFERENCE_NAME: the dir name for the project
  # DATA          : converted reads will live here
  # RUN           : intermediate files for the assembly
  # ASSEMBLIES    : it is hardcoded to assemblies you cannot control it
  # SUBDIR        : Your attempt of assembly will be here
  #
  OUTPUT_DIR="APLG.assembly.${SPECIE}" # that's the root dir for the assembly
  PREPARE_CMD=" ${COMMON_CMD}; \
  PrepareAllPathsInputs.pl \
   PICARD_TOOLS_DIR=${PICARD_TOOLS_DIR} \
   DATA_DIR=${PWD}/${OUTPUT_DIR}/data \
   PLOIDY=${PLOIDY} \
   IN_GROUPS_CSV=in_groups.csv \
   IN_LIBS_CSV=in_libs.csv \
   PHRED_64=$PHRED_64 \
   TMP_DIR=/space1/tmp \
   OVERWRITE=True"
  #
  ASSEMBLY_CMD=" ${COMMON_CMD}; \
  RunAllPathsLG \
   PRE=$PWD \
   REFERENCE_NAME=${OUTPUT_DIR} \
   DATA_SUBDIR=data \
   RUN=run \
   MAX_MEMORY_GB=${ASSEMBLY_MAX_MEM} \
   MIN_CONTIG=300 \
   THREADS=${ASSEMBLY_N_THREADS} \
   SUBDIR=data \
   TARGETS=standard \
   OVERWRITE=True"
}

dump_script()
{
  PBS_LOGS=${PWD}/logs
  PBS_QUEUE_PREP=gac
  PBS_QUEUE_ASSEMBLY=gac
  PBS_PREP_RESOURCES="nodes=1:ppn=3,mem=100000mb"
  PBS_ASSEMBLY_RESOURCES="nodes=1:ppn=${ASSEMBLY_N_THREADS},mem=${ASSEMBLY_MAX_MEM}G,feature=bigmem"
  PBS_JOB_NAME_SEED="APLG.${SPECIE}"

  mkdir -p ${PWD}/${OUTPUT_DIR}/data
  mkdir -p ${PBS_LOGS}
  RAND="${RANDOM}_${RANDOM}"
  prep_job_name="APLG_prepare_${SPECIE}_${RAND}"
  ass_job_name="APLG_assembly_${SPECIE}_${RAND}"

  cat<<-EOF
#!/bin/bash
#
# Prepare data step
#
id=\`echo "(${PREPARE_CMD})" | qsub \
  -N "${prep_job_name}" \
  -q "${PBS_QUEUE_PREP}" \
  -d "${PWD}" \
  -e ${PBS_LOGS}/${prep_job_name}.o -o ${PBS_LOGS}/${prep_job_name}.e \
  -l ${PBS_PREP_RESOURCES} -V\`

#
# Actuall assembly
#
echo "(${ASSEMBLY_CMD})" | qsub \
  -W "depend=afterok:\$id" \
  -N "${ass_job_name}" \
  -q "${PBS_QUEUE_ASSEMBLY}" \
  -d "${PWD}" \
  -e ${PBS_LOGS}/${ass_job_name}.o -o ${PBS_LOGS}/${ass_job_name}.e \
  -l ${PBS_ASSEMBLY_RESOURCES} -V
EOF
}

# Main
#
set_defaults
process_arguments $*
check_requirements
check_inputs
set_cmds
dump_script
