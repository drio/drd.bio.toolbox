#!/bin/bash
#
source "`dirname ${BASH_SOURCE[0]}`/common.sh"

usage()
{
  cat<<EOF
Usage: 
  $0 <cmd to run> <job_name> [number of cores]
EOF
}

error()
{
  echo "Ups!: $1."
  usage
  exit 1
}

log()
{
  echo "#`date`: $1"
}

dump_pbs_script()
{
  (
  cat <<-EOF
#!/bin/bash
cd `pwd`
$1
EOF
) > pbs.$2.sh
}

# Main
#
[ ".$1" == "." ] && error "Need cmd command"
cmd=$1
log "cmd: $cmd"

[ ".$2" == "." ] && error "Need job name"
j_name=$2 
log "job name: $j_name"

if [ ".$3" == "." ]
then 
  log "Using 1 core"
  n_cores=1
else
  log "Using $3 core(s)"
  n_cores=$3
fi

dump_pbs_script "$cmd" $j_name
chmod 755 ./pbs."$j_name".sh
echo "$submit_bin -N \"$j_name\" -q '$default_queue' -l \"nodes=1:ppn=$n_cores\" \"./pbs.$j_name.sh\""
