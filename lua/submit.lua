#!/usr/bin/env lua

-- Submit jobs to a moab cluster

require 'klib'

function error(defs, msg)
  if msg then print("ERROR: " .. msg) end
  print("Usage: submit -s job_name [<options> <cmd>]")
  print([[
  -m memory (]] .. defs.m .. [[)
  -c number of cores (]] .. defs.c .. [[)
  -q queue to use (]] .. defs.q .. [[)
  -l directory to use to dump logs (]] .. defs.l .. [[)
  -d list of dependencies. dep1:dep2:.... :depn ("")
  -f file to locate list of dependencies. 1 line one dep.
  -e dump in the stderr the cmd string

  <cmd> is the set of shell commands to execute
  If no <cmd> provided, we will read from stdin

Examples:
  # Iterate over a bunch of files and submit a jobs using each file
  $ ls *.fastq | xargs -i submit -s fmi.{} -m 8G -c 8 "sga index --no-reverse -d 5000000 -t 8 {}"

  # Similar to before but using shell's for
  $ for i in `ls ../input/reads.*.fastq`; do F=`basename $i .fastq`; submit -s pp.$F "sga preprocess -o $F.pp.fastq --pe-mode 2 $i"; done

  # Submit a two jobs, the second one has to run after the first one completes
  $ submit -s one "touch ./one.txt" | bash > /tmp/deps.txt ; submit -s two -f /tmp/deps.txt  "sleep 2;touch ./two.txt" | bash

  # Same as before but now we specify the jobid in the command line instead in a file
  $ submit -s filter -d 3678650.sug-moab -m 20G -c 6 "sga fm-merge -m 65 -t 6 final.filter.pass.fa"

  # And my favourite one.
  # The second command reads the jobid from the standard input and uses it as dep
  $ (submit -s one "sleep 15; touch ./one.txt" | bash) | submit -s two -f -  "sleep 2;touch ./two.txt" | bash
]])
  os.exit(1)
end

-- Load defaults
local defs = {}
defs.m = "4Gb"       -- mem
defs.c = "1"         -- cores
defs.q = "analysis"  -- queue
defs.l = "moab.logs" -- dir to dump logs
defs.d = ""          -- no deps by default
defs.e = false       -- by default, do not print cmd to stderr
local opts = {}
for k,v in pairs(defs) do opts[k] = v end

-- Process arguments
for o,a in os.getopt(arg, 'm:c:q:s:l:d:f:e') do
  opts[o] = a
  if o == 'e' then opts.e = true end
end
if not opts.s then error(defs,"I need a name for the job.") end

-- if -f, load deps in opts.d
if opts.f then
  local _ = {}
  local file = assert(io.xopen(opts.f), "Cannot find -f file: " .. opts.f)
  for line in file:lines() do table.insert(_,line) end
  opts.d = table.concat(_,':')
end

-- Get the cmd from arg ; data may be coming from stdin
local cmd = table.concat(arg, " ")
cmd = (cmd == "") and io.read("*all") or cmd

-- Dump the moab command
local output = {}
local ti = table.insert
ti(output, 'echo \'' .. cmd .. '\'')
ti(output, ' |')
ti(output, ' qsub -N ' .. opts.s)
ti(output, ' -W depend=afterok:' .. opts.d)
ti(output, ' -q ' .. opts.q)
ti(output, ' -d `pwd`')
ti(output, ' -o moab_logs/' .. opts.s .. ".o")
ti(output, ' -e moab_logs/' .. opts.s .. ".e")
ti(output, ' -l nodes=1:ppn=' .. opts.c .. ",mem=" .. opts.m)
ti(output, ' -V')

print("mkdir -p moab_logs;", table.concat(output, ''))

-- Dump to stderr the cmd if necessary
if opts.e then io.stderr:write(cmd .. "\n") end
