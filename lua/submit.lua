#!/usr/bin/env lua

--[[
  Submit jobs to a moab cluster
  Ex:
  $ submit -s test "sga preprocess -o $F.pp.fastq --pe-mode 2 $i"
  $ submit -s test2 -m 8G -c 8 "sga index --no-reverse -d 5000000 -t 8 lala.fa"
  $ submit -s test3 -m 8G -c 8 -q "sga index --no-reverse -d 5000000 -t 8 lala.fa"
--]]

require 'klib'

function error(defs, msg)
  if msg then print("ERROR: " .. msg) end
  print("Usage: " .. arg[0] .. " -s job_name [<options> <cmd>]")
  print([[
  -m memory (]] .. defs.m .. [[)
  -c number of cores (]] .. defs.c .. [[)
  -q queue to use (]] .. defs.q .. [[)
  -l directory to use to dump logs (]] .. defs.l .. [[)
  <cmd> is the set of shell commands to execute
  If no <cmd> provided, we will read from stdin
]])
  os.exit(1)
end

-- Load defaults
local defs = {}
defs.m = "4Gb"       -- mem
defs.c = "1"         -- cores
defs.q = "analysis"  -- queue
defs.l = "moab.logs" -- dir to dump logs
local opts = {}
for k,v in pairs(defs) do opts[k] = v end

-- Process arguments
for o,a in os.getopt(arg, 'm:c:q:s:l:') do opts[o] = a end
if not opts.s then error(defs,"I need a name for the job.") end

-- Get the cmd from arg ; data may be coming from stdin
local cmd = table.concat(arg, " ")
cmd = (cmd == "") and io.read("*all") or cmd

-- Dump the moab command
local output = {}
local ti = table.insert
ti(output, 'echo "' .. cmd .. '"')
ti(output, ' |')
ti(output, ' qsub -N ' .. opts.s)
ti(output, ' -q ' .. opts.q)
ti(output, ' -d `pwd`')
ti(output, ' -o moab_logs/' .. opts.s .. ".o")
ti(output, ' -e moab_logs/' .. opts.s .. ".e")
ti(output, ' -l nodes=1:ppn=' .. opts.c .. ",mem=" .. opts.m)
ti(output, ' -V')

print("mkdir -p moab_logs;", table.concat(output, ''))
