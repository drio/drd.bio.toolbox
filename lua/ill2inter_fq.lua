#!/usr/bin/env luajit
--
-- Given two files: read1 and read2 in fastq format, generate
-- and interlaced version of them.
-- Typical output will be sequence files from illumina
--
require 'klib'

local usage = function()
  io.stderr:write("Usage: ")
  io.stderr:write(arg[0] .. " -o <read1.fq> -t <read2.fq>" .. "\n")
  os.exit(1)
end

local r1_fn, r2_fn, r1, r2; -- fn for input files, streams for input files
for o, a in os.getopt(arg, 'o:t:') do
  if o == 'o' then r1_fn = a end
  if o == 't' then r2_fn = a end
end

-- Sanity checks
if not r1_fn or not r2_fn then usage() end

local r1 = io.xopen(r1_fn)
local r2 = io.xopen(r2_fn)

local i, e_one, e_two, index, line;
i = 0
e_one = {}; e_two = {}
while true do
  i = i + 1
  e_one[i] = r1:read("*line")
  e_two[i] = r2:read("*line")
  if not e_one[i] or not e_two[i] then break end
  if i % 4 == 0 then
    for index, line in ipairs(e_one) do print(e_one[index]) end
    for index, line in ipairs(e_one) do print(e_two[index]) end
    i = 0;
  end
end

r1:close()
r2:close()

