#!/usr/bin/env python
#
import sys, re

class Filter:
  PASS = "PASS"
  def md_string(self):
    return '##FILTER=<ID=%(id)s,Description="%(desc)s">' % {
      'id': self.fid(), 'desc': self.fdesc()
    }

  def passes(self, line):
    if self.logic_passes(line.split()):
      return self.PASS
    else:
      return self.fid()

class FilterQual(Filter):
  def __init__(self, q):
    self.min_qual = q

  def logic_passes(self, sp):
    col_qual = 5
    return float(sp[col_qual]) >= self.min_qual

  def fdesc(self):
    return "quality below " +  str(self.min_qual)

  def fid(self):
    return "q" + str(self.min_qual)

class FilterCov(Filter):
  COV_RE = re.compile(r"DP=(\d+);")
  def __init__(self, _min, _max):
    self._min, self._max = _min, _max

  def logic_passes(self, line):
    match = self.COV_RE.search('\t'.join(line))
    if match:
      dp = int(match.group(1))
      return dp >= self._min and dp <= self._max
    else:
      raise 'No RDP field in entry.'

  def fdesc(self):
    return "Coverage not within desired margins min:%(min)s max:%(max)s" % {
      'min': self._min, 'max': self._max
    }

  def fid(self):
    return "c(%(min)s|%(max)s)" % { 'min': self._min, 'max': self._max }

class FilterNumSamples(Filter):
  SF_RE = re.compile(r"SF=([\d,]+)")
  def __init__(self, _min):
    self.min_num_samples = _min

  def logic_passes(self, line):
    match = self.SF_RE.search('\t'.join(line))
    if match:
      num_sfs = len(match.group(1).split(","))
      return num_sfs >= self.min_num_samples
    else:
      raise 'Cannot match against the SF field.'

  def fdesc(self):
    return "Minimum number of samples with event: " +  str(self.min_num_samples)

  def fid(self):
    return "mns" + str(self.min_num_samples)

class FilterRawCoverage(Filter):
  SF_RDP = re.compile(r"RDP=([\d,]+)")
  def __init__(self, mc):
    self.min_cov = mc

  def logic_passes(self, line):
    match = self.SF_RDP.search('\t'.join(line))
    if match:
      rcov = match.group(1).split(",")
      for r in rcov:
        if int(r) < self.min_cov:
          return False
      return True
    else:
      return False
      #raise 'Cannot match against the RDP field.'

  def fdesc(self):
    return "Minimum raw coverage per sample: " +  str(self.min_cov)

  def fid(self):
    return "rcov" + str(self.min_cov)

def test():
  MD_QUAL = '##FILTER=<ID=q10,Description="quality below 10">'
  low_qual="20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3"
  ok_qual="20     17330   .         T      A       123    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3"

  f = FilterQual(10)
  assert(f.md_string() == MD_QUAL)
  assert(f.passes(low_qual) == "q10")
  assert(f.passes(ok_qual) == "PASS")

  MD_QUAL = '##FILTER=<ID=c(5|100),Description="Coverage not within desired margins min:5 max:100">'
  bad_cov="20     17330   .         T      A       123    q10    NS=3;DP=1;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3"
  ok_cov="20     17330   .         T      A       123    q10    NS=3;DP=21;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3"

  f = FilterCov(5, 100)
  assert(f.md_string() == MD_QUAL)
  assert(f.passes(ok_cov) == "PASS")
  assert(f.passes(bad_cov) == "c(5|100)")

  MD_QUAL = '##FILTER=<ID=mns2,Description="Minimum number of samples with event: 2">'
  bad="20     17330   .         T      A       123    q10    NS=3;DP=1;AF=0.017;SF=12              GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3"
  ok="20     17330   .         T      A       123    q10    NS=3;DP=21;AF=0.017;SF=21,34           GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3"

  f = FilterNumSamples(2)
  assert(f.md_string() == MD_QUAL)
  assert(f.passes(ok) == "PASS")
  assert(f.passes(bad) == "mns2")

  MD_QUAL = '##FILTER=<ID=rcov4,Description="Minimum raw coverage per sample: 4">'
  bad="20     17330   .         T      A       123    q10    NS=3;DP=1;AF=0.017;SF=12;RDP=43,1,34,4              GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3"
  ok="20     17330   .         T      A       123    q10    NS=3;DP=21;AF=0.017;SF=21,34;RDP=23,10,12           GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3"
  ok2="20     17330   .         T      A       123    q10    NS=3;DP=21;AF=0.017;SF=21,34;RDP=23           GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3"
  noRDP="20     17330   .         T      A       123    q10    NS=3;DP=21;AF=0.017;SF=21,34           GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3"

  f = FilterRawCoverage(4)
  assert(f.md_string() == MD_QUAL)
  assert(f.passes(ok) == "PASS")
  assert(f.passes(ok2) == "PASS")
  assert(f.passes(bad) == "rcov4")
  assert(f.passes(noRDP) == "rcov4")

class FilterEngine():
  MIN_QUAL      = 100
  MAX_COV       = 200
  MIN_COV       = 10
  MIN_N_SAMPLES = 2
  MIN_RAW_COV   = 4
  INDEL_RE      = re.compile(r"INDEL")
  FILTER_RE     = re.compile(r"FILTER")
  filters       = []
  metadata      = []
  cline         = "" # current_line

  def __init__(self):
    self.set_defaults()
    self.load_filters()

  def run(self):
    first_snp = True
    for self.cline in sys.stdin:
      self.cline = self.cline.rstrip('\n')
      if self.cline[0] == '#' or self.INDEL_RE.search(self.cline):
        self.metadata.append(self.cline)
      else:
        if first_snp:
          self.process_metadata()
          first_snp = False
        self.process_snp()

  # Print the metadata adding our filter entries if necessary
  def process_metadata(self):
    ftp = set() # filters to print
    for l in self.metadata[:-1]: # '^##' lines
      if self.FILTER_RE.search(l):
        ftp.add(l)
      else:
        print l

    for f in self.filters: ftp.add(f.md_string())
    for f in ftp: print f
    print self.metadata[-1] # the header

  def process_snp(self):
    split    = self.cline.split()
    not_pass = set()    # filter that don't pass
    fcv      = split[6] # the filter column value

    if fcv != '.' and fcv != 'PASS': # save the current filter values
      for f in ncf.split(';'):
        not_pass.add(f)

    # Now check our filters
    for f in self.filters:
      col_val = f.passes(self.cline)
      if col_val != 'PASS':
        not_pass.add(col_val)

    if len(not_pass) == 0:
      not_pass.add('PASS')

    # print the line adding the new FILTER field
    print '\t'.join(split[0:6]) + '\t' + ';'.join(not_pass) + '\t' + '\t'.join(split[7:-1])

  def set_defaults(self):
    if len(sys.argv) == 2 and re.search(r'wgs|wes', sys.argv[1]):
      if sys.argv[1] == "wes":
        self.MAX_COV = 200
        self.MIN_COV = 10
      else:
        self.MAX_COV = 100
        self.MIN_COV = 5
    else:
      sys.stderr.write("Usage: " + sys.argv[0] + " wgs or wes" + "\n")
      sys.exit(1)

  def load_filters(self):
    self.filters.append(FilterQual(self.MIN_QUAL))
    self.filters.append(FilterCov(self.MIN_COV, self.MAX_COV))
    self.filters.append(FilterNumSamples(self.MIN_N_SAMPLES))
    self.filters.append(FilterRawCoverage(self.MIN_RAW_COV))

def main():
  #test()
  FilterEngine().run()

if __name__=='__main__':
  sys.exit(main())
