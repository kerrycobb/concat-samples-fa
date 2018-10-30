let doc = """
Usage:
  concat-fasta-samples INPUT [OUTPUT | -f FORMAT]

Arguments:
  INPUT     Input file path
  OUTPUT    Output file path

Options:
  -h, --help   Show this screen
  -f FORMAT    Output file format, not needed if using OUTPUT argument,
                expects: 'nexus' or 'phylip'
"""

import strutils
import tables
import os
import docopt

# Parse command line arguments
let args = docopt(doc)
var
  input_path = $args["INPUT"]
  out_path: string
  format: string
if args["OUTPUT"]:
  out_path = $args["OUTPUT"]
  if out_path.endsWith(".phy"):
    format = "phylip"
  elif out_path.endsWith(".phylip"):
    format = "phylip"
  elif out_path.endsWith(".nex"):
    format = "nexus"
  elif out_path.endsWith(".nexus"):
    format = "nexus"
  else:
    echo "Error: Invalid output file extension. Expected: '.phy', '.phylip', '.nex', or 'nexus'"
    quit()
else:
  if args["-f"]:
    format = $args["-f"]
    if format == "phylip":
      out_path = changeFileExt(input_path, ".phy")
    elif format == "nexus":
      out_path = changeFileExt(input_path, ".nex")
    else:
      echo "Error: Invalid file format. Expected: 'nexus' or 'phylip'"
      quit()
  else:
    format = "phylip"
    out_path = changeFileExt($args["INPUT"], ".phy")

# Parse input into hash table
var
  in_handle = open(input_path)
  IUPAC = {"AC":"M", "CA":"M", "AG":"R", "GA":"R", "AT":"W", "TA":"W", "CG":"S", "GC":"S", "CT":"Y",  "TC":"Y", "GT":"K", "TG":"K"}.toTable
  samples = initTable[int, string]()
  loci = initOrderedTable[int, int]()
  sequences = initTable[int, Table[int, string]]()
while true:
  try:
    var line = readLine(in_handle)
    if line.startsWith(">"):
      var
        header = line.splitWhiteSpace()
        sample_id = header[1][1 .. ^2]
        header_split = header[0][1 .. ^10].split('_')
        locus_num = header_split[1].parseInt()
        sample_num = header_split[3].parseInt()
        sequence = readLine(in_handle)
      discard readLine(in_handle)
      var
        sequence2 = readLine(in_handle)
        seq_len = len(sequence)
        cons: seq[char] = @[]
      for i in 0 ..< seq_len:
        var
          a = sequence[i]
          b = sequence2[i]
        if a == 'N' and b == 'N':
          cons.add('N')
        elif a == b:
          cons.add(a)
        elif a == 'N' or a == '-':
          cons.add(b)
        elif b == 'N' or b == '-':
          cons.add(a)
        else:
          cons.add(IUPAC[a & b])
      var consensus = cons.join()
      if not samples.hasKey(sample_num):
        samples[sample_num] = sample_id
      if not loci.hasKey(locus_num):
        loci[locus_num] = seq_len
      if not sequences.hasKey(sample_num):
        sequences[sample_num] = {locus_num: consensus}.toTable
      else:
        sequences[sample_num][locus_num] = consensus
    else:
      discard
  except IOError:
    break
in_handle.close()

# Compute number of characters in alignment
var nchars = 0
for key, value in loci:
  nchars += value

# Output data paritions into nexus file
var out_handle = open(out_path, fmWrite)
if format == "nexus":
  out_handle.writeLine("NEXUS\n\nBEGIN SETS;")
  var cnt = 0
  for key, value in loci:
    var start = cnt + 1
    cnt = cnt + value
    out_handle.writeLine("charset " & $key & " = " & $start & "-" & $cnt & ";")
    out_handle.writeLine("END;")
    out_handle.writeLine("BEGIN DATA;\nDIMENSIONS NTAX=$1 NCHAR=$2;\nFORMAT DATATYPE=DNA MISSING=? GAP=-;\nMATRIX" % [$len(samples) , $nchars])
# Output data partitions into a partition file
elif format == "phylip":
  var
    part_path = out_path.changeFileExt("partition")
    part_handle = open(part_path, fmWrite)
    cnt = 0
  for key, value in loci:
    var start = cnt + 1
    cnt = cnt + value
    part_handle.writeLine("DNA, " & $key & " = " & $start & "-" & $cnt)
  part_handle.close()
  out_handle.writeLine($len(samples) & "\t" & $nchars)

for key, value in samples:
  var seqs: seq[string] = @[]
  for k, v in loci:
    if sequences[key].hasKey(k):
      seqs.add(sequences[key][k])
    else:
      seqs.add('N'.repeat(loci[k]))
  out_handle.writeLine(value & "\t" & seqs.join())

if format == "nexus":
  out_handle.writeLine(";\nEND;")

out_handle.close()
echo "Finished converting " & input_path
