#! /usr/bin/env python3

import gzip

input_maf_file_name = "hprc-v1.1-mc-grch38-full_fix.maf.gz"
output_maf_file_prefix = "hprc-v1.1-mc-grch38-"
output_maf = False

seen_chromosomes = {}
in_block = False

with gzip.open(input_maf_file_name) as input_file:
  for line in input_file:
    line = line.decode('UTF-8')
    if line and line[0] == 's':
      if not in_block:
        chromosome = line.split("\t")[1].split(".")[1]
        if chromosome in seen_chromosomes:
          output_maf.write("\na\n")
        else:
          seen_chromosomes[chromosome] = True
          print(f"New chromosome: {chromosome}")
          if output_maf:
            output_maf.close()
          output_maf = open(f"{output_maf_file_prefix}{chromosome}.maf", 'w')
          output_maf.write("##maf version=1 scoring=N/A\n\na\n")
        in_block = True
      output_maf.write(line)
    else:
      in_block = False
