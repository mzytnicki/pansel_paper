#! /usr/bin/env python3

import glob, gzip

input_maf_file_names = ["hprc-v1.1-mc-grch38-chr1.maf.gz", "hprc-v1.1-mc-grch38-full.maf.gz"]
output_maf_file_names = ["hprc-v1.1-mc-grch38-chr1_fix.maf", "hprc-v1.1-mc-grch38-full_fix.maf"]
input_dnd_file_name = "mashtree.dnd"
output_dnd_file_name = "mashtree_fix.dnd"
delims = "():;"

for input_maf_file_name, output_maf_file_name in zip(input_maf_file_names, output_maf_file_names):
  with gzip.open(input_maf_file_name, "r") as input_file:
    with open(output_maf_file_name, "w") as output_file:
      for line in input_file:
        line = line.decode('UTF-8')
        if line:
          if line[0] == 's':
            line = line.strip().split()
            if not (line[1].startswith("GRCh38") or line[1].startswith("CHM13")):
              # HG03516.1.JAGYYT010000183.1 -> HG03516_1.JAGYYT010000183.1
              line[1] = line[1].replace(".", "_", 1)
            output_file.write("%s\n" % "\t".join(line))
          else:
            output_file.write(line)
        else:
          output_file.write("\n")

with open(input_dnd_file_name) as input_file:
  with open(output_dnd_file_name, 'w') as output_file:
    for line in input_file:
      in_name = False
      current_name = ""
      for c in line:
        if c in delims:
          if in_name:
            current_name = current_name.replace("#", "_")[10:]
            if current_name == "grch38":
                current_name = "GRCh38"
            elif current_name == "chm13":
                current_name = "CHM13_1_1"
            output_file.write(current_name)
            in_name = False
            current_name = ""
          output_file.write(c)
        elif in_name:
          current_name += c
        elif c == "h":
          in_name = True
          current_name = "h"
        else:
          output_file.write(c)
