#!/usr/bin/python

'''
This script will count the number of Effector, CAZY, combined
effector genes and secreted genes in a coexpressed module compared to the rest
of the genome
'''

import argparse
import os

# -----------------------------------------------------
# Step 1
# Import command line arguments
# -----------------------------------------------------

ap = argparse.ArgumentParser()
ap.add_argument('--Module_Effector', required=True, type=str, help='List of Effectors \
in a coexpression module')
ap.add_argument('--Module_CAZY', required=True, type=str, help='List of CAZYs \
in a coexpression module')
ap.add_argument('--Module_Effectors_Combined', required=True, type=str, help='List of Combined Effectors \
in a coexpression module')
ap.add_argument('--Module_Secreted', required=True, type=str, help='List of \
Secreted genes in a coexpression module')
ap.add_argument('--Module_Genes', required=True, type=str, help='Complete list \
of genes in a coexpression module')
ap.add_argument('--Module_Name', required=True, type=str, help='Name of the \
module being analysed')
ap.add_argument('--Genome_Effectors', required=True, type=float, help='Number of \
total Effectors in the genome')
ap.add_argument('--Genome_CAZY', required=True, type=float, help='Number of \
total CAZY in the genome')
ap.add_argument('--Genome_Effectors_Combined', required=True, type=float, help='Number of \
total Effectors combined in the genome')
ap.add_argument('--Genome_Secreted', required=True, type=float, help='Number of \
total Secreted proteins in the genome')
ap.add_argument('--Genome_Genes', required=True, type=float, help='Number of \
total genes in the genome')
ap.add_argument('--OutDir', required=True, type=str, help='Directory to write \
output files to')
conf = ap.parse_args()

# -----------------------------------------------------
# Step 2
# Load inputs into python variables
# -----------------------------------------------------

Effectors = []
with open(conf.Module_Effector) as f:
    for line in f.readlines():
        Effectors.append(line)

Effectors_Set = set(Effectors)

CAZY = []
with open(conf.Module_CAZY) as f:
    for line in f.readlines():
        CAZY.append(line)

CAZY_Set = set(CAZY)

Effectors_Combined = []
with open(conf.Module_Effectors_Combined) as f:
    for line in f.readlines():
        Effectors.append(line)

Effector_Combined_Set = set(Effectors_Combined)

Secreted = []
with open(conf.Module_Secreted) as f:
    for line in f.readlines():
        Secreted.append(line)

Secreted_Set = set(Secreted)

Genes = []
with open(conf.Module_Genes) as f:
    for line in f.readlines():
        Genes.append(line)

Gene_Set = set(Genes)

Module_Name = conf.Module_Name
Genome_Effectors = conf.Genome_Effectors
Genome_CAZY = conf.Genome_CAZY
Genome_Effectors_Combined = conf.Genome_Effectors_Combined
Genome_Secreted = conf.Genome_Secreted
Genome_Genes = conf.Genome_Genes
OutDir = conf.OutDir
cwd = os.getcwd()

# -----------------------------------------------------
# Step 3
# Obtain counts for all gene types in module and rest of genome
# -----------------------------------------------------

# Modules

Mod_Effector_Num = len(Effector_Set)
Mod_CAZY_Num = len(CAZY_Set)
Mod_Effector_Combined_Num = len(Effector_Combined_Set)
Mod_Secreted_Num = len(Secreted_Set)
Mod_Gene_Num = len(Gene_Set)

# Genome

Gen_Effector_Num = Genome_Effectors - Mod_Effector_Num
Gen_CAZY_Num = Genome_CAZY - Mod_CAZY_Num
Gen_Effector_Combined_Num = Genome_Effectors_Combined - Mod_Effector_Combined_Num
Gen_Secreted_Num = Genome_Secreted - Mod_Secreted_Num
Gen_Gene_Num = Genome_Genes - Mod_Gene_Num

# -----------------------------------------------------
# Step 4
# Write output files
# -----------------------------------------------------

Effector_File = "_".join([Module_Name, "Effector_Fishertable.txt"])
Effector_Out = "/".join([cwd, OutDir, Effector_File])

with open(Effector_Out, "w") as o:
    Line1 = "\t".join(["Effector", str(Mod_Effector_Num), str(Gen_Effector_Num)]) + "\n"
    Mod_Genes = Mod_Gene_Num - Mod_Effector_Num
    Gen_Genes = Gen_Gene_Num - Mod_Gene_Num - Gen_Effector_Num
    Line2 = "\t".join(["Other Genes", str(Mod_Genes), str(Gen_Genes)]) + "\n"
    o.write("".join([Line1, Line2]))
    o.close()

CAZY_File = "_".join([Module_Name, "CAZY_Fishertable.txt"])
CAZY_Out = "/".join([cwd, OutDir, CAZY_File])

with open(CAZY_Out, "w") as o:
    Line1 = "\t".join(["CAZY", str(Mod_CAZY_Num), str(Gen_CAZY_Num)]) + "\n"
    Mod_Genes = Mod_Gene_Num - Mod_CAZY_Num
    Gen_Genes = Gen_Gene_Num - Mod_Gene_Num - Gen_CRN_Num
    Line2 = "\t".join(["Other Genes", str(Mod_Genes), str(Gen_Genes)]) + "\n"
    o.write("".join([Line1, Line2]))
    o.close()

Effector_Combined_File = "_".join([Module_Name, "Effector_Combined_Fishertable.txt"])
Effector_Combined_Out = "/".join([cwd, OutDir, Effector_Combined_File])

with open(Effector_Combined_Out, "w") as o:
    Line1 = "\t".join(["Effector_Combined", str(Mod_Effector_Num),
                      str(Gen_Effector_Num)]) + "\n"
    Mod_Genes = Mod_Gene_Num - Mod_Effector_Combined_Num
    Gen_Genes = Gen_Gene_Num - Mod_Gene_Num - Gen_Effector_Combined_Num
    Line2 = "\t".join(["Other Genes", str(Mod_Genes), str(Gen_Genes)]) + "\n"
    o.write("".join([Line1, Line2]))
    o.close()

Secreted_File = "_".join([Module_Name, "Secreted_Fishertable.txt"])
Secreted_Out = "/".join([cwd, OutDir, Secreted_File])

with open(Secreted_Out, "w") as o:
    Line1 = "\t".join(["Secreted", str(Mod_Secreted_Num),
                       str(Gen_Secreted_Num)]) + "\n"
    Mod_Genes = Mod_Gene_Num - Mod_Secreted_Num
    Gen_Genes = Gen_Gene_Num - Mod_Gene_Num - Gen_Secreted_Num
    Line2 = "\t".join(["Other Genes", str(Mod_Genes), str(Gen_Genes)]) + "\n"
    o.write("".join([Line1, Line2]))
    o.close()
