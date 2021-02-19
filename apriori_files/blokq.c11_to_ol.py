#!/usr/bin/python3
import sys

#Input file updated daily with lots of data including OL data
blokq = sys.argv[1]
#Output file which will contain only OL data
ol = sys.argv[2]

ol_f = open(ol,"w")

ol = False
for line in open(blokq):
    if "OCEAN LOADING CATALOG goes here:" in line:
        ol = True
    if ol:
        ol_f.write(line)
    if ol and ("END TABLE" in line):
        break


