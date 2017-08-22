#!/usr/bin/python

from __future__ import print_function
import argparse

!parser = argparse.ArgumentParser(description = "Description of the program to be written", epilog = "Name, date")

parser.add_argument("input", help = "A file input, specification not required")
parser.add_argument("-ref", dest = 'ref_file', help = "Requires the specification noted.")

args = parser.parse_args()
