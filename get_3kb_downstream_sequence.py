# !/usr/bin/python

from __future__ import print_function
import requests
import sys

animal = sys.argv[1]


if "horse" in animal:  
    server = "http://rest.ensembl.org/sequence/region/horse/"
    regions = ["9:62298597..62301597:1", "15:88547817..88550817:-1", "8:40892966..40895966:-1", "25:36814500..36817500:-1", "25:36804532..36807532:1", "2:28416935..28419935:1", "19:25141928..25144928:-1", "2:40443884..40446884:1", "1:88935579..88938579:1", "1:43242499..43245499:-1", "1:88924931..88927931:1", "1:88944864..88947864:-1"]

elif "cow" in animal:
    server = "http://rest.ensembl.org/sequence/region/cow/"
    regions = ["1:80647824..80650824:1",
    "8:112867044..112870044:-1",
    "11:106831644..106834644:1",
    "14:47260661..47263661:-1",
    "16:43477400..43480400:1",
    "24:35629346..35632346:-1",
    "26:6348278..6351278:1",
    "28:35602674..35605674:1",
    "28:35700096..35703096:1",
    "28:35722932..35725932:1",
    "28:35824384..35827384:1",
    "28:35837918..35840918:-1",
    "28:35848716..35851716:-1"
    ]

else:
    print("Try horse or cow")

for region in regions:
	r = requests.get(server+region+"?", headers={ "Content-Type" : "text/x-fasta"})
	print(r.text)
