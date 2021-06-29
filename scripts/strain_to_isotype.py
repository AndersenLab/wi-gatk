#!/usr/bin/env python
# -*- coding: utf-8 -*-
import requests
import sys

wi = requests.get("https://docs.google.com/spreadsheets/d/1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI/pub?output=tsv")
wi = [x.split("\t")[0:3] for x in wi.text.encode('utf-8').splitlines()]
wi = {x[0]:x[2] for x in wi[1:] if x[1] != ""}

for line in sys.stdin:
    l = line.split("\t")
    if l[0] in wi.keys():
        l[0] = wi[l[0]]
        print('\t'.join(l).strip())