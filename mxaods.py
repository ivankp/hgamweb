#!/usr/bin/env python3

import sys, json

if len(sys.argv)!=2:
    print('usage:',sys.argv[0],'mxaods.json')
    sys.exit(1)

def flatten(xs,d=''):
    for x in xs:
        if isinstance(x,list):
            flatten(x[1],d+x[0])
        else:
            print(d+x)

with open(sys.argv[1]) as f:
    flatten(json.load(f))
