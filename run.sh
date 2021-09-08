#!/bin/sh

./bin/hgamweb '{
  "set": "dat",
  "V": {
    "name": "pT_yy",
    "nbins": 100
  },
  "B": {
    "ndiv": 2,
    "deg": 3,
    "exp": false
  },
  "S": {
    "ndiv": 4,
    "pow": [false,false]
  },
  "edges": [
    0,5,10,15,20,25,30,35,45,60,80,100,120,140,170,200,250,350,450
  ]
}'
