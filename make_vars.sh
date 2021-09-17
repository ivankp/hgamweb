#!/bin/sh

./mxaods.py mxaods_h024.json \
| sed 's:^/eos/atlas/atlascerngroupdisk/phys-higgs/HSG1/MxAOD/:/msu/data/t3work9/ivanp/MxAODs/:' \
| xargs ./bin/make_vars h024
