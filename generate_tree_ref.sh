#!/bin/bash

#File: tree-md

tree=$(tree -d -tf --noreport -I '*~' --charset ascii $1 |
       sed -e 's/| \+/  /g' -e 's/[|`]-\+/ */g' -e 's:\(* \)\(\(.*/\)\([^/]\+\)\):\1[\4](\2):g' -e 's/\.\//https:\/\/github.com\/davidtrestini\/PNpedia\/tree\/main\//g'  -e ':loop' -e 's/(\(.*\) \(.*\))/(\1%20\2)/g' -e 't loop')
echo "${tree}" > "./generated_tree.md"
