#!/bin/bash

#File: tree-md

tree=$(tree -d -tf --noreport -I '*~' --charset ascii $1 |
       sed    -e 's:\.\/:\/:g' -e 's:\/.*\/:\/:g' -e 's:\/::g')
echo "${tree}" > "./generated_tree.md"
