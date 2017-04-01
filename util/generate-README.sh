#!/bin/bash

## this files is used to generate a README.md file from the jupyter notebook.
## it converts the notebook to markdown, and outputs it to README.md .

jupyter nbconvert --to markdown ../HIP-fitting-usage.ipynb 
rm -rf HIP-fitting-usage_files
mv ../HIP-fitting-usage_files ./

## replace pngs path
cat ../HIP-fitting-usage.md |  sed -e 's/HIP-fitting-usage_files/util\/HIP-fitting-usage_files/g' > ../README.md
rm ../HIP-fitting-usage.md
