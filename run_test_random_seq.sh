#!/bin/bash

# Path to EasyQuant repository/master folder
EASYQUANT_PATH=$1

rm -rf $EASYQUANT_PATH/test_out/

python $EASYQUANT_PATH/random_sequence_generator.py \
    -c 1000 \
    -p 0,200,800,1000 \
    -s 10 \
    -r 50 \
    -n 100 \
    -i 40 \
    -t $EASYQUANT_PATH/test.tsv \
    -1 $EASYQUANT_PATH/test_R1.fastq.gz \
    -2 $EASYQUANT_PATH/test_R2.fastq.gz

python $EASYQUANT_PATH/easy_quant.py \
    -i $EASYQUANT_PATH/test_R1.fastq.gz \
    $EASYQUANT_PATH/test_R2.fastq.gz \
    -s $EASYQUANT_PATH/test.tsv \
    -d 20 \
    -o $EASYQUANT_PATH/test_out/
