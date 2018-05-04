#!/bin/sh


qsub -cwd -e io.trash/ -o io.trash/ -t 1-300:1 -tc 20 ./call_univ.sh

