#!/bin/sh


qsub -cwd -e io.trash/ -o io.trash/ -t 1-450:1 -tc 20 ./call_sim.sh
