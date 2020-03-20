#!/bin/bash

#tmux attach -t 0

# desktop: 5/6 cores
for i in $(seq 5);
	do { /usr/bin/time -o time-$i ./KLP ; } > result-$i &
done
