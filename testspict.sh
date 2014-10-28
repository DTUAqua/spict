#!/bin/bash
./package.sh
echo "library(spict); example(spict)" | R --slave
