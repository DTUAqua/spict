#!/bin/bash
echo "library(shiny); runApp('spictapp', port=5000, display.mode='normal')" | R --slave

