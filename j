#!/bin/bash

# This program compiles and runs an OSP project

# Installation: Copy to your home directory and change the marked line below
# Usage: ~/osp ch99/MyJavaApp
# Note: Do not include the ".java" extension on your program's name

# The following line should point to your osp_project directory:
# cd ~/home/mnt756/osp_project
java -classpath classes/ org/opensourcephysics/sip/$1
