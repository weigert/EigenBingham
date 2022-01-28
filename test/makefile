################################################################################
# dirstat
# Author: Nicholas McDonald
# Version: 1.0
################################################################################

# Compilation Configuration
CC = g++ -std=c++17
CFLAGS = -Wfatal-errors -O3
LFLAGS = -I/usr/local/include -L/usr/local/lib

# Make All
all: dirstat
.PHONY: all

dirstat: main.cpp
	$(CC) main.cpp $(CFLAGS) $(LFLAGS) -o main
