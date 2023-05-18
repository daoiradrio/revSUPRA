#!/bin/bash
gawk '/TOTAL ENERGY/ { val=$4; printf("%18.12f", val); }' $1 >> $2