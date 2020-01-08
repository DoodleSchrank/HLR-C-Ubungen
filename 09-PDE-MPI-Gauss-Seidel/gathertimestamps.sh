#!/bin/sh

for f in output/*
do
	cat $f | grep Berechnungszeit >> timestamps
done
