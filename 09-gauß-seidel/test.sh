#!/bin/bash

BINARY="./partdiff-par"

diff -s <($BINARY 1 1 0 2 2 85 | grep -A9 Matrix) <(cat referenz/GaussSeidel.f2 | grep -A9 Matrix)
diff -s <($BINARY 1 1 0 2 1 1e-4 | grep -A9 Matrix) <(cat referenz/GaussSeidel.f2 | grep -A9 Matrix)
diff -s <($BINARY 1 1 0 1 2 82 | grep -A9 Matrix) <(cat referenz/GaussSeidel.f1 | grep -A9 Matrix)
diff -s <($BINARY 1 1 0 1 1 1e-4 | grep -A9 Matrix) <(cat referenz/GaussSeidel.f1 | grep -A9 Matrix)
