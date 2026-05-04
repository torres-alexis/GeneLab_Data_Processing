#!/usr/bin/env bash
set -e

# Usage: clean_paths.sh <file>
# Removes local filesystem paths from processing-info files

file="$1"

sed -i \
    -e 's|/[A-Za-z0-9._/-]*/\(OSD-[0-9][A-Za-z0-9._-]*/RawData/[^" )]*\)|\1|g' \
    -e 's|/[A-Za-z0-9._/-]*/\(OSD-[0-9][A-Za-z0-9._-]*/results/[^" )]*\)|\1|g' \
    -e 's|/[A-Za-z0-9._/-]*/\(GLDS-[0-9][A-Za-z0-9._-]*/RawData/[^" )]*\)|\1|g' \
    -e 's|/[A-Za-z0-9._/-]*/\(GLDS-[0-9][A-Za-z0-9._-]*/results/[^" )]*\)|\1|g' \
    -e 's|/[A-Za-z0-9._/-]*/workflow_code/\([^" )]*\)|workflow_code/\1|g' \
    -e 's|/[A-Za-z0-9._/-]*/main.nf|main.nf|g' \
    -e 's|/[A-Za-z0-9._/-]*/work/[A-Za-z0-9][A-Za-z0-9._/-]*|<workdir>|g' \
    -e 's|"/home/[^"]*"|"<path-removed-for-security-purposes>"|g' \
    -e "s|'/home/[^']*'|'<path-removed-for-security-purposes>'|g" \
    -e 's|/home/[^" )]*|<path-removed-for-security-purposes>|g' \
    -e 's|"/global/[^"]*"|"<path-removed-for-security-purposes>"|g' \
    -e "s|'/global/[^']*'|'<path-removed-for-security-purposes>'|g" \
    -e 's|/global/[^" )]*|<path-removed-for-security-purposes>|g' \
    "$file"