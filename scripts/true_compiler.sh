#!/usr/bin/env bash
#

# If the argument is mpicc or mpif90, prints the true compiler that is
# being called. For the moment it only works if the -show argument is
# accepted by the wrapper (mpich, openmpi and derivatives do).
# It also works for the Cray compiler wrappers cc and ftn.

# MPI
if echo $1 | grep mpi > /dev/null; then
  if $1 -show &> /dev/null; then
      printf "("`$1 -show | cut -f 1 -d" "`")"
  else
      printf "(unknown)"
  fi
# Cray compiler wrappers
elif [ x$1 == xcc -o x$1 == xftn -o x$1 == xCC ]; then
  if $1 -show &> /dev/null; then
      printf "("`$1 -show | grep DRIVERNAME | cut -f 2 -d"="`")"
  else
      printf "(unknown)"
  fi
fi
