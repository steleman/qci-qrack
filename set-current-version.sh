#!/bin/bash

export HERE="`pwd`"
export QRACK="`basename ${HERE}`"

if [ "${QRACK}" != "qrack" ] ; then
  echo "You are in the wrong directory."
  exit 1
fi

git checkout vm6502qv9.5.2

