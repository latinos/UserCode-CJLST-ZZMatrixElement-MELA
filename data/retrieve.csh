#!/bin/tcsh -f
# Script to retrieve the libmcfm.so library from the link specified in download.url

cd `dirname $0`

if (! -e libmcfm.so) then
  cat download.url | xargs wget -q
endif

