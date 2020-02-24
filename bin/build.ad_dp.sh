#!/bin/bash
# This script can be used to build the ad_dp filter

nim c -d:release -o:ad_dp_macos ad_dp.nim
docker run --rm -v `pwd`:/usr/src/app -w /usr/src/app nimlang/nim nim c -d:release -o:ad_dp_linux ad_dp.nim
