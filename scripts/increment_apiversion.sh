#!/bin/bash

cmakefile=`dirname $0`/../CMakeLists.txt
dummyconfig=`dirname $0`/../src/papilo/Config.hpp

lineapiversion=$(sed -n '/PAPILO_API_VERSION/p' $cmakefile)

version=$(echo "$lineapiversion" | grep -o -E '[0-9]+')

echo "Raising to API version $((version + 1))!"

sed -i "s/^set(PAPILO_API_VERSION .*/set(PAPILO_API_VERSION $((version + 1)))/g" $cmakefile

sed -i "s/^#define PAPILO_API_VERSION .*/#define PAPILO_API_VERSION $((version + 1))/g" $dummyconfig
