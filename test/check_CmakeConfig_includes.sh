#!/bin/bash
#Project DIR
path=$1
expression=$2

# shellcheck disable=SC2086
FILES_WITH_EXPRESSION="$(grep -ilR $expression "$path"/src)"
FILES_WITH_IMPORT_CONFIG="$(grep -ilR "papilo/Config.hpp" "$path"/src)"

for FILE_WITHEXPRESSION in ${FILES_WITH_EXPRESSION[@]}; do
  echo "$FILE_WITHEXPRESSION"
  found=false
  if [[ "$FILE_WITHEXPRESSION" = *in ]]; then
      continue
  fi
  if [[ "$FILE_WITHEXPRESSION" = *Config.hpp ]]; then
      continue
  fi
  for i in ${FILES_WITH_IMPORT_CONFIG[@]}; do
    if [ "$i" = "$FILE_WITHEXPRESSION" ]; then
        found=true
        break
    fi
  done
  if [ $found = true ]; then
    echo "found"
  else
    echo "\"#import papilo/Config.hpp\" is missing in the last printed file!"
    exit 1
  fi
done

