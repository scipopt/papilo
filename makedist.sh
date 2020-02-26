#!/bin/bash -e

VERSION=$(grep project.papilo CMakeLists.txt | cut -d ' ' -f 3 | cut -d ')' -f 1)
NAME="papilo-$VERSION"
rm -f $NAME.tgz
rm -f $NAME.tar

echo ">>> Packaging $NAME."

# echo "store git hash"
GITHASH=`git describe --always --dirty  | sed 's/^.*-g//'`
sed -i "s/undef PAPILO_GITHASH_AVAILABLE/define PAPILO_GITHASH_AVAILABLE/g" src/papilo/Config.hpp
sed -i "s/undef PAPILO_GITHASH/define PAPILO_GITHASH \"$GITHASH\"/g" src/papilo/Config.hpp

# Before we create a tarball change the directory and file rights in a command way
echo "adjust file modes"
git ls-files | xargs dirname | sort -u | xargs chmod 750
git ls-files | xargs chmod 640
git ls-files "*.sh" "*.py" "scripts/*" | grep -v external | xargs chmod 750

# pack files tracked by git and append $NAME to the front
git ls-files -c | xargs tar --transform "s|^|${NAME}/|" -cvhf $NAME.tar \
--exclude="*~" \
--exclude=".*"

# compress the archive
gzip -c $NAME.tar > $NAME.tgz

# remove temporary archive
rm -f $NAME.tar

echo ""
echo "check version numbers ($VERSION):"
grep -H "VERSION" src/papilo/Config.hpp
grep -H project.papilo CMakeLists.txt
# tail src/presol/git_hash.cpp
