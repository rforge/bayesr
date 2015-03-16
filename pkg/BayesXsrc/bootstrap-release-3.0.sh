REPOS=http://svn.gwdg.de/svn/bayesx/trunk
USER=guest
PASSWD=guest
DIRS="adaptiv alex andrea bib dag graph leyre mcmc psplines samson structadd"
FILES="export_type.h main.cpp values.h"
mkdir -p src/bayesxsrc
cd src/bayesxsrc
for i in $DIRS ; do
  svn checkout -r1474 --username "${USER}" --password "${PASSWD}" $REPOS/$i $i
done
for i in $FILES ; do
  svn export --username "${USER}" --password "${PASSWD}" $REPOS/$i $i
done
cd ..
cp dev-Makefile Makefile
cp dev-Makefile.win Makefile.win
