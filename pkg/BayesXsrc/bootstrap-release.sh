REPOS=http://svn.gwdg.de/svn/bayesx/trunk
USER=guest
PASSWD=guest
DIRS="adaptiv alex andrea bib dag graph leyre mcmc psplines samson structadd"
FILES="export_type.h main.cpp values.h"
mkdir -p src/bayesxsrc
cd src/bayesxsrc
for i in $DIRS ; do
  svn checkout --revision 953 --username "${USER}" --password "${PASSWD}" $REPOS/$i $i
done
for i in $FILES ; do
  svn export --revision 953 --username "${USER}" --password "${PASSWD}" $REPOS/$i $i
done
svn update --revision 1233 --username "${USER}" --password "${PASSWD}" bib/Random.cpp
cd ..
cp rel-Makefile Makefile
cp rel-Makefile.win Makefile.win

