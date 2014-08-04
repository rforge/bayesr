REPOS=http://svn.gwdg.de/svn/bayesx/trunk
USER=guest
PASSWD=guest
DIRS="adaptiv alex andrea bib dag graph leyre mcmc psplines samson structadd"
FILES="export_type.h main.cpp values.h"
mkdir -p src/bayesxsrc
cd src/bayesxsrc
for i in $DIRS ; do
  svn checkout -r1333 --username "${USER}" --password "${PASSWD}" $REPOS/$i $i
done
for i in $FILES ; do
  svn export --username "${USER}" --password "${PASSWD}" $REPOS/$i $i
done
svn update --revision r1321 --username "${USER}" --password "${PASSWD}" structadd/superbayesreg.cpp
svn update --revision r1321 --username "${USER}" --password "${PASSWD}" structadd/superbayesreg.h
svn update --revision r1336 --username "${USER}" --password "${PASSWD}" mcmc/fullcond_merror.cpp
cd ..
cp dev-Makefile Makefile
cp dev-Makefile.win Makefile.win
