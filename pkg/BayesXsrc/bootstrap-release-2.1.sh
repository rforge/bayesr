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
svn update --revision 1233 --username "${USER}" --password "${PASSWD}" bib/Random.h
svn update --revision 1246 --username "${USER}" --password "${PASSWD}" bib/remlreg.cpp
svn update --revision 1246 --username "${USER}" --password "${PASSWD}" bib/bayesreg.cpp
svn update --revision 1246 --username "${USER}" --password "${PASSWD}" bib/stepwisereg.cpp
svn update --revision 1246 --username "${USER}" --password "${PASSWD}" bib/mapobject.cpp
svn update --revision 1246 --username "${USER}" --password "${PASSWD}" bib/dataobj.cpp
svn update --revision 1247 --username "${USER}" --password "${PASSWD}" mcmc/fullcond.h
svn update --revision 1247 --username "${USER}" --password "${PASSWD}" mcmc/mcmc_const.h
svn update --revision 1247 --username "${USER}" --password "${PASSWD}" mcmc/distribution.h
svn update --revision 1249 --username "${USER}" --password "${PASSWD}" bib/clstring.h
svn update --revision 1249 --username "${USER}" --password "${PASSWD}" leyre/nbinomial.cpp
cd ..
cp rel-Makefile Makefile
cp rel-Makefile.win Makefile.win
