R=http://svn.gwdg.de/svn/bayesx/trunk
U=guest
P=guest
DIRS="adaptiv alex andrea bib dag graph leyre mcmc psplines samson structadd"
mkdir -p src/bayesxsrc
cd src/bayesxsrc
for i in $DIRS ; do
  svn co --username "${U}" --password "${P}" $R/$i $i
done
FILES="export_type.h main.cpp values.h"
for i in $FILES ; do
  svn export --username "${U}" --password "${P}" $R/$i $i
done

