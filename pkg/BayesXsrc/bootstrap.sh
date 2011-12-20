CVS_RSH=ssh
CVSROOT=:ext:lang@penelope.stat.uni-muenchen.de:/home/lang/cvs
export CVS_RSH
export CVSROOT
FILES="adaptiv alex andrea bib dag export_type.h graph leyre main.cpp mcmc psplines samson structadd values.h"
for i in $FILES ; do
  MODULES="${MODULES} bayesx/$i"
done
cd src
cvs checkout ${MODULES}
mv bayesx bayesxsrc

