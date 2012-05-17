Installation:
=============

Build & Install
---------------
cd BayesXsrc
sh ./bootstrap.sh
cd ..
R CMD build BayesXsrc
R CMD INSTALL BayesXsrc

Build Windows Multiarch
-----------------------
R CMD build BayesXsrc
R CMD INSTALL --merge-multiarch --build BayesXsrc_2.1-0.tar.gz

Build OS X Universal
--------------------
rm -rf /tmp/pkg
mkdir /tmp/pkg
mkdir -p builds/osx/universal
R --arch=i386 CMD INSTALL -c -l /tmp/pkg BayesXsrc
R --arch=x86_64 CMD INSTALL -c -l /tmp/pkg --libs-only BayesXsrc
tar fcvz builds/osx/universal/BayesXsrc_2.1-0.tgz -C /tmp/pkg BayesXsrc
rm -rf /tmp/pkg


Build status:
=============
Success:
- R 2.14.0 / gcc-4.2.1 / Mac OS X 10.6 / x86_64
- R 2.14.0 / gcc-4.4.5 / Debian 6      / x86_64

Issues:
- R 2.12.0 / gcc-4.5   / FreeBSD 8.2,  / i386    (Permission Problems, no install)

- Rtools 2.14, gcc-4.4.5 mingw does not like #include without following a white-space.

Debugging:
- R 2.14.0, gcc-4.5,   Windows XP,    i386


History:
========

BayesX SVN revision | BayesX version | BayesXsrc version | Comment
--------------------+----------------+-------------------+--------
                948 |            2.1 |             2.1-0 | First GPL-2 release of BayesX
                953 |              - |             2.1-1 | Fixes to comply with gcc-4.7		
		    |                |                   |
