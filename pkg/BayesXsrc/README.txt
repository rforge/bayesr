installation:
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
R CMD INSTALL --merge-multiarch --build bayesxsrc_0.1.tar.gz

Build OS X Universal
--------------------
rm -rf /tmp/pkg
mkdir /tmp/pkg
mkdir -p builds/osx/universal
R --arch=i386 CMD INSTALL -c -l /tmp/pkg bayesxsrc
R --arch=x86_64 CMD INSTALL -c -l /tmp/pkg --libs-only bayesxsrc
tar fcvz builds/osx/universal/bayesxsrc_0.1.tgz -C /tmp/pkg bayesxsrc
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

