## Before checking do the following:
## 1. cd in the BayesXsrc directory.
## 2. remove all files with
##      rm -rf *
## 3. svn up
## 4. Obtain latest BayesX sources with
##      sh ./bootstrap-devel.sh
## 5. Build the .tar.gz with
##      R CMD build BayesXsrc

library("rhub")

## Validate your e-mail address.
validate_email()

## See all platforms that can be used for checking.
platforms()

## Run check with Rdevel.
check_with_rdevel("BayesXsrc_3.0-0.tar.gz")

## Run checks with address sanitizers.
check_with_sanitizers("BayesXsrc_3.0-0.tar.gz")
