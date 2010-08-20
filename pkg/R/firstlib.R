.onLoad <- function(...) 
	{ 
    	check <- suppressPackageStartupMessages(require(fields,quietly=TRUE))
    	if(!check)
		{
		cat("BayesR Info: Installing required package fields in default library location!\n")
    		download.packages(pkgs="fields",destdir=.Library)
		install.packages("fields")
		require("fields")
		}
    	check <- suppressPackageStartupMessages(require(geometry,quietly=TRUE))
    	if(!check)
		{
		cat("BayesR Info: Installing required package geometry in default library location!\n")
    		download.packages(pkgs="geometry",destdir=.Library)
		install.packages("geometry")
		require("geometry")
		}
    	check <- suppressPackageStartupMessages(require(msm,quietly=TRUE))
    	if(!check)
		{
		cat("BayesR Info: Installing required package msm in default library location!\n")
    		download.packages(pkgs="msm",destdir=.Library)
		install.packages("msm")
		require("msm")
		}
    	check <- suppressPackageStartupMessages(require(mgcv,quietly=TRUE))
    	if(!check)
		{
		cat("BayesR Info: Installing required package mgcv in default library location!\n")
    		download.packages(pkgs="mgcv",destdir=.Library)
		install.packages("mgcv")
		require("mgcv")
		}
    	}

print.bayesr.version <- function()
    	{ 
    	library(help=BayesR)$info[[1]] -> version
    	version <- version[pmatch("Version",version)]
    	um <- strsplit(version," ")[[1]]
    	version <- um[nchar(um)>0][2]
    	cat(paste("This is BayesR",version,"\n"))
    	}

.onAttach <- function(...) 
	{ 
    	# options(help_type="html") 
    	print.bayesr.version()
    	}

.onUnload <- function(libpath) library.dynam.unload("BayesR", libpath)

.First.lib <- function(lib, pkg) 
    	{
    	library.dynam("BayesR", pkg, lib)
    	print.bayesr.version()
   	}
