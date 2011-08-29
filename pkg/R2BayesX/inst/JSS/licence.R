library("tools")


## Suchen und setzen von Code-Files, z.B. mit Dateiendungen ('exts')
## .cpp und .h. Der Lizenz 'text' wird dann automatisch an den Anfang jeder
## betreffenden Datei gestellt.
## dir: die directory des Projekts.
## text: der Lizenz Vermerk, der am Anfang der jeweiligen Datei eingefuegt werden soll.
## exts: nur in Dateien mit den in exts spezifizierten Endungen werden modifiziert.
licence <- function(dir, text, exts = c("cpp", "h"))
{
  dir <- path.expand(dir)
  files <- list.files(dir)
  isdir <- file.info(file.path(dir, files))[["isdir"]]
  nodir <- list_files_with_exts(dir = dir, exts = exts, full.names = FALSE)
  if(length(nodir)) {
    for(i in nodir) {
      con <- file.path(dir, i)
      Sys.chmod(con)
      file <- readLines(con = con)
      file <- c(text, file)
      writeLines(text = file, con = con)
    }
  }
  if(length(files[isdir])) {
    for(i in files[isdir]) {
      licence(dir = file.path(dir, i), text = text, exts = exts)
    }
  }
}


## Der Vermerk auf die GPL2 Lizenz, wird an den Anfang von
## bspw. den .cpp und .h Dateien gestellt.
licence.text <- "/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2011  Christiane Belitz, Andreas Brezger,
Thomas Kneib (Georg-August-Universitaet Goettingen),
Stefan Lang (Leopold-Franzens-Universitaet Innsbruck)
stefan.lang@uibk.ac.at

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */\n\n"


## Kleines Beispiel mit dem entpackten sourcecode von der BayesX Homepage.
licence("~/bayesx", text = licence.text)
