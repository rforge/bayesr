require("tools", quietly = TRUE)


## Write 00check.log in HTML, thanks to Kurt.
write_check_log_as_HTML <- function(log, out = "", subsections = FALSE)
{
    if(out == "") 
        out <- stdout()
    else if(is.character(out)) {
        out <- file(out, "wt")
        on.exit(close(out))
    }
    if(!inherits(out, "connection")) 
        stop("'out' must be a character string or connection")
    
    lines <- readLines(log, warn = FALSE, skipNul = TRUE)[-1L]
    ## The first line says
    ##   using log directory '/var/www/R.check/......"
    ## which is really useless ...

    ## Encoding ...
    ## Get the session charset from the log (added in R 2.7.0).
    pos <- grep("^[*] using session charset:", lines, useBytes = TRUE)
    charset <- if(length(pos))
        sub("^[*] using session charset: *([^ ]*) *$", "\\1",
            lines[pos[1L]], useBytes = TRUE)
    else ""
    ## Re-encode to UTF-8.
    lines <- iconv(lines, from = charset, to = "UTF-8", sub = "byte")
    
    ## HTML charset:
    codepoints <- c(1L : 8L, 11L, 12L, 14L : 31L, 127L : 159L)
    lines <- gsub(sprintf("[%s]", intToUtf8(codepoints)), " ", lines,
                  perl = TRUE)
    ## HTML escapes:
    lines <- gsub("&", "&amp;", lines, fixed = TRUE)
    lines <- gsub("<", "&lt;", lines, fixed = TRUE)
    lines <- gsub(">", "&gt;", lines, fixed = TRUE)

    ## Fancy stuff:
    ind <- grep("^\\*\\*? ", lines)
    lines[ind] <- sub("(\\.\\.\\.( \\[.*\\])?) (WARNING|ERROR)$",
                      "\\1 <span class=\"boldred\">\\3</span>",
                      lines[ind])

    ## Convert pointers to install.log:
    ind <- grep("^See ['‘]http://.*['’] for details.$", lines)
    if(length(ind))
        lines[ind] <- sub("^See ['‘](.*)['’] for details.$",
                          "See <a href=\"\\1\">\\1</a> for details.",
                          lines[ind])

    ## Handle footers.
    footer <- character()
    len <- length(lines)

    ## SU seems to add refs and notes about elapsed time.
    if(grepl("^\\* elapsed time", lines[len])) {
        len <- len - 1L
        lines <- lines[seq_len(len)]
    }
    if(all(lines[c(len - 2L, len)] == c("See", "for details."))) {
        len <- len - 3L
        lines <- lines[seq_len(len)]
    }

    ## Handle NOTE/WARNING log summary footers.
    num <- length(grep("^(NOTE|WARNING): There",
                       lines[c(len - 1L, len)]))
    if(num > 0L) {
        pos <- seq.int(len - num + 1L, len)
        footer <- sprintf("<p>\n%s\n</p>",
                          paste(lines[pos], collapse = "<br/>\n"))
        lines <- lines[-pos]
    }

    if(subsections) {
        ## In the old days, we had bundles.
        ## Now, we have multiarch ...
        ## Sectioning.
        ## Somewhat tricky as we like to append closing </li> to the
        ## lines previous to new section starts, so that we can easily
        ## identify the "uninteresting" OK lines (see below).
        count <- rep.int(0L, length(lines))
        count[grep("^\\* ", lines)] <- 1L
        count[grep("^\\*\\* ", lines)] <- 2L
        ## Hmm, using substring() might be faster than grepping.
        ind <- (count > 0L)
        ## Lines with count zero are "continuation" lines, so the ones
        ## before these get a line break.
        pos <- which(!ind) - 1L
        if(length(pos))
            lines[pos] <- paste(lines[pos], "<br/>", sep = "")
        ## Lines with positive count start a new section.
        pos <- which(ind)
        lines[pos] <- sub("^\\*{1,2} ", "<li>", lines[pos])
        ## What happens to the previous line depends on whether a new
        ## subsection is started (bundles), and old same-level section
        ## or subsection is closed, or both a subsection and section are
        ## closed: these cases can be distinguished by looking at the
        ## count differences (values 1, 0, and -1, respectively).
        delta <- c(0, diff(count[pos]))
        pos <- pos - 1L
        if(length(p <- pos[delta > 0]))
            lines[p] <- paste(lines[p], "\n<ul>", sep = "")
        if(length(p <- pos[delta == 0]))
            lines[p] <- paste(lines[p], "</li>", sep = "")
        if(length(p <- pos[delta < 0]))
            lines[p] <- paste(lines[p], "</li>\n</ul></li>", sep = "")
        ## The last line always ends a section, and maybe also a
        ## subsection.
        len <- length(lines)
        lines[len] <- sprintf("%s</li>%s", lines[len],
                              if(count[pos[length(pos)] + 1L] > 1L)
                              "\n</ul></li>" else "")
    } else {
        ind <- substring(lines, 1L, 2L) == "* "
        if(!any(ind))
            lines <- character()
        ## Lines not starting with '* ' are "continuation" lines, so the
        ## ones before these get a line break.
        pos <- which(!ind) - 1L
        if(length(pos))
            lines[pos] <- paste(lines[pos], "<br/>", sep = "")
        ## Lines starting with '* ' start a new block, and end the old
        ## one unless first.
        pos <- which(ind)
        if(length(pos)) {
            ## Replace star by <li>.
            lines[pos] <- sprintf("<li>%s", substring(lines[pos], 3L))
            ## Append closing </li> to previous line unless first.
            pos <- (pos - 1L)[-1L]
            lines[pos] <- sprintf("%s</li>", lines[pos])
        }
        ## The last line always ends the last block.
        len <- length(lines)
        lines[len] <- paste(lines[len], "</li>", sep = "")
    }

    grayify <- function(lines, subsections = FALSE) {
        ## Turn all non-noteworthy parts into gray.

        ## Handle 'checking extension type ... Package' directly.
        ind <- lines == "<li>checking extension type ... Package</li>"
        if(any(ind))
            lines[ind] <-
                "<li class=\"gray\">checking extension type ... Package</li>"
        ## Same for 'DONE'
        ind <- lines == "<li>DONE</li>"
        if(any(ind))
            lines[ind] <-
                "<li class=\"gray\">DONE</li>"

        foo_simple <- function(lines) {
            chunks <-
                split(lines, cumsum(substring(lines, 1L, 4L) == "<li>"))
            unlist(lapply(chunks, function(s) {
                s <- paste(s, collapse = "\n")
                sub("^<li>( *([^\n]*)\\.\\.\\.( \\[.*\\])? OK(<br/>.*)?)</li>",
                    "<li class=\"gray\">\\1</li>",
                    s)
            }),
                   use.names = FALSE)
        }

        foo_tricky <- function(lines) {
            chunks <-
                split(lines, cumsum(substring(lines, 1L, 3L) == "<li"))
            unlist(lapply(chunks, function(s) {
                s <- paste(s, collapse = "\n")
                s <- sub("^<li class=\"black\">( *([^\n]*)\\.\\.\\.( \\[.*\\])? OK(<br/>.*)?)</li>\n/ul></li>$",
                         "<li class=\"gray\">\\1</li>\n</ul></li>",
                         s)
                sub("^<li class=\"black\">( *([^\n]*)\\.\\.\\.( \\[.*\\])? OK(<br/>.*)?)</li>",
                    "<li class=\"gray\">\\1</li>",
                    s)
            }),
                   use.names = FALSE)
        }
        
        if(!subsections) {
            foo_simple(lines)
        } else {
            ## Determine lines starting and ending subsections.
            ind_ss_s <- grepl("\n<ul>$", lines)
            ind_ss_e <- grepl("</li>\n</ul></li>$", lines)
            ## The former (currently?) give subsection titles only, and
            ## hence always get grayified.  However, apparently this
            ## results in nested <li> to be grayified as well: hence,
            ## tag everthing as black first.
            lines[ind_ss_s] <-
                sub("<li>(.*)(\n<ul>)$",
                    "<li class=\"gray\">\\1\\2",
                    lines[ind_ss_s])
            lines <- sub("^<li>", "<li class=\"black\">", lines)
            ## Split into subsection-related blocks.
            blocks <- split(lines,
                            cumsum(ind_ss_s +
                                   c(0L, head(ind_ss_e, -1L))))
            unlist(lapply(blocks, foo_tricky), use.names = FALSE)
        }
    }

    lines <- grayify(lines, subsections)

    ## Header.
    writeLines(c("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">",
                 "<html xmlns=\"http://www.w3.org/1999/xhtml\">",
                 "<head>",
                 sprintf("<title>Check results for '%s'</title>",
                         sub("-00check.(log|txt)$", "", basename(log))),
                 "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"/>",
                 "<link rel=\"stylesheet\" type=\"text/css\" href=\"../R_check_log.css\"/>",
                 "</head>",
                 "<body>"),
               out)
    ## Body.
    if(!length(lines))
        writeLines("<p>\ncheck results unavailable\n</p>", out)
    else
        writeLines(c("<ul>", lines, "</ul>", footer), out)
    ## Footer.
    writeLines(c("</body>",
                 "</html>"),
               out)
}


## Search R packages in directory.
search_Rpkgs <- function(dir = ".")
{
  files <- dir(dir, full.names = TRUE, include.dirs = TRUE)
  packages <- NULL
  for(f in files) {
    if(file.info(f)$isdir) {
      ff <- dir(f, include.dirs = TRUE)
      if(all(c("DESCRIPTION", "R", "NAMESPACE", "man") %in% ff))
        packages <- c(packages, f)
    }
  }
  return(packages)
}


## Run build, install or check.
run_R_bic <- function(package, R = "R", check = FALSE, suggests = FALSE, rm = FALSE)
{
  owd <- getwd()
  on.exit(setwd(owd))
  dir.create(tdir <- tempfile())
  on.exit(unlink(tdir), add = TRUE)
  depends <- readLines(file.path(package, "DESCRIPTION"))
  if(length(i <- grep("Depends:", depends)) > 0L) {
    depends <- depends[i]
    depends <- gsub("Depends:", "", depends, fixed = TRUE)
    depends <- strsplit(depends, ",")[[1]]
    depends <- depends[!(grepl("R", depends) & grepl("(", depends, fixed = TRUE))]
    if(length(depends) > 0L) {
      depends <- gsub(" ", "", depends)
      lib_script <- '
        lp <- grep("home", .libPaths(), fixed = TRUE, value = TRUE)
        writeLines(lp, "libPaths.txt")
      '
      writeLines(lib_script, "lib_script.R")
      cmd <- paste(R, "CMD BATCH --no-save lib_script.R")
      system(cmd)
      libs <- readLines("libPaths.txt")
      install <- list()
      for(j in depends) {
        for(i in libs) {
          install[[j]] <- c(install[[j]], j %in% dir(i))
        }
      }
      install <- sapply(install, function(x) { !any(x) })
      inst_script <- NULL
      for(j in names(install)) {
        if(install[j]) {
          inst_script <- c(inst_script, paste('install.packages("', j, '")', sep = ''))
        }
      }
      if(!is.null(inst_script)) {
        writeLines(inst_script, "inst_script.R")
        cmd <- paste(R, "CMD BATCH --no-save inst_script.R")
        system(cmd)
      }
    }
  }
  setwd(tdir)
  if(rm) {
    cmd <- paste(R, "CMD REMOVE", basename(package), ">> /dev/null 2>&1")
    try(system(cmd))
  }
  cmd <- paste(R, "CMD build", package, ">> pkg_build.log 2>&1")
  ok <- !inherits(try(system(cmd)), "try-error")
  logs <- list("build" = if(ok) readLines("pkg_build.log") else NULL)
  tar <- grep("tar.gz", dir(tdir), value = TRUE)
  if(!check) {
    cmd <- paste(R, "CMD INSTALL", tar, ">> pkg_install.log 2>&1")
    system(cmd, intern = FALSE)
    logs[["install"]] = readLines("pkg_install.log")
  } else {
    cmd <- paste(R, "CMD check --as-cran", tar, ">> /dev/null 2>&1")
    if(!suggests)
      cmd <- paste("export _R_CHECK_FORCE_SUGGESTS_=FALSE;", cmd)
    system(cmd, intern = FALSE)
    check_dir <- paste(basename(package), "Rcheck", sep = ".")
    logs[["check"]] = readLines(file.path(check_dir, "00check.log"))
    logs[["install"]] = readLines(file.path(check_dir, "00install.out"))
    if(file.exists(test_dir <- file.path(check_dir, "tests"))) {
      if(length(Rout <- grep(".Rout", dir(test_dir), value = TRUE))) {
        logs[["tests"]] <- list()
        for(j in Rout)
          logs$tests[[j]] <- readLines(file.path(test_dir, j))
      }
    }
  }
  return(logs)
}


## Function to run checks on all packages in dir
## with all flavors of R.
run_R_checks <- function(dir = ".",
  R = c("/var/tmp/BayesR/R-clang/install/bin/R", "/var/tmp/BayesR/R-gcc/install/bin/R"))
{
  packages <- search_Rpkgs(dir)
  packages <- packages[-grep("dev", packages)]
  if(!length(packages)) stop(paste("cannot find any packages in", dir, "for checking!"))
  res <- list()
  for(pkg in packages) {
    cat(paste("** checking package", basename(pkg), "\n"))
    res[[basename(pkg)]] <- list()
    for(rf in R) {
      cat("**** using R flavor", rf, "\n")
      res[[basename(pkg)]][[rf]] <- run_R_bic(pkg, R = rf, check = TRUE)
    }
  }
  return(res)
}


## Create html ckeck results.
checks2html <- function(x, dir = ".", name = NULL)
{
  owd <- getwd()
  setwd(dir)
  on.exit(setwd(owd))
  dir2 <- paste(name, "logs", sep = "_")
  if(!file.exists(dl <- file.path(dir, dir2)))
    dir.create(dl)
  packages <- names(x)
  tab <- c(
    '<table>',
    '<tr>',
    '<th>Package</th>',
    '<th>R flavor</th>',
    '<th>Build</th>',
    '<th>Install</th>',
    '<th>Check</th>',
    '<th>Tests</th>',
    '</tr>'
  )
  for(pkg in packages) {
    pkg2 <- basename(pkg)
    tab <- c(tab, '<tr>',
      paste('<td rowspan="', length(x[[pkg]]), '"><font size="4"><code>',
        pkg, '</code></font></td>', sep = ''))
    k <- 1
    for(rf in names(x[[pkg]])) {
      if(k > 1) tab <- c(tab, '<tr>')
      tab <- c(tab,
        paste('<td>', rf, '</td>', sep = ''),
        make_check_cell(x[[pkg]][[rf]]$build, name = paste(pkg2, paste("R", k, sep = ""), "build", sep = "_"), dir = dir2),
        make_check_cell(x[[pkg]][[rf]]$install, name = paste(pkg2, paste("R", k, sep = ""), "install", sep = "_"), dir = dir2),
        make_check_cell(x[[pkg]][[rf]]$check, name = paste(pkg2, paste("R", k, sep = ""), "check", sep = "_"), dir = dir2, check = TRUE),
        make_check_cell(x[[pkg]][[rf]]$tests, name = paste(pkg2, paste("R", k, sep = ""), "tests", sep = "_"), dir = dir2)
      )
      if(k > 1) tab <- c(tab, '</tr>')
      k <- k + 1
    }
    tab <- c(tab, '</tr>')
  }
  tab <- c(tab, '</table>')

  html <- c(
    '<!DOCTYPE html>',
    '<html>',
    '<head>',
    '<style>',
    'table, th, td {',
    'border: 1px solid black;',
    'border-collapse: collapse;',
    '}',
    'th, td {',
    'padding: 5px;',
    'text-align: left;',  
    '}',
    'div#tab {',
    'margin-top: 2%;',
    'margin-left: 2%;',
    '}',
    '</style>',
    '</head>',
    '<body>',
    '<div id="tab">',
    paste('<h2>R-devel package checks ', Sys.Date(), '</h2>', sep = ''),
    tab,
    '</div>',
    '</body>'
  )

  writeLines(html, file.path(dir, paste(name, "pkg_checks.html", sep = "_")))
}


## Create a check cell.
make_check_cell <- function(x, name, dir, check = FALSE) {
  for(char in c(".", "/", ":"))
    name <- gsub(char, "", name, fixed = TRUE)
  if(is.null(x)) return('<td><strong>na</strong></td>')
  if(!is.list(x)) {
    if(length(i <- grep("WARNING: ignoring environment value of R_HOME", x, fixed = TRUE)))
      x <- x[-i]
  }
  x0 <- x
  x <- unlist(x)
  note <- if(check) {
    any(grepl("NOTE:", x, fixed = TRUE)) | any(grepl("Note:", x, fixed = TRUE))
  } else any(grepl("note", x, ignore.case = TRUE))
  warn <- if(check) {
    any(grepl("WARNING:", x, fixed = TRUE)) | any(grepl("Warning:", x, fixed = TRUE))
  } else any(grepl("warning", x, ignore.case = TRUE))
  error <- if(check) {
    any(grepl("ERROR:", x, fixed = TRUE)) | any(grepl("Error:", x, fixed = TRUE))
  } else any(grepl("error:", x, ignore.case = TRUE))
  ok <- !any(c(note, warn, error))
  status <- c("ok", "note", "warn", "error")[c(ok, note, warn, error)]
  status <- status[length(status)]
  href <- file.path(dir, paste(name, "html", sep = "."))
  if(check) {
    tf <- tempfile()
    writeLines(x0, tf)
    x0 <- write_check_log_as_HTML(tf, out = href)
  } else {
    x0l <- FALSE
    if(is.list(x0)) {
      x0l <- TRUE
      x2 <- NULL
      for(j in names(x0))
        x2 <- c(x2, '<p>', paste("<h2>", j, '</h2>', sep = ""), paste(x0[[j]], '<br>'), '</p>')
      x0 <- x2
    }
    x0 <- c(
      '<!DOCTYPE html>',
      '<html>',
      '<head>',
      '</head>',
      '<body>',
      paste(x0, if(!x0l) '<br>' else NULL),
      '</body>'
    )
    writeLines(x0, con = href)
  }
  acol <- switch(status,
    "ok" = "009900",
    "note" = "FFCC00",
    "warn" = "FF6600",
    "error" = "FF0000"
  )
  acol <- paste('#', acol, sep = '')
  cell <- paste('<a href="', href,'" style="color: ', acol, '; font-weight: bold">', status, '</a>', sep = '')
  cell <- paste('<td>', cell, '</td>', sep = '')

  return(cell)
}


## Check the BayesR project.
check_BayesR <- function(dir = ".",
  R = c("/var/tmp/BayesR/R-clang/install/bin/R", "/var/tmp/BayesR/R-gcc/install/bin/R"),
  load = FALSE)
{
  owd <- getwd()
  if(dir == ".") dir <- owd
  tdir <- tempfile()
  dir.create(tdir)
  on.exit(unlink(tdir))
  setwd(tdir)
  on.exit(setwd(owd), add = TRUE)
  Rok <- NULL
  writeLines('print(runif(10))', 'Rcheck.R')
  for(j in R) {
    cmd <- paste(j, "CMD BATCH --no-save Rcheck.R >> /dev/null 2>&1")
    log <- system(cmd, intern = FALSE)
    if(log < 1)
      Rok <- c(Rok, j)
  }
  R <- Rok
  if(!length(R)) stop("no valid R installations!")
  if(file.exists("~/R_pkg_checks.rda") & load) {
    load("~/R_pkg_checks.rda")
  } else {
    cmd <- "svn checkout svn://scm.r-forge.r-project.org/svnroot/bayesr/"
    system(cmd)
    if(any(grepl("BayesXsrc", dir(file.path(tdir, "bayesr/pkg"))))) {
      setwd(file.path(tdir, "bayesr/pkg/BayesXsrc"))
      cmd <- "sh ./bootstrap-devel.sh"
      system(cmd)
      setwd(tdir)
    }
    system("rm -rf BayesXdev")
    R_check_logs <- run_R_checks(dir = file.path(tdir, "bayesr/pkg"), R = R)
    save(R_check_logs, file = "~/R_pkg_checks.rda")
  }
  dir.create(www <- file.path(tdir, "www"))
  checks2html(R_check_logs, dir = www, name = "BayesR")
  setwd(www)
  cmd <- paste("cp -r *", dir)
  system(cmd)
#  cmd <- paste("scp -i /home/c403129/.ssh/id_dsa -r *", "eeecon:public_html/Rpkgs/")
#  ok <- try(system(cmd), silent = TRUE)
#  if(inherits(ok, "try-error")) {
#    cmd <- paste("scp -i /home/nik/.ssh/id_dsa -r *", "eeecon:public_html/Rpkgs/")
#    system(cmd)
#  }
  return(invisible(NULL))
}


## Show results in browser.
showChecks <- function(url = "http://eeecon.uibk.ac.at/~umlauf/Rpkgs/BayesR_pkg_checks.html")
{
  browseURL(url)
}

build_R <- function(...)
{
  dir <- "/var/tmp/BayesR"
  owd <- getwd()
  on.exit(setwd(owd))
  if(!file.exists(dir))
    dir.create(dir)
  setwd(dir)

  ## Checkout sources.
  if(!file.exists("R-gcc")) {
    cat("Checking out R sources for gcc..\n")
    script <- c(
      "#!/bin/bash",
      "cd /var/tmp/BayesR",
      "mkdir R-gcc",
      "mkdir R-gcc/install",
      "mkdir R-gcc/log",
      "cd R-gcc",
      "svn checkout https://svn.r-project.org/R/trunk/ src > /dev/null"
    )
    writeLines(script, "Rgcc.sh")
    system("sh ./Rgcc.sh")
    file.remove("Rgcc.sh")
  }

  if(!file.exists("R-clang")) {
    cat("Checking out R sources for clang..\n")
    script <- c(
      "#!/bin/bash",
      "cd /var/tmp/BayesR",
      "mkdir R-clang",
      "mkdir R-clang/install",
      "mkdir R-clang/log",
      "cd R-clang",
      "svn checkout https://svn.r-project.org/R/trunk/ src > /dev/null"
    )
    writeLines(script, "Rclang.sh")
    system("sh ./Rclang.sh")
    file.remove("Rclang.sh")
  }

  setwd(dir)

  ## Build R with gcc.
  cat("Building R-gcc.\n")
  script <- c(
    '#!/bin/bash',
    'cd /var/tmp/BayesR/R-gcc/src/',
    'rm -f /var/tmp/BayesR/R-gcc/log/build.out /var/tmp/BayesR/R-gcc/log/build.err',
    'make uninstall >> /var/tmp/BayesR/R-gcc/log/build.out 2>> /var/tmp/BayesR/R-gcc/log/build.err',
    'svn update > /dev/null',
    'tools/rsync-recommended > /dev/null',
    'export CC="gcc -std=gnu99 -fsanitize=address -fno-omit-frame-pointer"',
    'export CXX="g++ -fsanitize=address -fno-omit-frame-pointer"',
    'export F77="gfortran -fsanitize=address"',
    'export FC="gfortran -fsanitize=address"',
    './configure --prefix /var/tmp/BayesR/R-gcc/install --enable-R-shlib> /var/tmp/BayesR/R-gcc/log/build.out 2> /var/tmp/BayesR/R-gcc/log/build.err',
    'make #>> /var/tmp/BayesR/R-gcc/log/build.out 2>> /var/tmp/BayesR/R-gcc/log/build.err',
    'make install >> /var/tmp/BayesR/R-gcc/log/build.out 2>> /var/tmp/BayesR/R-gcc/log/build.err'
  )
  writeLines(script, "build.R-gcc.sh")
  system("sh ./build.R-gcc.sh >> /dev/null")

  ## Build R with clang.
  cat("Building R-clang.\n")
  symbolizer <- grep("llvm-", dir("/usr/lib/"), fixed = TRUE, value = TRUE)
  script <- c(
    '#!/bin/bash',
    'cd /var/tmp/BayesR/R-clang/src/',
    'rm -f /var/tmp/BayesR/R-clang/log/build.out /var/tmp/BayesR/R-clang/log/build.err',
    'make uninstall >> /var/tmp/BayesR/R-clang/log/build.out 2>> /var/tmp/BayesR/R-clang/log/build.err',
    'svn update > /dev/null',
    'tools/rsync-recommended > /dev/null',
    'export CC="clang -fsanitize=undefined,address -fno-sanitize=float-divide-by-zero -fno-omit-frame-pointer"',
    'export CXX="clang++ -fsanitize=undefined,address -fno-sanitize=float-divide-by-zero -fno-omit-frame-pointer"',
    if(length(symbolizer)) {
      paste0('export ASAN_SYMBOLIZER_PATH="/usr/lib/', symbolizer[1], '/bin/llvm-symbolizer"')
    } else NULL,
    './configure --prefix /var/tmp/BayesR/R-clang/install --enable-R-shlib> /var/tmp/BayesR/R-clang/log/build.out 2> /var/tmp/BayesR/R-clang/log/build.err',
    'make #>> /var/tmp/BayesR/R-clang/log/build.out 2>> /var/tmp/BayesR/R-clang/log/build.err',
    'make install >> /var/tmp/BayesR/R-clang/log/build.out 2>> /var/tmp/BayesR/R-clang/log/build.err'
  )
  writeLines(script, "build.R-clang.sh")
  system("sh ./build.R-clang.sh >> /dev/null")
}


if(FALSE) {
  build_R()
  if(!file.exists("/var/tmp/BayesR/www"))
    dir.create("/var/tmp/BayesR/www")
  check_BayesR(
    dir = "/var/tmp/BayesR/www",
    R = c("/var/tmp/BayesR/R-clang/install/bin/R", "/var/tmp/BayesR/R-gcc/install/bin/R")
  )
}

