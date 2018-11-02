# Using this Package

## Installing Prerequisite Software Packages
This software assumes that various tools are already available in your
environment.  If you have not done so already, it will be useful for you to do
the following (for Debian-based environments):

* sudo apt-get install samtools
* sudo apt-get install bamtools
* sudo apt-get install python-pip  # (note, use python 2.7)
* sudo apt-get install r-base r-recommended

## Configuring Your Local Environment
In the top-level directory of your Git clone, there is a file named
'configuration.tmpl'.  You should copy this file locally to one named
'configuration', edit the contents as you prefer (to indicate local directories
where generated files will be stored, and the like), and execute the commmand
'source configuration' before using the scripts in this package.  That will
ensure that necessary shell environment variables are set.

## First Steps
From the top-level directory of your Git project clone, the first script to
invoke, after 'source configuration', is 'RNA-Seq/STAR-index-and-Alignment'.
That script will download the sources for a specific version of STAR, and
compile them on your machine.  Then it will download and sequence a reference
human genome.  Then it will process various project-specific files against That
reference genome.

The STAR-index-and-Alignment script attempts to Do The Right Thing if it is
invoked multiple times.  It creates output directories with sensible permissions
if they are absent.  If STAR is already built and installed it does not repeat
the download and compilation.  If human genome files are already downloaded it
does not download them again.  If a reference genome has already been sequenced,
it is reused.  Etc.  Thus, it should be safe to re-invoke the script if it was
interrupted.  However note that the detection of previously-completed work is not
foolproof.  If in doubt, it is best to rm -rf ${generatedDataRoot}/* and then
re-execute.


