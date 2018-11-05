# Using this Package

## System Requirements
This project makes use of the STAR aligner by Alex Dobin to sequence human
genome data.  STAR requires more than 30G of RAM for this purpose.  In general
genome sequencing carries relatively high demands for compute, RAM, and
disk I/O resources (but not for graphics resources, at least not with the
software used in this project).

This project also makes use of various resource-intensive R packages.

The following systems were used during development of this work:

1. System76 "Serval" mobile workstation
    * Intel i7-8700k CPU
    * 64GB RAM (DDR4-3000, non-ECC)
    * Samsung NVMe Pro SSD
    * Pop!_OS 17.10
2. PowerSpec G460 desktop computer
    * Intel i7-8700k CPU
    * 48GB RAM (DDR4-3200, non-ECC)
    * Intel M.2 SATA SSD
    * Windows 10 Pro ; Windows Subsystem for Linux, "WLinux" flavor

Note that WSL is still evolving and rough around the edges.  For example,
the "atom" editor has some bugs that prevent plug-in installation under WSL.

In Q4 2018, Apple announced a new Mac Mini lineup.  A fully upgraded Mini
would appear to have comparable specs and in theory should be suitable, but
no attempt has been made to evaluate this.


## Installing Prerequisite Software Packages
This software assumes that various tools are already available in your
environment.  If you have not done so already, it will be useful for you to do
the following (for Debian-based environments):

* sudo apt-get install git-lfs  # needed for large .fastq raw data files
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


