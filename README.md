---
output:
  html_document:
    toc: yes
---


# Using this Package

## System Requirements
This project makes use of the [STAR aligner](https://github.com/alexdobin/STAR)
to sequence human genome data.  STAR requires more than 30G of RAM for this
purpose.  In general genome sequencing carries relatively high demands for
compute, RAM, and disk I/O resources (but not for graphics resources, at
least not with the software used in this project).

This project also makes use of various resource-intensive R packages.

Nonetheless, the necessary hardware is attainable in high-end
consumer-grade systems.

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

System (1) is capable, starting from scratch, to download and compile
the STAR source code, download human genome data, align a reference
human genome, and translate 8 .fastq.gz files into .bam counterparts,
all in about 75 minutes.


## Installing Prerequisite Software Packages
This software assumes that various tools are already available in your
environment.  If you have not done so already, it will be useful for you to do
the following (for Debian-based environments):

* `sudo apt-get install libssl-dev  # needed by various R packages`
* `sudo apt-get install git-lfs  # needed for large .fastq raw data files`
* `sudo apt-get install samtools`
* `sudo apt-get install bamtools`
* `sudo apt-get install python-pip  # (note, use python 2.7)`
* `sudo apt-get install r-base r-recommended`

## Cloning this Package from Github

In a directory of your choosing, execute:

`git clone git@github.com:a3609640/DNFA-gene-analysis.git`

Note, both here and in "Configuring Your Local Environment", that you
may wish to consider locations carefully.  In particular, you should
make sensible plans for file system free space, and for underlying disk
speed.

## Configuring Your Local Environment
In the top-level directory of your Git clone, there is a file named
'configuration.tmpl'.  You should copy this file locally to one named
'configuration', edit the contents as you prefer (to indicate local directories
where generated files will be stored, and the like), and execute the commmand
`source configuration` before using the scripts in this package.  That will
ensure that necessary shell environment variables are set.

## First Steps
(Note: References in this section to STAR-Index-and-Alignment may soon be
replaced by instructions to `make bamfiles`.)

From the top-level directory of your Git project clone, the first script to
invoke, after `source configuration`, is `RNA-Seq/STAR-Index-and-Alignment`.
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

# Software Engineering Discussion for Biologists

## Build Systems

One aspect of this project is a file named `Makefile`.  To a non-initiate,
its contents undoubtedly appear as though, periodically during its creation,
a cat walked across a keyboard with a stuck `shift` key, taking particular
care to step only on the number keys.  However the file in fact specifies
many specific buildable units, and the relationships between them.

### A Primer on Makefiles

Putting aside the meanings of various symbols, the basic idea is
straightforward.  The Makefile provides a set of rules with the following
form (as described in the
[documentation for GNU Make](https://www.gnu.org/software/make/manual/html_node/Rule-Introduction.html#Rule-Introduction)):
```
target : prerequisites
    recipe
```

The basic idea is: to 'make' a target, first verify the existence
of its prerequisites, then execute a recipe.

The value starts to become more apparent with just a slightly adjusted example:
```
target2 : target1
    recipe2

target1 :
    recipe1
```

In this case, `make target2` will check for the existence of target1.  If
target1 exists, then 'make' will execute recipe2.  If target1 does not exist,
'make' will notice that target1 is *itself* a target.  Since target1 has no
prerequisites, it will execute recipe1 and thereby generate target1.  Then,
with target2's prerequisites satisfied, it will execute recipe2.

In this way it is possible to specify a detailed protocol to produce a
particular target.  When multiple dependencies exist, 'make' can keep track of
which ones have been updated, and how they affect the rest of the target
hierarchy, and rebuild other targets as needed.  When targets are well
specified, it is often possible to delete and regenerate a very specific output
in isolation.

### The Wider World of Build Systems

In fact, 'make' is just one of many things like 'make'.  It is one of the
oldest, and its Makefile syntax produces some of the most arcane project
descriptions.

Alternatives include Python's 'snakemake' and Google's 'bazel'
([don't call it 'blaze'](https://bazel.build/faq.html#whats-up-with-the-word-blaze-in-the-codebase)).

This project uses 'make' because it is very well-established.  Users
of this package most likely will not need to undertake any special steps
to install it, nor does the package itself need to make any provision for
such installation.  The downside is notoriously arcane Makefile syntax, but
that is mostly a burden for project contributors rather than users.  Without
caring about the actual contents of the Makefile, users can easily
understand and execute an instruction like this:

`make /opt/dnfa_genfiles/data/Analysis/Samsort/test2_S3_L004Aligned.sorted.bam`.

If you consider using 'make' for one of your own projects, good advice
would be: don't start by copying this project's Makefile and trying to
adjust it for your purposes.  Rather, just start out with some very simple
targets, then start learning and using automatic variables, one or two at
a time. Eventually you will find that your project has a gibberish Makefile
too, but to your surprise it will work, and to your even greater surprise
you will understand it. (Well, that's not a promise, on either count.  YMMV.)

### Features Provided by This Project's Makefile
.... <This section in progress> ...
