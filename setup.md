---
layout: page
title: Setup 
permalink: /setup/
---

# Software Requirements & Setup Instructions
To fully participate in this boot camp, you will need access to the software described below on **your own laptop** (N.B. You will likely need _Administer privileges/permissions_ to install some of these). 

**Access to the KOBIC computer cluster:**
In addition to the software listed further below, you will require
access to the KOBIC computer cluster. Each user will be informed
the login information to the KOBIC server by the instructors. 

<br>

### Terminal Utilities
Bash is a commonly-used shell that gives you the power to do simple tasks more quickly.

**Windows:** Install [MobaXterm](http://mobaxterm.mobatek.net), an enhanced terminal with bash for Windows. Note that the default 'Personal Edition install' typically places the MobaXterm executable in `C:\Program Files (x86)\Mobatek\MobaXterm Personal Edition`.  

Please also install the [**Plugin CygUtils**](http://mobaxterm.mobatek.net/CygUtils.plugin). Once downloaded please move the `CygUtils.plugin` file to the folder `C:\Program Files (x86)\Mobatek\MobaXterm Personal Edition`. Launching mobaxterm will complete the install.  

**Mac OS X:** You do not need to install anything. You can access bash from the **Terminal** (found in **/Applications/Utilities**). You may want to keep Terminal in your dock for this class.

**Linux:** There is no need to install anything.

<br>

### Text Editor
When you're writing code, it's nice to have a text editor that is optimized for writing code, with features like automatic color-coding of key words. The default text editor on Mac OS X and Linux is usually set to Vim, which is not famous for being intuitive. if you accidentally find yourself stuck in it, try typing the **escape key**, followed by **:q!** (colon, lower-case 'q', exclamation mark), then hitting Return to return to the shell.

**Windows:** nano is a basic editor and the default that instructors use in the workshop. Nano should be installed as a plugin to mobaxterm (see above **Git** instructions for windows).

**Mac OS X:** nano should be pre-installed.

**Linux:** nano should be pre-installed.


<br>

### R and RStudio
R Binaries for Windows, MacOSX and Linux can be downloaded and
installed from [CRAN](http://cran.r-project.org/index.html)
(Comprehensive R Archive Network). If possible download the latest
binary version of R for your operating system. As of course launch
(Aug 16) the latest release (2018/07/02, "Feather Sprat") is R-3.5.1.

After installing R itself we recommend installing the *preview* version of [RStudio](https://www.rstudio.com/products/rstudio/download/preview/) desktop (v1.1.453 or close), a slick visual interface for R.



<br>

### Python
Python is a popular language for scientific computing. Installing all of its scientific packages individually can be a bit difficult, so we recommend an all-in-one installer.

Regardless of how you choose to install it, please make sure you install Python version 2.x and not version 3.x (e.g., 2.7 is fine but not 3.4). Python 3 introduced changes that will break some of the code we teach during the workshop.

We will teach Python using the IPython notebook, a programming environment that runs in a web browser. For this to work you will need a reasonably up-to-date browser. The current versions of the Chrome, Safari and Firefox browsers are all supported (some older browsers, including Internet Explorer version 9 and below, are not).

**Windows:** Download and install [Anaconda](http://continuum.io/downloads.html).
Download the default *Python 2 graphical installer* installer (do not follow the link to version 3). Use all of the defaults for installation except make sure to check _Make Anaconda the default Python_.

**Mac OS X:** Download and install [Anaconda](http://continuum.io/downloads.html).
Download the default "Mac OS X **Python 2.7** Graphical Installer" (do not follow the link to version 3). Use all of the defaults for installation.

**Linux:** As above, we recommend the all-in-one scientific Python installer for linux [Anaconda](http://continuum.io/downloads.html).
