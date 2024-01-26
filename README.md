# 51EriProject
Project description. Two sentences about the goal of the project and the paper that uses this code.

Link to paper here: [Paper](https://arxiv.org/abs/2401.01468)

PARSEC Isochrone website: [CMD 3.7 input form](http://stev.oapd.inaf.it/cgi-bin/cmd)
# Required Software
## Git
Git and GitHub are standard tools across astronomy, computer science, computer engineering,
and other fields. This is a toolset that everyone should learn for code management. It is 
always useful and often expected for contributors to a project of any size.

https://git-scm.com/downloads

This installation will include _gitbash_ on Windows machines. This is a terminal application
that simulates a unix terminal on Windows. It is useful for managing code and working with
git and github. If this is the first terminal/shell you have installed on your
Windows computer, you will be able to use this to enter all of the terminal
commands present in the rest of this README.md file.

A 15-minute introduction video for *git* is available at  https://www.youtube.com/watch?v=USjZcfj8yxE


## Python3
Get the latest version of Python3

https://www.python.org/downloads/

This project was tested with Python 3.8.5 to 3.10.x. It probably works with older
versions of Python, but it is not guaranteed.

Some python 
installations will need to call `python3` instead of `python` from the terminal. 
Check your python version with `python --version` if it comes back with 
2.7.x, try `python3 --version` and expect to get 3.x.x. 

Modern Python installations require you to sign an SSL certificate after
installation. This is done by simply running a script in the installation 
directory. Programs like, [Homebrew](https://brew.sh/) will do this step for
you. If you see SSL errors in Python, it is probably because you skip the
certificate step.

## Installation

### Clone the repository, and navigate to the directory.

```
git clone https://github.com/aaelli99/51EriProject
cd 51EriProject
```

### Remember to update the repository periodically with:

`git pull`

### Virtual Environment Setup and Activation (recommended)

Configure the virtual environment for DetMap or from the terminal.

### Virtual Environment Setup

This step is done one time to initialize the virtual environment.

Window and Linux Tested in Windows PowerShell and OSX Terminal. Some Python 
installations will need to call `python3` instead of `python`. Check your Python version with
`python --version` if it comes back with 2.7.x, try `python3 --version` and expect to get 3.x.x. 

If you have many Python installations, or if you just downloaded the latest version
of python then you should replace `python` with the full path to the python executable.

```
python -m pip install --upgrade pip
pip install virtualenv
virtualenv venv
```

### Activating a Virtual Environment
This step needs to be done every time you want to use the virtual environment.

For unix-like environment activation:

```source venv/bin/activate```

Windows CMD environment activation:

```.\venv\Scripts\activate```

for Windows powershell environment activation:

```.\venv\Scripts\activate.ps1```

After activation, the term ail will add a `(venv)` to the command prompt. At this point
you can drop the full paths to the Python and pip executables into the terminal, 
and the `python3` in place of `python` commands.

I like test that everything went as expected with the following:

```
python --version
pip list
