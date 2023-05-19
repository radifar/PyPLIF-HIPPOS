Installation
============

Quick Installation (Recommended)
--------------------------------

The easiest way to install HIPPOS is using `Conda <https://docs.anaconda.com/anaconda/install/>`_, you
can choose either Anaconda or Miniconda. If you never used Conda before then most likely Miniconda
is more suitable for you.

After you installed Conda in your machine, install HIPPOS using the following command:

``conda install -c conda-forge pyplif-hippos``

Conda will deal with the dependencies like Open Babel, Bitarray, and Numpy libraries. Therefore you
won't need the extra step below.

To check if HIPPOS installed correctly try these commands:

|  ``hippos``
|  ``hippos-genref``

If HIPPOS installed correctly you should get a message that inform you the configuration file not found.
Next, you can just jump to Getting Started with PLANTS or Vina tutorial. But if you prefer not to use Conda
then you can install HIPPOS using the instructions below.

Manual Installation
-------------------

The following instruction is for manual installation step which is more tedious, but you can use it in case
you would like to use the latest feature PyPLIF-HIPPOS by installing the newest code from Github.
Another advantage of manual installation is the visibility of PyPLIF-HIPPOS which makes it easier to find,
in case you would like to do some modification to PyPLIF-HIPPOS code.

Requirement
-----------

* python >= 3.5
* python3-openbabel >= 2.1
* python3-bitarray
* python3-numpy >= 1.3

Getting HIPPOS
--------------

You can get HIPPOS at Github by cloning the repository or download the code `here <https://github.com/radifar/PyPLIF-HIPPOS>`_.

Installing on Linux
-------------------

The following instructions were tested and work well on Ubuntu 20.04, Debian 11 and Fedora 30,
so they will most likely work in newer version too.

First make sure that you have Python 3 installed, and it is your default when you run `python`
command. If it is not then you can install ``python-is-python3`` package or use ``update-alternatives``
command:

| ``sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 1``

The next step is to install the requirements for PyPLIF-HIPPOS. And depending on your OS or
package manager that you use, the command could be different.

In Ubuntu or Debian you can simply enter these commands to install the requirements above:

| ``sudo apt install python3-openbabel``
| ``sudo apt install python3-bitarray``
| ``sudo apt install python3-numpy``

In Fedora you can enter these commands instead:

| ``sudo dnf install python3-openbabel``
| ``sudo dnf install python3-bitarray``
| ``sudo dnf install python3-numpy``

After all of the dependencies installed, you can install HIPPOS by opening
the terminal and enter the PyPLIF-HIPPOS directory and run setup.sh:

``./setup.sh``

This will run the script that will create an `alias` for the PyPLIF-HIPPOS
program in this directory. So every time you call `hippos` or `hippos-genref` it
will refer to this directory. Therefore it is advised **not to delete this directory**.

Notice that this process does not add anything to your machine, so this step is
just registering a shortcut to your command line environment.

If HIPPOS installed successfully then a message like 'PyPLIF HIPPOS is ready
for service' should appear.

If you would like to install newer version of HIPPOS and overwrite the old
one then all you need to do is to run the setup.sh script again.

To test if HIPPOS had been installed, *open new command line window* and type the following:

``hippos``

Then press enter. Please note that it is imperative to open the new command
line window right after the installation to allow the alias in your system get updated.
Next type the following:

``hippos-genref``

Then press enter

If there are no error message then the installation is success.

Uninstalling
------------

To uninstall PyPLIF-HIPPOS you can simply run setup.sh script with `uninstall`
option

``./setup.sh uninstall``

And a message 'PyPLIF HIPPOS has been successfully uninstalled' will appear.

Notice that as the installation process is just registering alias for `hippos`
and `hippos-genref`, the uninstallation process is actually just removing those
alias and therefore nothing is really removed from your machine. If you would
like to remove PyPLIF-HIPPOS you need to manually remove the PyPLIF-HIPPOS from
the directory where you install it.

Installing on Windows 10
------------------------

Windows 10 provides Ubuntu in their `Microsoft Store <https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6>`_. 
After you successfully add Ubuntu on your Windows, the next step is to do the same installation steps on Ubuntu.

