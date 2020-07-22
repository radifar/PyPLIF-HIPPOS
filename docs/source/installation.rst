Installation
============


Requirement
-----------

* Python >= 2.6 or >= 3.6
* Open Babel (library) >= 2.2
* Python-OpenBabel >= 2.1
* Python-BitArray
* Python-Numpy >= 1.3

Getting HIPPOS
--------------

You can get HIPPOS at Github by cloning the repository or download the code `here <https://github.com/radifar/PyPLIF-HIPPOS>`_.

Installing on Linux
-------------------

In Ubuntu you can simply enter these commands to install the requirements above (Python is already installed in Linux):

| ``sudo apt-get install openbabel``
| ``sudo apt-get install python-openbabel``
| ``sudo apt-get install python-bitarray``
| ``sudo apt-get install python-numpy``

In Fedora you can enter these commands instead:

| ``sudo yum install openbabel``
| ``sudo yum install python-openbabel``
| ``sudo yum install python-bitarray``
| ``sudo yum install numpy``

After the requirements fulfilled

After all of the dependencies installed, you can install HIPPOS by opening
the terminal and enter the HIPPOS directory and run setup.sh like so:

``./setup.sh``

You will be prompted with a question whether to install HIPPOS in ``[HOME_DIRECTORY]/.hippos``
which is a hidden directory inside your ``HOME_DIRECTORY``. Or would you rather
choose your own directory. If you want to keep the default then just press `enter`
or if you would like to choose your own installation directory you can type `y`
and press enter, then you have to type the installation directory for example
``/home/radifar/apps/hippos`` then press enter.

If HIPPOS installed successfully then a message like 'HIPPOS successfully
installed' should appear. When HIPPOS is already installed and you're running
setup.sh a message like 'HIPPOS is already installed' will appear, and the
installation process will stop and exit.

If you would like to install newer version of HIPPOS and overwrite the old
one then all you need to do is by adding 'force' option to setup.sh like so:

``./setup.sh force``

To test if HIPPOS had been installed, *open new command line window* and type the following:

``hippos``

Then press enter. Please note that it is imperative to open the new command
line window right after the installation to allow the alias in your system get updated.
Next type the following:

``hippos-genref``

Then press enter

If there are no error message then the installation is success.

Installing on Windows 10
------------------------

Windows 10 provides Ubuntu in their `Microsoft Store <https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6>`_. 
After you successfully add Ubuntu on your Windows, the next step is to do the same installation steps on Ubuntu.

