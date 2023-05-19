#!/bin/sh

if [ $1 > 0 ] ; then
    if [ $1 = 'uninstall' ]; then
        cp ~/.bashrc .bashrc-copy
        grep -Ev hippos .bashrc-copy > ~/.bashrc
        rm .bashrc-copy
        echo "PyPLIF HIPPOS has been successfully uninstalled"
        echo "PyPLIF HIPPOS will not be able to run the next time"
        echo "you open new terminal"
    else
        echo "$1 command is not recognized"
    fi
else
    sourcedir=`pwd`
    cp ~/.bashrc ~/.bashrc-ori
    grep -Ev hippos ~/.bashrc-ori > ~/.bashrc
    echo "alias hippos='$sourcedir/src/pyplif_hippos/hippos.py'" >> ~/.bashrc
    echo "alias hippos-genref='$sourcedir/src/pyplif_hippos/hippos_genref.py'" >> ~/.bashrc 
    chmod u+x $sourcedir/src/pyplif_hippos/hippos.py
    chmod u+x $sourcedir/src/pyplif_hippos/hippos_genref.py
    echo "PyPLIF HIPPOS is ready for service"
    echo "To use PyPLIF HIPPOS, you have to open new terminal or"
    echo "run 'source ~/.bashrc' in this terminal."
fi
