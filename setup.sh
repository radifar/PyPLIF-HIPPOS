#!/bin/bash

# A script to copy the hippos to hidden folder inside
# HOME folder then add a 'command-line shortcut'
# to HIPPOS inside that hidden folder.
if [ $1 > 0 ]  ; then
    if  [ $1 = 'uninstall' ]; then
        if [ -e $HOME/.hippos ]; then
            rm -r $HOME/.hippos
            gawk '$0 !~/alias hippos/ { print $0 }' $HOME/.bashrc > .bashrc.tmp
            cat .bashrc.tmp > $HOME/.bashrc
            echo "HIPPOS has been successfully uninstalled"
        else
            echo "HIPPOS never been installed before"
        fi
    elif [ $1 = 'force' ]; then
        if [ -e $HOME/.hippos ]; then
            rm -r $HOME/.hippos
        else
            echo "alias hippos='$HOME/.hippos/hippos.py'" >> $HOME/.bashrc
            echo "alias hippos-genref='$HOME/.hippos/hippos-genref.py'" >> $HOME/.bashrc
        fi
        cp -r hippos $HOME/.hippos
        chmod +x $HOME/.hippos/hippos.py
        chmod +x $HOME/.hippos/hippos-genref.py
        alias hippos='$HOME/.hippos/hippos.py'
        alias hippos-genref='$HOME/.hippos/hippos-genref.py'
        echo "HIPPOS successfully installed"
    else
        echo "Wrong Argument"
    fi
else  
    if  cat "$HOME/.bashrc" | grep "hippos" > /dev/null ; then
        echo "HIPPOS installation failed because HIPPOS is already installed"
        echo "Use option 'force' to overwrite the old installation:"
        echo "./setup.sh force"
        exit
    else
        if [ -e $HOME/.hippos ]; then
            rm -r $HOME/.hippos
        fi
        cp -r hippos $HOME/.hippos
        chmod +x $HOME/.hippos/hippos.py
        chmod +x $HOME/.hippos/hippos-genref.py
        echo "alias hippos='$HOME/.hippos/hippos.py'" >> $HOME/.bashrc
        echo "alias hippos-genref='$HOME/.hippos/hippos-genref.py'" >> $HOME/.bashrc
        alias hippos='$HOME/.hippos/hippos.py'
        alias hippos-genref='$HOME/.hippos/hippos-genref.py'
        echo "HIPPOS successfully installed"
    fi
fi
