#!/bin/bash

# If argument present do this

if [ $1 > 0 ]  ; then
    if  [ $1 = 'uninstall' ]; then
        INSTALLDIR=$(cat ~/.bashrc | grep hippos.py | cut -d"'" -f 2 | sed "s/hippos.py//")
        if [ $INSTALLDIR ] ; then
            rm -r $INSTALLDIR
            gawk '$0 !~/alias hippos/ { print $0 }' $HOME/.bashrc > .bashrc.tmp
            cat .bashrc.tmp > $HOME/.bashrc
            echo "HIPPOS has been successfully uninstalled"
        else
            echo "HIPPOS never been installed before"
        fi
    elif [ $1 = 'force' ]; then
        INSTALLDIR=$(cat ~/.bashrc | grep hippos.py | cut -d"'" -f 2 | sed "s/hippos.py//")
        if [ $INSTALLDIR ]; then
            rm -r $INSTALLDIR
        else
            echo "alias hippos='$INSTALLDIR/hippos.py'" >> $HOME/.bashrc
            echo "alias hippos-genref='$INSTALLDIR/hippos_genref.py'" >> $HOME/.bashrc
        fi
        mkdir -p $INSTALLDIR && rmdir $INSTALLDIR && cp -r pyplif_hippos $INSTALLDIR
        chmod +x $INSTALLDIR/hippos.py
        chmod +x $INSTALLDIR/hippos_genref.py
        # alias hippos='$HOME/.hippos/hippos.py'
        # alias hippos-genref='$HOME/.hippos/hippos_genref.py'
        echo "HIPPOS successfully installed"
    else
        echo "Wrong Argument"
    fi
# When no argument present do this
else  
    if  cat "$HOME/.bashrc" | grep "hippos" > /dev/null ; then
        echo "HIPPOS installation failed because HIPPOS is already installed"
        echo "Use option 'force' to overwrite the old installation:"
        echo "./setup.sh force"
        exit
    else
        read -p "HIPPOS will be installed in $HOME/.hippos, do you want to install it in another directory?\
`echo $'\n> '`[y]es/[n]o, default: n `echo $'\n> '`" NODEFAULT
        if [ "$NODEFAULT" = '' ] || [ "$NODEFAULT" = 'n' ]; then
            if [ -e $HOME/.hippos ]; then
                rm -r $HOME/.hippos
            fi
            cp -r pyplif_hippos $HOME/.hippos
            chmod +x $HOME/.hippos/hippos.py
            chmod +x $HOME/.hippos/hippos_genref.py
            echo "alias hippos='$HOME/.hippos/hippos.py'" >> $HOME/.bashrc
            echo "alias hippos-genref='$HOME/.hippos/hippos_genref.py'" >> $HOME/.bashrc
            # alias hippos='$HOME/.hippos/hippos.py'
            # alias hippos-genref='$HOME/.hippos/hippos_genref.py'
            echo "HIPPOS successfully installed"
        elif [ "$NODEFAULT" = 'y' ]; then
            read -p "Enter the installation directory: " INSTALLDIR
            if [ "$INSTALLDIR" = '' ]; then
                echo "You are entering empty directory"
                echo "Exiting now"
                exit
            elif [ -e $INSTALLDIR ]; then
                echo "Installation directory already exist"
                echo "Aborting installation"
                exit
            else 
                mkdir -p $INSTALLDIR && rmdir $INSTALLDIR && cp -r pyplif_hippos $INSTALLDIR
                chmod +x $INSTALLDIR/hippos.py
                chmod +x $INSTALLDIR/hippos_genref.py
                echo "alias hippos='$INSTALLDIR/hippos.py'" >> $HOME/.bashrc
                echo "alias hippos-genref='$INSTALLDIR/hippos_genref.py'" >> $HOME/.bashrc
                echo "HIPPOS successfully installed"
                exit
            fi
        else
            echo "Wrong input, please enter 'y' or 'n'"
            echo "exiting now"
            exit
        fi
    fi
fi
