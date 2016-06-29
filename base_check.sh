#!/bin/bash

echo "Installing base dependencies for SSUnique..."

# SOFTWARE:
# FastTree - http://www.microbesonline.org/fasttree/FastTree
# infernal - http://eddylab.org/infernal/infernal-1.1.1.tar.gz 
# ssualign - http://eddylab.org/software/ssu-align/ssu-align-0.1.1.tar.gz
# base R - sudo apt-get install r-base r-base-dev
# HMMER: - http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz

mkdir -p downloads

cd downloads
# Installing FastTree
echo "Installing FastTree"
if hash FastTree 2>/dev/null; then
    echo "FastTree already installed"
else
    echo "Installing from http://www.microbesonline.org/fasttree/FastTree"
    wget "http://www.microbesonline.org/fasttree/FastTree"
    cp FastTree /usr/local/bin/
    chmod a+x /usr/local/bin/FastTree
    rm FastTree
fi

# Installing Infernal
echo "Installing Infernal"
if hash cmalign 2>/dev/null; then
    echo "Infernal already installed"
else
    echo "Installing from http://eddylab.org/infernal/infernal-1.1.1.tar.gz"
    wget "http://eddylab.org/infernal/infernal-1.1.1.tar.gz"
    tar -xf infernal-1.1.1.tar.gz
    cd infernal-1.1.1
    ./configure
    make
    make install
    cd ..
fi

# Installing ssu-align
echo "Installing ssu-align"
if hash ssu-align 2>/dev/null; then
    echo "ssu-align already installed"
else
    echo "Installing from http://eddylab.org/software/ssu-align/ssu-align-0.1.1.tar.gz"
    wget "http://eddylab.org/software/ssu-align/ssu-align-0.1.1.tar.gz"
    tar -xf ssu-align-0.1.1.tar.gz
    cd ssu-align-0.1.1
    ./configure
    make
    make install
    export PATH="$PATH:/usr/local/bin"
    export MANPATH="$MANPATH:/usr/local/share/man"
    export SSUALIGNDIR="/usr/local/share/ssu-align-0.1.1"
    cd ..
fi

# Installing HMMER
echo "Installing HMMER"
if hash hmmbuild 2>/dev/null; then
    echo "HMMER already installled"
else
    echo "Installing from http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz"
    wget "http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz"
    tar -xf hmmer-3.1b2-linux-intel-x86_64.tar.gz
    cd hmmer-3.1b2-linux-intel-x86_64
    ./configure
    make
    make install
    cd ..
fi

# Installing R
echo "Installing base R"
if has R 2>/dev/null; then
    echo "R already installed"
else
    echo "Installing R"
    sudo apt-get install r-base r-base-dev
fi

cd ..

rm -rf downloads
