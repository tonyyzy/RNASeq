if [ ! -d "./software" ]; then
    mkdir ./software
    # install STAR
    wget https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz
    tar -xzf 2.6.0a.tar.gz
    cp ./STAR-2.6.0a/bin/Linux_x86_64/STAR ./software
    # install samtools
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar -vxjf samtools-1.9.tar.bz2
    cd samtools-1.9
    make
    make prefix="$PWD" install
    cp ./bin/samtools ../software
    cd ..
    # install stringtie
    git clone -n https://github.com/gpertea/stringtie
    cd stringtie
    git checkout 2584c142d7e636f3a9ce3f0546511a0974dc0609
    make release
    cp stringtie ../software
    cd ..
    # clean up
    echo "Clean up"
    rm ./2.6.0a.tar.gz ./samtools-1.9.tar.bz2
    rm -rf ./STAR-2.6.0a ./samtools-1.9 ./stringtie
fi
