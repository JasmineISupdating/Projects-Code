############Run gistic2.0 in virtual machine##############
wget -c ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTIC_2_0_23.tar.gz
mkdir GISTIC2
mv GISTIC_2_0_23.tar.gz GISTIC2/ && cd GISTIC2/
tar -zxvf GISTIC_2_0_23.tar.gz
# $ ls
# examplefiles                          LICENSE.txt
# example_results                       MATLAB_Compiler_Runtime
# gistic2                               MCR_Installer
# GISTIC_2_0_23.tar.gz                  README.txt
# GISTICDocumentation_standalone_files  refgenefiles
# GISTICDocumentation_standalone.htm    run_gistic_example
# gp_gistic2_from_seg                   source
# INSTALL.txt

## install Matlab environment
cd MCR_Installer/
unzip MCRInstaller.zip 
./install -mode silent -agreeToLicense yes -destinationFolder /home/jasmine97/GISTIC2/MATLAB_Compiler_Runtime/

## configure environment variables
export LD_LIBRARY_PATH=/home/jasmine97/GISTIC2/MATLAB_Compiler_Runtime/v83/runtime/glnxa64:/home/jasmine97/GISTIC2/MATLAB_Compiler_Runtime/v83/bin/glnxa64:/home/jasmine97/GISTIC2/MATLAB_Compiler_Runtime/v83/sys/os/glnxa64:
export XAPPLRESDIR=/home/jasmine97/GISTIC2/MATLAB_Compiler_Runtime/v83/X11/app-defaults

## error solution
sudo find / -name "libstdc++.so*"
cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /home/jasmine97/GISTIC2/MATLAB_Compiler_Runtime/v83/sys/os/glnxa64/

#### run example data ####
cd
cd GISTIC2
./run_gistic_example

#########run my data############
cd
cd GISTIC2
./run_gistic_st

### write the script "run_gistic_st" ###
#!/bin/sh
## run example GISTIC analysis

## output directory
echo --- creating output directory ---
basedir=`pwd`/st2_results
mkdir -p $basedir 

echo --- running GISTIC ---
## input file definitions
segfile=`pwd`/stfiles/St2_segemnt_file.txt
markersfile=`pwd`/stfiles/Marker_file.txt
refgenefile=`pwd`/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat

## call script that sets MCR environment and calls GISTIC executable 
./gistic2 -b $basedir -seg $segfile -mk $markersfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme
