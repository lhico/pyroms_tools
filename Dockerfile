FROM ubuntu:18.04
MAINTAINER Danilo A. Silva <nilodna@gmail.com>, Dalton K. Sasaki <dalton.sasaki@gmail.com>

# install dependencies for ROMS
RUN apt-get update && apt-get -y install sudo git wget bzip2 curl nano build-essential gfortran && \
    apt-get -y install libnetcdf-dev libnetcdff-dev libhdf5-serial-dev  &&\
    apt-get -y install libkernlib1-gfortran netcdf-bin hdf5-tools &&\
    apt-get install -y apt-utils &&\
    apt-get -y install python3-pip &&\
    apt-get -y install python-setuptools &&\
    apt-get -y install nco
#RUN pip3 install distribute && pip3 install nose
RUN pip3 install nose

# add user lhico with no password, add to sudo group
RUN adduser --disabled-password --gecos '' lhico &&\
    adduser lhico sudo                           &&\
    echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

USER lhico

# set a working directory to open when run the container
WORKDIR /home/lhico/

# Install miniconda to /miniconda
RUN chmod a+rwx /home/lhico/ && cd /home/lhico/ &&\
    curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
    bash Miniconda3-latest-Linux-x86_64.sh -p /home/lhico/miniconda3 -b &&\
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH=/home/lhico/miniconda3/bin:${PATH}

# clonning repository
RUN cd /home/lhico/ && git clone -b  python3 https://github.com/ESMG/pyroms pyroms3 &&\
    cd pyroms3 &&\
    conda update -y conda

RUN conda install -y python=3.7 \
                     numpy \
                     netcdf4 \
                     matplotlib=3.2 \
                     basemap \
                     scipy \
                     basemap-data-hires \
                     ipython \
                     xarray

RUN conda install -y -c conda-forge \
                      pygridgen

RUN conda install -y --channel https://conda.anaconda.org/conda-forge esmf

# setting environment variables
# PROJ_LIB is needed to run basemap on python
# PYTHONPATH to import packages installed on site-packages
# CURDIR is only to make things easier during the installation process
ENV PROJ_LIB="/home/lhico/miniconda3/share/proj" \
PYTHONPATH=/home/lhico/miniconda3/lib/python3.7/site-packages/:${PYTHONPATH} \
CURDIR=/home/lhico/pyroms3 \
SITE_PACKAGES=/home/lhico/miniconda3/lib/python3.7/site-packages

# installing PyROMS and other packages, such as bathy_smoother and pyroms_toolbox
# and dependencias (csa, nn, gridutils, gridgen, natgrid)
RUN cd $CURDIR/pyroms && python setup.py build && python setup.py install &&\
    cp -r $CURDIR/pyroms/pyroms/* /home/lhico/miniconda3/lib/python3.7/site-packages/pyroms/ &&\
    cd $CURDIR/pyroms_toolbox && python setup.py build && python setup.py install &&\
    cd $CURDIR/bathy_smoother && python setup.py build && python setup.py install

#-- changing to root privileges to avoid nosuid related error --#
USER root

RUN cd $CURDIR/pyroms/external/nn && ./configure && sudo make install &&\
    cd $CURDIR/pyroms/external/csa && ./configure && sudo make install &&\
    cd $CURDIR/pyroms/external/gridutils && ./configure && sudo make && sudo make install &&\
    cd $CURDIR/pyroms/external/gridgen && ./configure && sudo make &&\
        sudo make lib && sudo make shlib && sudo make install

ENV PREFIX=/home/lhico/miniconda3
RUN cd $CURDIR/pyroms/external/scrip/source &&\
    make DEVELOP=1 PREFIX=$PREFIX install &&\
    cp scrip*so ../../../

RUN cd $CURDIR &&\
cp -v  ./pyroms/external/scrip/source/scrip.*.so ${SITE_PACKAGES}/pyroms &&\
cp -v ./pyroms/build/lib.linux-x86_64-3.7/pyroms/_iso.cpython-37m-x86_64-linux-gnu.so ${SITE_PACKAGES} &&\
cp -v ./pyroms/build/lib.linux-x86_64-3.7/pyroms/_obs_interp.cpython-37m-x86_64-linux-gnu.so ${SITE_PACKAGES} &&\
cp -v ./pyroms/build/lib.linux-x86_64-3.7/pyroms/_interp.cpython-37m-x86_64-linux-gnu.so ${SITE_PACKAGES} &&\
cp -v ./pyroms/build/lib.linux-x86_64-3.7/pyroms/_remapping.cpython-37m-x86_64-linux-gnu.so ${SITE_PACKAGES} &&\
cp -v ./pyroms_toolbox/build/lib.linux-x86_64-3.7/pyroms_toolbox/_average.cpython-37m-x86_64-linux-gnu.so ${SITE_PACKAGES} &&\
cp -v ./pyroms/build/lib.linux-x86_64-3.7/pyroms/_remapping_fast.cpython-37m-x86_64-linux-gnu.so ${SITE_PACKAGES} &&\
cp -v ./pyroms_toolbox/build/lib.linux-x86_64-3.7/pyroms_toolbox/creep.cpython-37m-x86_64-linux-gnu.so ${SITE_PACKAGES} &&\
cp -v ./pyroms/build/lib.linux-x86_64-3.7/pyroms/_remapping_fast_weighted.cpython-37m-x86_64-linux-gnu.so ${SITE_PACKAGES}
#RUN    cp -v ./bathy_smoother/build/lib.linux-x86_64-3.7/bathy_smoother/lpsolve55.cpython-37m-x86_64-linux-gnu.so ${SITE_PACKAGES}

# copying *.so installed on /usr/local/lib and installing natgrid
RUN cd ${SITE_PACKAGES}/pyroms &&\
    sudo chmod ugo+x /usr/local/lib/libgridgen.so &&\
    sudo chmod ugo+x /usr/local/lib/libgu.so &&\
    sudo ln -s /usr/local/lib/libgridgen.so . &&\
    sudo ln -s /usr/local/lib/libgu.so . &&\
    sudo ln -s /usr/local/lib/scrip*.so scrip.so &&\
    cd ${HOME} && git clone https://github.com/matplotlib/natgrid.git && cd natgrid/ &&\
    python setup.py install

# lpsolve55 must be installed after bathy_smoother
RUN conda install -y -c conda-forge lpsolve55

ENV PYROMS_GRIDID_FILE=/home/lhico/data/gridid.txt

#-- assigning to use the lhico user when entering the container--#
USER lhico

##------------ CHANGES BELOW --------------###
# setting project's name
#ENV project_name=data

# creating projects folder, on the the $HOME (guest) folder
#RUN mkdir ${HOME}/projects/${project_name}

# copying project's files
# COPY scrips /home/lhico/projects/$scripts
