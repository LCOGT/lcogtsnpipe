FROM continuumio/miniconda3:25.3.1-1
RUN conda create -y -n lcogtsnpipe python=2
SHELL ["conda", "run", "-n", "lcogtsnpipe", "--live-stream", "/bin/bash", "-c"]

ENV iraf=/iraf/iraf/
# To make this a 32 bit version linux64 -> linux
ENV IRAFARCH=linux64

RUN sed -i 's/http:/https:/g' /etc/apt/sources.list.d/debian.sources 2>/dev/null || true \
        && sed -i 's/http:/https:/g' /etc/apt/sources.list 2>/dev/null || true \
        && apt-get --allow-releaseinfo-change update \
        && apt -y install gcc make flex git gfortran zlib1g-dev bison \
        && apt -y install libcurl4-openssl-dev libexpat-dev libreadline-dev gettext \
        && apt-get autoclean \
        && rm -rf /var/lib/apt/lists/*

RUN mkdir -p $iraf \
        && cd /iraf \
        && git clone https://github.com/iraf-community/iraf.git \
        && cd $iraf \
        && git checkout a12a399 \
        && make \
        && make install

RUN apt-get --allow-releaseinfo-change update \
        && apt-get -y install libx11-dev libcfitsio-bin wget x11-apps libtk8.6 source-extractor procps g++ \
        default-mysql-client libmariadb-dev default-libmysqlclient-dev openssh-client wcstools libxml2 vim zip pkg-config \
        libpng-dev libfreetype6-dev libcfitsio-dev libffi-dev libopenblas-dev libssl-dev libfftw3-dev libopenblas-dev \
        && apt-get autoclean \
        && rm -rf /var/lib/apt/lists/*

RUN apt-get --allow-releaseinfo-change update \
        && apt-get -y install autoconf \
        && apt-get -y install automake \
        && apt-get -y install libtool-bin \
        && cd / \
        && git clone https://github.com/astromatic/swarp.git \
        && cd swarp \
        && ./autogen.sh \
        && ./configure \
        && make \
        && ln -s /swarp/src/swarp /usr/bin/

RUN ln -s /usr/bin/source-extractor /usr/bin/sex

RUN pip install numpy>=1.12

# Running this line to assign the instance reconnect
RUN sed  '/st_mysql_options options;/a unsigned int reconnect;' /usr/include/mysql/mysql.h -i.bkp 

RUN pip install cryptography==2.4.1 astropy matplotlib==2.2.5 pyraf mysql-python scipy astroquery==v0.4 statsmodels==0.10 cython reproject

RUN pip install sep==1.0.3 git+https://github.com/dguevel/PyZOGY.git

RUN wget http://ds9.si.edu/download/debian12x86/ds9.debian12x86.8.6.tar.gz \
        && tar -xzvf ds9.debian12x86.8.6.tar.gz -C /usr/local/bin \
        && rm -rf ds9.debian12x86.8.6.tar.gz

RUN wget http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz \
        && tar -xzvf cdsclient.tar.gz -C /usr/src && rm cdsclient.tar.gz \
        && cd /usr/src/cdsclient-* && ./configure && make && make install

RUN cd / \
        && git clone https://github.com/acbecker/hotpants.git \
        && cd hotpants \
        && sed -i 's/^COPTS = .*/& -fcommon/' Makefile \
        && make \
        && ln -s /hotpants/hotpants /usr/bin/

ENV LCOSNPIPE=/lcogtsnpipe

RUN mkdir -p /home/supernova/iraf && /usr/sbin/groupadd -g 20000 "domainusers" \
        && /usr/sbin/useradd -g 20000 -d /home/supernova -M -N -u 10197 supernova \
        && chown -R supernova:domainusers /home/supernova \
        && mkdir -p $LCOSNPIPE

RUN chown -R supernova:domainusers $LCOSNPIPE /opt/conda/envs/lcogtsnpipe/

USER supernova

ENV USER=supernova

COPY --chown=supernova:domainusers . $LCOSNPIPE

WORKDIR $LCOSNPIPE/trunk

RUN python setup.py build -f && python setup.py install -f

WORKDIR /home/supernova/iraf

RUN mkiraf -i -d -c

WORKDIR /home/supernova

COPY bashrc /home/supernova/.bashrc

ENTRYPOINT ["/bin/bash"]
