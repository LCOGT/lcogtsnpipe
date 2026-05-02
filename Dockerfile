FROM continuumio/miniconda3:26.3.2-2
RUN conda create -y -n lcogtsnpipe python=3.11 pip
RUN conda install -y -n lcogtsnpipe pyqt=6.11.0
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
        && git checkout d980a65 \
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

RUN python -m pip install --upgrade pip==25.3 setuptools==80.9.0 wheel==0.45.1

RUN python -m pip install numpy==2.4.4

# Running this line to assign the instance reconnect
RUN sed  '/st_mysql_options options;/a unsigned int reconnect;' /usr/include/mysql/mysql.h -i.bkp 

RUN python -m pip install cryptography==47.0.0 astropy==7.2.0 matplotlib==3.10.8 pyraf==2.2.4 mysqlclient==2.2.8 scipy==1.17.1 astroquery==0.4.11 statsmodels==0.14.6 cython==3.2.4 reproject==0.19.0

RUN python -m pip install sep==1.4.1

RUN git clone https://github.com/dguevel/PyZOGY.git /tmp/PyZOGY \
        && sed -i 's/PyZOGY.__main__ : main/PyZOGY.__main__:main/' /tmp/PyZOGY/setup.py \
        && python -m pip install /tmp/PyZOGY \
        && rm -rf /tmp/PyZOGY

RUN case "$(uname -m)" in \
        aarch64) DS9_PKG="debian12arm64" ;; \
        x86_64)  DS9_PKG="debian12x86" ;; \
        *) echo "Unsupported architecture: $(uname -m)" && exit 1 ;; \
    esac \
        && wget "http://ds9.si.edu/download/${DS9_PKG}/ds9.${DS9_PKG}.8.6.tar.gz" \
        && tar -xzvf "ds9.${DS9_PKG}.8.6.tar.gz" -C /usr/local/bin \
        && rm -f "ds9.${DS9_PKG}.8.6.tar.gz"

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