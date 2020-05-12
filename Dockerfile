FROM python:2.7.16-slim-stretch

ENV iraf /iraf/iraf/
ENV IRAFARCH linux64
ENV LCOSNDIR /supernova
ENV LCOSNPIPE $LCOSNDIR/github/lcogtsnpipe
ENV DISPLAY host.docker.internal:0

RUN apt-get update \
        && apt -y install gcc make flex git \
        && apt -y install libcurl4-openssl-dev libexpat-dev libreadline-dev \
        && apt-get autoclean \
        && rm -rf /var/lib/apt/lists/*

RUN mkdir -p $iraf \
        && cd /iraf \
        && git clone https://github.com/iraf-community/iraf.git \
        && cd $iraf \
        && git checkout 567961f \
        && ./install < /dev/null \
        && make linux64 \
        && make sysgen

RUN apt-get update \
        && apt-get -y install libx11-dev libcfitsio-bin wget x11-apps libtk8.6 sextractor \
        mysql-client libmariadbclient-dev openssh-client wcstools libxml2 vim libssl1.0.2 zip \
        && apt-get autoclean \
        && rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/sextractor /usr/bin/sex

RUN pip install numpy>=1.12 astropy matplotlib pyraf mysql-python scipy astroquery statsmodels==0.10 cython

RUN pip install git+git://github.com/kbarbary/sep.git@master git+git://github.com/dguevel/PyZOGY.git && rm -rf ~/.cache/pip

RUN wget http://ds9.si.edu/download/debian9/ds9.debian9.8.0.1.tar.gz \
        && tar -xzvf ds9.debian9.8.0.1.tar.gz -C /usr/local/bin \
        && rm -rf ds9.debian9.8.0.1.tar.gz

RUN mkdir -p /home/supernova/iraf && /usr/sbin/groupadd -g 20000 "domainusers" \
        && /usr/sbin/useradd -g 20000 -d /home/supernova -M -N -u 10197 supernova \
        && chown -R supernova:domainusers /home/supernova \
        && mkdir -p $LCOSNDIR

RUN wget http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz \
        && tar -xzvf cdsclient.tar.gz -C /usr/src && rm cdsclient.tar.gz \
        && cd /usr/src/cdsclient-* && ./configure && make && make install

COPY ./configure $LCOSNDIR

COPY . $LCOSNPIPE

RUN chown -R supernova:domainusers $LCOSNDIR /usr/local

USER supernova

WORKDIR $LCOSNPIPE/trunk

RUN mv src/lsc/standard/cat/apass/*.cat $LCOSNDIR/data/apass/ 2> /dev/null || true

RUN mv src/lsc/standard/cat/sloan/*.cat $LCOSNDIR/data/sloan/ 2> /dev/null || true

RUN python setup.py build -f && python setup.py install -f

RUN ln -s $LCOSNDIR/data/apass/ $LCOSNDIR/data/sloan/ /usr/local/lib/python2.7/site-packages/lsc/standard/cat/

WORKDIR /home/supernova/iraf

RUN mkiraf --term=xgterm -i

WORKDIR /home/supernova

ENTRYPOINT /bin/bash
