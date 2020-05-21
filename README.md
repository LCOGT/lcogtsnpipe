#  lcogtsnpipe
This is the pipeline that ingests and reduces new data from the lcogt key project

# Table of Contents
- [ Pipeline Documentation](#pipeline-documentation)
- [Installing the Pipeline and Database](#installing-the-pipeline-and-database)
- [Testing your installation](#testing-your-installation)
- [Appendix A: Expected output from show tables](#appendix-a-expected-output-from-show-tables)
- [Appendix B: Installing 64 bit IRAF on Catalina ](#appendix-b-installing-64-bit-iraf-on-catalina)
- [Appendix C: Other packages you may need to install](#appendix-c-other-packages-you-may-need-to-install)
- [Appendix D: Installing Source Extractor](#appendix-d-installing-source-extractor)
- [Appendix E: Using the Pipeline with Docker](#appendix-e-using-the-pipeline-with-docker)

# Pipeline Documentation:
1. Image Subtraction: https://www.authorea.com/users/75900/articles/96044-image-subtraction-with-lcogtsnpipe
2. [Manual (beta)](manual.md) 

# Installing the Pipeline and Database
1. Install msql  
    1. Install MySQL server from: https://dev.mysql.com/downloads/mysql/  
    **Note:** you don’t need an Oracle account to sign up, there is a just start download link at the bottom of the page. Links from here were very useful: https://dev.mysql.com/doc/refman/5.6/en/installing.html
    2. Additional notes for Mac users: I used the DMG for OS X to install. This installs things in /usr/local. You can see if mysql is installed correctly by going to System Preferences. There should be mysql icon. Clicking on this you can start a mysql server
2. Set mysql to only throw a warning (rather than an error) when inserting or updating empty values into fields without default values:

    ```
    mysql -uroot -p

    mysql> set session sql_mode = 'NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_AUTO_CREATE_USER,NO_ENGINE_SUBSTITUTION';

    mysql> set global sql_mode = 'NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_AUTO_CREATE_USER,NO_ENGINE_SUBSTITUTION';

    mysql> exit 
    ```

3. Create database called Supernova
    1. download supernova database: https://www.dropbox.com/s/b28k5dkbs3i48nd/supernova_2020-01-31.sql?dl=0
    2. Connect to mysql as root: 
        ```
        mysql -u root -p
        ```
    3. Create the supernova database:
        ```
        mysql> CREATE DATABASE supernova;
        ```
        Expected output: Query OK, 1 row affected (0.01 sec)
    4. Create user supernova and grant privileges
        ```
        mysql> GRANT ALL PRIVILEGES ON *.* TO 'supernova'@'localhost' IDENTIFIED BY 'supernova';

        mysql> exit
        ```
        Expected output: Query OK, 0 rows affected, 1 warning (0.01 sec)
    5. Create the supernova database and choose a password:
        ```
        mysql -u root -p *your-password* < *database-including-path*
        
        e.g. mysql -u root -p zenith < /Users/valenti/Desktop/supernova_2020-01-31.sql
        ```
    6. Check that the database structure has been loaded
        * Start mysql connecting to the supernova database as the supernova user and entering the password you just created
            ```
            mysql -u supernova -D supernova -p
            ```
            Enter password: \*your-password*
        * Check that the list of tables matches the list of tables matches the list in Appendix A
            ```
            mysql> show tables;
            
            mysql> exit
            ```
4. Prepare your pipeline environment using Anaconda
    1. If you do not already have it, download anaconda for Python 3 (https://www.anaconda.com/distribution/)
    2. Add the astroconda channel to your conda search path:
        ```
        conda config --add channels http://ssb.stsci.edu/astroconda
        ```
    3. Create a conda environment called lcogtsnpipe for the pipeline that uses Python 2.7
        ```
        conda create -n lcogtsnpipe python=2.7
        ```
    4. Switch to your new lcogtsnpipe environment
        ```
        conda activate lcogtsnpipe
        ``` 
    5. install astroconda (for details https://astroconda.readthedocs.io/en/latest/installation.html)
        * If you are running mac OS <10.15 (catalina) then install with iraf:
            ```
            conda install iraf-all pyraf-all stsci
            ```
        * If you are running macOS >= 10.15 (catalina) then install without iraf and see Appendix B for 64bit IRAF installation instructions. Then:
            ```
            conda install stsci
            ```
    6. Install astroquery
        ```
        conda install astroquery
        ```
    7. Install mysqldb
        ```
        conda install MySQL-python
        ```
5. Install the pipeline
    1. Install git if you do not already have it (https://ucdavis.github.io/DS4S/#setup)
    2. Download the pipeline:
        ```
        git clone https://github.com/svalenti/lcogtsnpipe
        ```
    3. Switch to the correct version of the pipeline
        ```
        cd lcogtsnpipe
        git checkout exportpipe
        ```
    4. Install the pipeline:
        ```
        cd trunk
        python setup.py install
        ```
6. Set up directory for data and configuration file
    ```
    mkdir <your directory name>
    e.g. mkdir /Users/valenti/lco
    ```
7. Set LCOSNDIR environment variable either on the command line for a single session or in your .bashrc or .bash_profile file (this assumes you are using a bash shell) that points to where your data will live
    ```
    export LCOSNDIR=<your directory name>
    e.g. export LCOSNDIR='/Users/valenti/lco'
    ```
8. Set LCOSNPIPE environment variable either on the command line for a single session or in your .bashrc or .bash_profile file (this assumes you are using a bash shell) that points to your lcogtsnpipe git repository
    ```
    export LCOSNPIPE=<your directory name>
    e.g. export LCOSNPIPE='/Users/valenti/lco/lcogtsnpipe'
    ```
9. Create configuration file named configure in LCOSNDIR with the following lines:
    ```
    hostname           127.0.0.1
    database           supernova
    mysqluser          supernova
    mysqlpasswd        *your-password*
    proposal           ['']
    users              ['']
    triggerpass        ''
    extraobject        ['']
    skipobjects        ['auto_focus']
    proposalingestion  ['']
    ptfhost            ''
    ptfdatabase        ''
    ptfpasswd          '''
    ptfuser            ''
    ```
    **Don't forget** to fill in your password on the third line
10. Install vizquery: https://vizier.u-strasbg.fr/vizier/doc/cdsclient.html  
    **Note:** I had to sudo make install on my mac

# Testing your installation
1. Download test data: 
https://www.dropbox.com/s/o3ls0zcqd64f49f/snexdata_2020-01-30_07%3A37%3A02.455163_23_6139.tar.gz?dl=0
2. Add your test data to the Supernova database:
    ```
    ingesttar.py -f snexdata_2020-01-30_07-37-02.455163_23_6139.tar
    ```
    **Note:** if you run this more than once, it may delete your entry. Check that the photlco table has rows. If not, run a third time.
3. Create a standard star catalog for apass and sdss catalog for AT2020oi field
    ```
    comparecatalogs.py
    ```
    **Note:** if you run this more than once, you may need to use the -F option to force the script to look again for a catalog
4. Run cosmic ray rejection on the images:
    ```
    lscloop.py -e 20200101-20200129 -s cosmic
    ```
5. Generate a PSF model for images:
    ```
    lscloop.py -e 20200101-20200130 -f apass --catalog=<path-to-your-conda>/envs/lcogtsnpipe/lib/python2.7/site-packages/lsc/standard/cat/apass/AT2020oi_apass.cat -s psf
    ```
    Where \<path-to-your-conda> is the path to your anaconda installation
6. Calculate instrumental magnitudes with PSF photometry, displaying output in DS9
    1. Open DS9
    ```
    ds9& 
    lscloop.py -e 20200101-20200129 -s psfmag --show
    ```
7. Find the zeropoint of the image:
    ```
    lscloop.py -e 20200101-20200130 -f apass --catalog=<path-to-your-conda>/envs/lcogtsnpipe/lib/python2.7/site-packages/lsc/standard/cat/apass/AT2020oi_apass.cat -s zcat
    ```
    Where \<path-to-your-conda> is the path to your anaconda installation
8. Find the apparent magnitude using zeropoint
    ```
    lscloop.py -e 20200101-20200129 -s mag  -F -f apass --type fit
    ```
    **Note:** this step may give you a warning at the end: Error 1364: Field 'targetid' doesn't have a default value
9. Check the magnitudes you just calculated by plotting them
    ```
    lscloop.py -e 20200101-20200129 -s getmag --type mag --show
    ```

# Appendix A: Expected output from show tables
```
mysql> show tables;
+---------------------+
| Tables_in_supernova |
+---------------------+
| aseatide            |
| atels               |
| classifications     |
| contention          |
| datarequests        |
| dbsyncs             |
| eseatide            |
| favorites           |
| glade               |
| glade_2             |
| groups              |
| headerdefaults      |
| hitslogger          |
| iaunames            |
| instruments         |
| interests           |
| lvc_galaxies        |
| lvc_triggers        |
| notes               |
| obslog              |
| obsrequests         |
| obsrequests_tags    |
| papers              |
| permissionlog       |
| phot_deprecated     |
| photlco             |
| photlcoraw          |
| photpairing         |
| programs            |
| psns                |
| reference_status    |
| scheduling_run      |
| schedulinglog       |
| spec                |
| speclcoguider       |
| speclcoraw          |
| tags                |
| targetnames         |
| targets             |
| telescopes          |
| timecharged         |
| useractionlog       |
| userrequests        |
| users               |
| voevent_amon        |
| voevent_lvc         |
| ztf_alerts          |
| ztf_targets         |
+---------------------+
48 rows in set (0.00 sec)	
```

# Appendix B: Installing 64 bit IRAF on Catalina 
Liberally adopted from: https://iraf-community.github.io/install
1. If you haven’t already, install Xcode compilers
    ```
    xcode-select —install
    ```
2. Set iraf environment variables (I don’t know if this is strictly necessary):
    ```
    export TERM='xterm'
    export IRAFARCH='macintel' 
    export OS_VERS='catalina'
    ```
3. Disable System Integrity Protection 
    1. Restart computer in recovery mode by holding down <code>command+R</code> during restart
    2. From the Utilities menu select terminal
        ```
        csrutil disable
        ```
4. Restart your computer
5. Disable read-only access to your root directory
    ```
    sudo mount -uw /
    ```
6. Create iraf directory in /
    ```
    sudo mkdir /iraf
    ```
7. Download IRAF
    ```
    sudo git clone https://github.com/iraf-community/iraf.git
    ```
8. Install IRAF
    1. Install IRAF
        ```
        cd iraf
        sudo git checkout 567961f
        sudo ./install
        ```
        set term to xterm
    2. Add scripts install bin to path (this will be printed to the screen during the install process. For me it was ~/.iraf/bin
        ```
        e.g. export PATH=‘/Users/bostroem/.iraf/bin':$PATH
        ```
    3. build iraf
        ```
        sudo make macintel
        sudo make sysgen 2>&1 | tee build.log
        ```
    4. Edit the last line of <code>~/.iraf/setup.sh</code> replacing <code>xgterm -e</code> with <code>xterm</code>
    5. Install Pyraf
        Make sure you are in your lcogtsnpipe environment (conda activate lcogtsnpipe)
        conda install pyraf

# Appendix C: Other packages you may need to install
* Source Extractor:
    * http://www.astromatic.net/software/sextractor
    **Note:** This was non-trivial on my Mac. See Appendix D for directions
* For difference image photometry:
    * *hotpants*  https://github.com/acbecker/hotpants  
        **Note:** There is a pull request submitted to this repo to fix issues with mac installation. I ended up cloning https://github.com/exowanderer/hotpants.git, checking out the remote 
    * *PyZOGY* https://github.com/dguevel/PyZOGY

# Appendix D: Installing Source Extractor
In the directions below, <your username> is used repeatedly as a place holder for youractual username (e.g. bostroem)
1. Install *FFTW* from source:  
    * Download: http://www.fftw.org/download.html  
    * Install:  
    ```
    ./configure --prefix=/Users/<your username>
    make
    make install
    ```
2. Install *OpenBLAS* from source:
    * Download tar.gz file from link at top of page: https://www.openblas.net
    * Install:
    ```
    make PREFIX=/Users/<your username> install
    ```
3. Install *Source Extractor*:
    * Download: http://www.astromatic.net/software/sextractor
    * Install:
    ```
    sh autogen.sh
    ./configure  --with-fftw-incdir=/Users/<your username>/include --with-fftw-libdir=/Users/<your username>/lib --with-openblas-incdir=/Users/<your username/include --with-openblas-libdir=/Users/<your username>/lib --with-atlas-incdir=/Users/<your username>/include --with-atlas-libdir=/Users/<your username>/lib --enable-openblas
    make -j
    sudo make PREFIX=/Users/<your username> install
    ```
    **Note:** with these direction you do **NOT** need ATLAS
4. Install *vizquery*:
    * Download and follow installation directions: http://cdsarc.u-strasbg.fr/vizier/doc/cdsclient.html

# Appendix E: Using the Pipeline with Docker
You can also run the pipeline from within a Docker container.
In the following instructions, the database server and image directories will live locally on your computer, so they will persist outside the Docker.
The pipeline itself will run inside the Docker container and forward graphics to your local computer.

## Installation
These instructions only need to be run once, when you set up the pipeline.

   1. Install [MySQL](https://dev.mysql.com/downloads/mysql/). You'll set a root password that you'll need to remember later.
       * Make sure to use Legacy Password Encryption. On Mac, this is System Preferences &rarr; MySQL &rarr; Initialize Database &rarr; Use Legacy Password Encryption.
       * Also make sure MySQL is correctly added to your `PATH`.
   2. Install [Docker](https://docs.docker.com/get-docker/).
   3. (MacOS only) Install [XQuartz](https://www.xquartz.org).
   4. (MacOS only) Install [socat](http://www.dest-unreach.org/socat/). If you have [Homebrew](https://brew.sh) installed, you can just run `brew install socat`.
   5. Clone this repository: `git clone https://github.com/svalenti/lcogtsnpipe`
   6. Initialize the database: `mysql -u root -p < lcogtsnpipe/supernova.sql`. You'll need to type the MySQL root password.
   7. Build the Docker image: `docker build -t lcogtsnpipe lcogtsnpipe`
   8. Create directories on your local machine (outside the Docker) where the images and data products will be stored:
      ```
      mkdir /your/data/directory
      cd /your/data/directory
      mkdir lsc fts 0m4 floyds extdata apass sloan
      ```

## Use
Follow these instructions each time you want to use the pipeline.

   1. Make sure the MySQL server and Docker daemon are running.
   2. (MacOS only) Run XQuartz from the Finder.
   3. (MacOS only) Run this hack in the background to get the X11 forwarding to work: `socat TCP-LISTEN:6000,reuseaddr,fork UNIX-CLIENT:\"$DISPLAY\" &`
   4. Run the Docker container: `docker run -it --rm -v /your/data/directory:/supernova/data lcogtsnpipe`.
      Replace `/your/data/directory` with the directory you created above.
   5. When you're done, type `exit` to leave the Docker.
