# ************************************************************
# Sequel Pro SQL dump
# Version 4541
#
# http://www.sequelpro.com/
# https://github.com/sequelpro/sequelpro
#
# Host: 127.0.0.1 (MySQL 5.7.12)
# Database: supernova
# Generation Time: 2020-02-01 00:34:12 +0000
# ************************************************************


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

SET SESSION SQL_MODE = 'NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION';
SET GLOBAL SQL_MODE = 'NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION';
CREATE DATABASE supernova;
CREATE USER supernova;
GRANT ALL PRIVILEGES ON supernova.* TO supernova;
USE supernova;

# Dump of table aseatide
# ------------------------------------------------------------

DROP TABLE IF EXISTS `aseatide`;

CREATE TABLE `aseatide` (
  `GRA` double DEFAULT NULL,
  `GDEC` double DEFAULT NULL,
  `OBJID` bigint(11) DEFAULT NULL,
  `E_BV` double DEFAULT NULL,
  `MAG_FUV` double DEFAULT NULL,
  `MAG_NUV` double DEFAULT NULL,
  `S2N_NUV` double DEFAULT NULL,
  `MAGERR_FUV` double DEFAULT NULL,
  `MAGERR_NUV` double DEFAULT NULL,
  `MAG_D3P0_FUV` double DEFAULT NULL,
  `MAG_D3P0_NUV` double DEFAULT NULL,
  `MAG_D12P0_FUV` double DEFAULT NULL,
  `MAG_D12P0_NUV` double DEFAULT NULL,
  `FLUX20_RADIUS_NUV` double DEFAULT NULL,
  `FLUX50_RADIUS_NUV` double DEFAULT NULL,
  `FLUX90_RADIUS_NUV` double DEFAULT NULL,
  `FLUX20_RADIUS_FUV` double DEFAULT NULL,
  `FLUX50_RADIUS_FUV` double DEFAULT NULL,
  `FLUX90_RADIUS_FUV` double DEFAULT NULL,
  `Z` double DEFAULT NULL,
  `RA0` double DEFAULT NULL,
  `DEC0` double DEFAULT NULL,
  `SPECOBJID` bigint(11) DEFAULT NULL,
  `COLUMN1` bigint(11) DEFAULT NULL,
  `SN_MEDIAN` double DEFAULT NULL,
  `H_ALPHA_EQW` double DEFAULT NULL,
  `H_ALPHA_EQW_ERR` double DEFAULT NULL,
  `LICK_HD_A` double DEFAULT NULL,
  `LICK_HD_A_ERR` double DEFAULT NULL,
  `H_ALPHA_FLUX` double DEFAULT NULL,
  `H_ALPHA_FLUX_ERR` double DEFAULT NULL,
  `D4000_N` double DEFAULT NULL,
  `D4000_N_ERR` double DEFAULT NULL,
  `MODELMAG_U` double DEFAULT NULL,
  `MODELMAG_G` double DEFAULT NULL,
  `MODELMAG_R` double DEFAULT NULL,
  `MODELMAG_I` double DEFAULT NULL,
  `MODELMAG_Z` double DEFAULT NULL,
  `FIBERMAG_U` double DEFAULT NULL,
  `FIBERMAG_G` double DEFAULT NULL,
  `FIBERMAG_R` double DEFAULT NULL,
  `FIBERMAG_I` double DEFAULT NULL,
  `FIBERMAG_Z` double DEFAULT NULL,
  `LGM_TOT_P50` double DEFAULT NULL,
  `SFR_TOT_P50` double DEFAULT NULL,
  `W1MPRO` double DEFAULT NULL,
  `W2MPRO` double DEFAULT NULL,
  `W3MPRO` double DEFAULT NULL,
  `W4MPRO` double DEFAULT NULL,
  `W1RSEMI` double DEFAULT NULL,
  `W2RSEMI` double DEFAULT NULL,
  `W3RSEMI` double DEFAULT NULL,
  `W4RSEMI` double DEFAULT NULL,
  `W1SNR` double DEFAULT NULL,
  `W2SNR` double DEFAULT NULL,
  `W3SNR` double DEFAULT NULL,
  `W4SNR` double DEFAULT NULL,
  `SEARCH_ID` int(11) DEFAULT NULL,
  `MATCHED_ID` bigint(11) DEFAULT NULL,
  `NAME` varchar(255) DEFAULT NULL,
  `num_of_images` int(11) NOT NULL DEFAULT '0'
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table atels
# ------------------------------------------------------------

DROP TABLE IF EXISTS `atels`;

CREATE TABLE `atels` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `atelid` bigint(20) NOT NULL,
  `ateldate` datetime DEFAULT NULL,
  `ra_orig` text COMMENT 'Original RA string in ATel',
  `dec_orig` text COMMENT 'Original dec string in ATel',
  `ra_str` text COMMENT 'Parsed RA string',
  `dec_str` text COMMENT 'Parsed dec string',
  `ra0` double DEFAULT NULL,
  `dec0` double DEFAULT NULL,
  `title` text,
  `keywords` text,
  `authors` text,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table classifications
# ------------------------------------------------------------

DROP TABLE IF EXISTS `classifications`;

CREATE TABLE `classifications` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `name` text NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table contention
# ------------------------------------------------------------

DROP TABLE IF EXISTS `contention`;

CREATE TABLE `contention` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `starttime` datetime DEFAULT NULL,
  `program_id` varchar(20) DEFAULT '',
  `program_title` text,
  `hoursfromnow` float DEFAULT NULL,
  `pressure` float DEFAULT NULL,
  `site` varchar(5) DEFAULT NULL,
  `instrument` varchar(20) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `starttime` (`starttime`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table datarequests
# ------------------------------------------------------------

DROP TABLE IF EXISTS `datarequests`;

CREATE TABLE `datarequests` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(11) DEFAULT NULL,
  `type` char(11) DEFAULT NULL,
  `userid` bigint(11) DEFAULT NULL,
  `formats` char(20) DEFAULT 'fits',
  `ids` text,
  `maxfilesize` float DEFAULT '2000',
  `status` varchar(11) DEFAULT 'new',
  `datecreated` timestamp NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table dbsyncs
# ------------------------------------------------------------

DROP TABLE IF EXISTS `dbsyncs`;

CREATE TABLE `dbsyncs` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `origindb` varchar(100) NOT NULL,
  `origintable` varchar(100) NOT NULL,
  `originid` bigint(20) NOT NULL,
  `destdb` varchar(100) NOT NULL,
  `desttable` varchar(100) NOT NULL,
  `destid` bigint(20) NOT NULL,
  `datecreated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  KEY `datecreated` (`datecreated`),
  KEY `originid` (`originid`),
  KEY `destid` (`destid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Syncing data to/from iPTF and others';



# Dump of table eseatide
# ------------------------------------------------------------

DROP TABLE IF EXISTS `eseatide`;

CREATE TABLE `eseatide` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `ra0` double DEFAULT NULL,
  `dec0` double DEFAULT NULL,
  `redshift` double DEFAULT NULL,
  `source` varchar(30) DEFAULT '',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table favorites
# ------------------------------------------------------------

DROP TABLE IF EXISTS `favorites`;

CREATE TABLE `favorites` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `userid` bigint(20) NOT NULL,
  `favorited` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `unfavorited` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table glade
# ------------------------------------------------------------

DROP TABLE IF EXISTS `glade`;

CREATE TABLE `glade` (
  `ra0` double DEFAULT NULL,
  `dec0` double DEFAULT NULL,
  `dist` double DEFAULT NULL,
  `bmag` double DEFAULT NULL,
  `bt_hyp` double DEFAULT NULL,
  `e_bt_hyp` double DEFAULT NULL,
  `it_hyp` double DEFAULT NULL,
  `e_it_hyp` double DEFAULT NULL,
  `modz_hyp` double DEFAULT NULL,
  `mod0_hyp` double DEFAULT NULL,
  `logd25_hyp` double DEFAULT NULL,
  `e_logd25_hyp` double DEFAULT NULL,
  `logr25_hyp` double DEFAULT NULL,
  `e_logr25_hyp` double DEFAULT NULL,
  `logdc_hyp` double DEFAULT NULL,
  `pa_hyp` double DEFAULT NULL,
  `btc_hyp` double DEFAULT NULL,
  `itc_hyp` double DEFAULT NULL,
  `ubtc_hyp` double DEFAULT NULL,
  `bvtc_hyp` double DEFAULT NULL,
  `jmag_2mass` double DEFAULT NULL,
  `errjmag_2mass` double DEFAULT NULL,
  `hmag_2mass` double DEFAULT NULL,
  `errhmag_2mass` double DEFAULT NULL,
  `kmag_2mass` double DEFAULT NULL,
  `errkmag_2mass` double DEFAULT NULL,
  `ab_ratio_2mass` double DEFAULT NULL,
  `pa_in_kmag_2mass` double DEFAULT NULL,
  `bmag_gwgc` double DEFAULT NULL,
  `majdiam_gwgc` double DEFAULT NULL,
  `errmd_gwgc` double DEFAULT NULL,
  `mindiam_gwgc` double DEFAULT NULL,
  `errmid_gwgc` double DEFAULT NULL,
  `pa_gwgc` double DEFAULT NULL,
  `dist_gwgc` double DEFAULT NULL,
  `errdist_gwgc` double DEFAULT NULL,
  `errbmag_gwgc` double DEFAULT NULL,
  `kmag_2mpz` double DEFAULT NULL,
  `errkmag_2mpz` double DEFAULT NULL,
  `bmag_2mpz` double DEFAULT NULL,
  `errbmag_2mpz` double DEFAULT NULL,
  `errbmag_2_2mpz` double DEFAULT NULL,
  `zspec_2mpz` double DEFAULT NULL,
  `zphot_2mpz` double DEFAULT NULL,
  `errzphot_2mpz` double DEFAULT NULL,
  `errzphot_2_2mpz` double DEFAULT NULL,
  `flag` double DEFAULT NULL,
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`id`),
  KEY `ra0` (`ra0`,`dec0`,`dist`,`bmag`),
  KEY `dec0` (`dec0`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table glade_2
# ------------------------------------------------------------

DROP TABLE IF EXISTS `glade_2`;

CREATE TABLE `glade_2` (
  `pgc` varchar(255) DEFAULT NULL,
  `gwgc_name` varchar(255) DEFAULT NULL,
  `hyperleda_name` varchar(255) DEFAULT NULL,
  `2mass_name` varchar(255) DEFAULT NULL,
  `sdss_dr12_name` varchar(255) DEFAULT NULL,
  `flag1` char(2) DEFAULT NULL,
  `ra0` double DEFAULT NULL,
  `dec0` double DEFAULT NULL,
  `dist` double DEFAULT NULL,
  `dist_err` double DEFAULT NULL,
  `z` double DEFAULT NULL,
  `b` double DEFAULT NULL,
  `b_err` double DEFAULT NULL,
  `j` double DEFAULT NULL,
  `j_err` double DEFAULT NULL,
  `h` double DEFAULT NULL,
  `h_err` double DEFAULT NULL,
  `k` double DEFAULT NULL,
  `k_err` double DEFAULT NULL,
  `flag2` int(11) DEFAULT NULL,
  `flag3` int(11) DEFAULT NULL,
  `b_err_min` double DEFAULT NULL,
  `b_err_max` double DEFAULT NULL,
  `flag4` int(11) DEFAULT NULL,
  `z_err_min` double DEFAULT NULL,
  `z_err_max` double DEFAULT NULL,
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`id`),
  KEY `ra0` (`ra0`),
  KEY `dec0` (`dec0`),
  KEY `dist` (`dist`),
  KEY `b` (`b`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table groups
# ------------------------------------------------------------

DROP TABLE IF EXISTS `groups`;

CREATE TABLE `groups` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(20) NOT NULL,
  `description` text,
  `idcode` bigint(20) NOT NULL,
  `sharing` int(11) DEFAULT NULL,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='User groups';



# Dump of table headerdefaults
# ------------------------------------------------------------

DROP TABLE IF EXISTS `headerdefaults`;

CREATE TABLE `headerdefaults` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `tablename` varchar(20) DEFAULT NULL,
  `columnname` varchar(20) DEFAULT NULL,
  `keywords` text,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table hitslogger
# ------------------------------------------------------------

DROP TABLE IF EXISTS `hitslogger`;

CREATE TABLE `hitslogger` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Logging times that objects are viewed for most popular stars';



# Dump of table iaunames
# ------------------------------------------------------------

DROP TABLE IF EXISTS `iaunames`;

CREATE TABLE `iaunames` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `name` text NOT NULL,
  `ra0` double DEFAULT NULL,
  `dec0` double DEFAULT NULL,
  `type` text,
  `discoverydate` datetime DEFAULT NULL,
  `discoverymag` double DEFAULT NULL,
  `hostname` text,
  `hostra` double DEFAULT NULL,
  `hostdec` double DEFAULT NULL,
  `offsetew` text,
  `offsetns` text,
  `discoveryref` text,
  `positionref` text,
  `discoverer` text,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table instruments
# ------------------------------------------------------------

DROP TABLE IF EXISTS `instruments`;

CREATE TABLE `instruments` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `name` text NOT NULL,
  `type` text,
  `webpage` text,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table interests
# ------------------------------------------------------------

DROP TABLE IF EXISTS `interests`;

CREATE TABLE `interests` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `userid` bigint(20) NOT NULL,
  `interested` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `uninterested` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table lvc_galaxies
# ------------------------------------------------------------

DROP TABLE IF EXISTS `lvc_galaxies`;

CREATE TABLE `lvc_galaxies` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `voeventid` bigint(20) DEFAULT NULL,
  `glade_ver` double DEFAULT NULL,
  `gladeid` bigint(20) DEFAULT NULL,
  `score` double DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table lvc_triggers
# ------------------------------------------------------------

DROP TABLE IF EXISTS `lvc_triggers`;

CREATE TABLE `lvc_triggers` (
  `id` bigint(11) unsigned NOT NULL AUTO_INCREMENT,
  `voeventid` bigint(11) DEFAULT NULL,
  `gladeid` bigint(20) DEFAULT NULL,
  `obsrequestsid` bigint(11) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table notes
# ------------------------------------------------------------

DROP TABLE IF EXISTS `notes`;

CREATE TABLE `notes` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `note` text NOT NULL,
  `tablename` text COMMENT 'Table to which the comment refers',
  `tableid` bigint(20) DEFAULT NULL COMMENT 'ID in the table to which the comment refers',
  `posttime` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `userid` bigint(20) NOT NULL,
  `groupidcode` bigint(20) DEFAULT NULL COMMENT 'Which groups can see this',
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Notes on targets, images, spectra, atels, etc.';



# Dump of table obslog
# ------------------------------------------------------------

DROP TABLE IF EXISTS `obslog`;

CREATE TABLE `obslog` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `user` varchar(20) DEFAULT NULL,
  `targetid` bigint(20) DEFAULT NULL,
  `triggerjd` double DEFAULT NULL,
  `windowstart` double DEFAULT NULL,
  `windowend` double DEFAULT NULL,
  `filters` varchar(30) DEFAULT NULL,
  `exptime` varchar(39) DEFAULT NULL,
  `numexp` varchar(30) DEFAULT NULL,
  `proposal` varchar(30) DEFAULT NULL,
  `site` varchar(10) DEFAULT NULL,
  `instrument` varchar(30) DEFAULT NULL,
  `sky` float DEFAULT '9999',
  `seeing` float DEFAULT '9999',
  `airmass` float DEFAULT '9999',
  `slit` varchar(15) DEFAULT '9999',
  `acqmode` varchar(20) DEFAULT NULL,
  `priority` varchar(50) DEFAULT NULL,
  `ipp` float DEFAULT NULL,
  `reqnumber` int(11) DEFAULT NULL,
  `tracknumber` int(11) DEFAULT NULL,
  `requestsid` bigint(20) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `obslog_targid_idx` (`targetid`),
  KEY `tracknumber` (`tracknumber`),
  KEY `requestsid` (`requestsid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table obsrequests
# ------------------------------------------------------------

DROP TABLE IF EXISTS `obsrequests`;

CREATE TABLE `obsrequests` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `sequencestart` timestamp NULL DEFAULT NULL,
  `sequenceend` timestamp NULL DEFAULT NULL,
  `userstart` varchar(20) NOT NULL DEFAULT '',
  `userend` varchar(20) DEFAULT NULL,
  `cadence` double DEFAULT NULL COMMENT 'Cadence in days',
  `window` double DEFAULT NULL COMMENT 'Window length in days',
  `filters` varchar(100) DEFAULT NULL,
  `exptimes` varchar(100) DEFAULT NULL,
  `expnums` varchar(100) DEFAULT NULL,
  `blocknums` varchar(100) DEFAULT NULL,
  `proposalid` text,
  `ipp` float DEFAULT '1',
  `site` varchar(10) DEFAULT NULL,
  `instrument` varchar(20) DEFAULT NULL,
  `airmass` float DEFAULT '9999' COMMENT 'Requested airmass limit',
  `moondistlimit` float DEFAULT NULL,
  `slit` float DEFAULT '9999',
  `acqradius` double DEFAULT '0',
  `guidermode` varchar(50) DEFAULT NULL,
  `guiderexptime` int(11) DEFAULT '10',
  `priority` varchar(50) DEFAULT NULL,
  `approved` int(11) DEFAULT '0',
  `nextreminder` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `groupidcode` bigint(20) DEFAULT NULL COMMENT 'Which groups can see this',
  `dismissed` int(11) NOT NULL DEFAULT '0',
  `autostop` int(11) NOT NULL DEFAULT '0',
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `combo` (`targetid`,`sequenceend`,`cadence`,`filters`,`exptimes`,`expnums`,`blocknums`,`site`,`instrument`,`airmass`,`slit`,`priority`,`groupidcode`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Requested observations';



# Dump of table obsrequests_tags
# ------------------------------------------------------------

DROP TABLE IF EXISTS `obsrequests_tags`;

CREATE TABLE `obsrequests_tags` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `tagid` bigint(20) NOT NULL,
  `requestsid` bigint(20) NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Connecting observations requests to tags';



# Dump of table papers
# ------------------------------------------------------------

DROP TABLE IF EXISTS `papers`;

CREATE TABLE `papers` (
  `id` bigint(11) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(11) NOT NULL,
  `reference` text,
  `status` text,
  `contents` text,
  `datecreated` timestamp NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table permissionlog
# ------------------------------------------------------------

DROP TABLE IF EXISTS `permissionlog`;

CREATE TABLE `permissionlog` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) DEFAULT NULL,
  `groupname` bigint(20) DEFAULT NULL,
  `jd` float DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table phot_deprecated
# ------------------------------------------------------------

DROP TABLE IF EXISTS `phot_deprecated`;

CREATE TABLE `phot_deprecated` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `photlcoid` bigint(20) DEFAULT NULL,
  `telescopeid` bigint(20) DEFAULT NULL,
  `instrumentid` bigint(20) DEFAULT NULL,
  `exptime` float DEFAULT NULL,
  `mjd` double DEFAULT NULL,
  `mag` double DEFAULT NULL,
  `dmag` double DEFAULT NULL,
  `limmag` double DEFAULT NULL,
  `filter` text,
  `issub` int(11) DEFAULT NULL COMMENT '0: Not template subtracted, 1: Template subtracted',
  `refsys` text,
  `observer` text,
  `reducer` text,
  `pipeline` text,
  `groupidcode` bigint(20) DEFAULT NULL,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  KEY `targetid` (`targetid`,`photlcoid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='All reduced photometry (including external sources)';



# Dump of table photlco
# ------------------------------------------------------------

DROP TABLE IF EXISTS `photlco`;

CREATE TABLE `photlco` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `objname` text COMMENT 'Object name used in the header',
  `dayobs` text COMMENT 'Date used to identify the night',
  `dateobs` date DEFAULT NULL COMMENT 'Date of observation',
  `ut` time DEFAULT NULL COMMENT 'UT time of observation',
  `mjd` double DEFAULT NULL COMMENT 'MJD of observation',
  `exptime` float DEFAULT NULL,
  `filter` varchar(10) DEFAULT NULL,
  `telescopeid` bigint(20) DEFAULT NULL,
  `instrumentid` bigint(20) DEFAULT NULL,
  `telescope` text,
  `instrument` text,
  `mag` double DEFAULT '9999',
  `dmag` double DEFAULT '9999',
  `airmass` float DEFAULT NULL,
  `wcs` float DEFAULT NULL COMMENT '0: Good solution, >=1: Error in solution',
  `psf` varchar(100) DEFAULT 'X' COMMENT 'Fits file containing the PSF',
  `apmag` double DEFAULT '9999' COMMENT 'Aperture photometry magnitude (instrumental)',
  `psfx` double DEFAULT '9999' COMMENT 'x-position of object center in image',
  `psfy` double DEFAULT '9999' COMMENT 'y-position of object center in image',
  `psfmag` double DEFAULT '9999' COMMENT 'PSF magnitude (instrumental)',
  `psfdmag` double DEFAULT '9999' COMMENT 'Error in psfmag',
  `z1` float DEFAULT '9999' COMMENT 'Zeropoint using zcol1',
  `z2` float DEFAULT '9999' COMMENT 'Zeropoint using zcol2',
  `zn` float DEFAULT '9999' COMMENT 'Zeropoint for the natural system',
  `c1` float DEFAULT '9999' COMMENT 'Colorterm for z1',
  `c2` float DEFAULT '9999' COMMENT 'Colorterm for z2',
  `znnum` float DEFAULT NULL COMMENT 'Number of stars used for zn',
  `dz1` float DEFAULT '9999' COMMENT 'Error for z1',
  `dz2` float DEFAULT '9999' COMMENT 'Error for z2',
  `dzn` float DEFAULT NULL COMMENT 'Error for zn',
  `dc1` float DEFAULT '9999' COMMENT 'Error for c1',
  `dc2` float DEFAULT '9999' COMMENT 'Error for c2',
  `zcol1` varchar(2) DEFAULT NULL COMMENT 'Color for z1',
  `zcol2` varchar(2) DEFAULT NULL COMMENT 'Color for z2',
  `quality` tinyint(1) DEFAULT '127' COMMENT '0: Bad, 1: Good',
  `zcat` varchar(100) DEFAULT 'X' COMMENT 'Catalog file used to determine the zeropoint',
  `abscat` varchar(100) DEFAULT 'X' COMMENT 'Catalog file used to calibrate mag',
  `fwhm` float DEFAULT '9999' COMMENT 'Units: pixels',
  `magtype` int(11) DEFAULT '1' COMMENT '-1: Non-detection limit, 1: None, 2: PSF fitting, 3: Aperture photometry',
  `ra0` double DEFAULT '9999' COMMENT 'RA in the image header',
  `dec0` double DEFAULT '9999' COMMENT 'Dec in the image header',
  `tracknumber` int(11) DEFAULT NULL,
  `filename` varchar(100) DEFAULT NULL COMMENT 'Image fits file',
  `difftype` tinyint(4) DEFAULT NULL COMMENT 'Type of difference algorithm used (0: Hotpants, 1: Optimal)',
  `filepath` text COMMENT 'Path to image file',
  `filetype` text,
  `groupidcode` bigint(20) DEFAULT NULL COMMENT 'Which groups can see this',
  `datecreated` timestamp NOT NULL DEFAULT '2000-01-01 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `apflux` double DEFAULT '9999' COMMENT 'Aperture flux',
  `dapflux` double DEFAULT '9999' COMMENT 'Aperture flux error',
  `dapmag` double DEFAULT '9999' COMMENT 'Aperture magnitude error',
  `limmag` double DEFAULT NULL COMMENT 'Three sigma limiting magnitude of the image as a whole',
  `lastunpacked` timestamp NULL DEFAULT '1970-01-01 00:00:01',
  `apercorr` float DEFAULT NULL COMMENT 'Aperture correction to correct psfmag to apmag',
  PRIMARY KEY (`id`),
  UNIQUE KEY `filename` (`filename`),
  KEY `targetid` (`targetid`),
  KEY `telescopeid` (`telescopeid`),
  KEY `instrumentid` (`instrumentid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Reduced LCOGT photometry';



# Dump of table photlcoraw
# ------------------------------------------------------------

DROP TABLE IF EXISTS `photlcoraw`;

CREATE TABLE `photlcoraw` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `objname` text COMMENT 'Object name used in the header',
  `dayobs` text COMMENT 'Date used to identify the night',
  `dateobs` date DEFAULT NULL COMMENT 'Date of observation',
  `ut` time DEFAULT NULL COMMENT 'UT time of observation',
  `mjd` double DEFAULT NULL COMMENT 'MJD of observation',
  `exptime` float DEFAULT NULL,
  `filter` varchar(10) DEFAULT NULL,
  `telescopeid` bigint(20) DEFAULT NULL,
  `instrumentid` bigint(20) DEFAULT NULL,
  `telescope` text,
  `instrument` text,
  `airmass` float DEFAULT NULL,
  `ra0` double DEFAULT '9999' COMMENT 'RA in the image header',
  `dec0` double DEFAULT '9999' COMMENT 'Dec in the image header',
  `cat_ra` double DEFAULT NULL,
  `cat_dec` double DEFAULT NULL,
  `temperature` float DEFAULT NULL,
  `propid` text,
  `obid` int(11) DEFAULT NULL,
  `rotskypa` float DEFAULT NULL,
  `userid` varchar(30) DEFAULT NULL,
  `fwhm` float DEFAULT '9999',
  `tracknumber` int(11) DEFAULT NULL,
  `filename` varchar(100) DEFAULT NULL COMMENT 'Image fits file',
  `filepath` text COMMENT 'Path to image file',
  `groupidcode` bigint(20) DEFAULT NULL COMMENT 'Which groups can see this',
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `moonfrac` float DEFAULT NULL,
  `moondist` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `filename` (`filename`),
  KEY `targetid` (`targetid`),
  KEY `tracknumber` (`tracknumber`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Raw LCOGT photometry';



# Dump of table photpairing
# ------------------------------------------------------------

DROP TABLE IF EXISTS `photpairing`;

CREATE TABLE `photpairing` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `namein` text COMMENT 'Filename of single image',
  `tablein` text COMMENT 'Tablename of single image',
  `nameout` text COMMENT 'Filename of combined/subtracted image',
  `tableout` text COMMENT 'Tablename of combined/subtracted image',
  `nametemplate` text COMMENT 'Filename of template image',
  `tabletemplate` text COMMENT 'Tablename of template image',
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Connecting combined / subtracted images to their original fi';



# Dump of table programs
# ------------------------------------------------------------

DROP TABLE IF EXISTS `programs`;

CREATE TABLE `programs` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idcode` varchar(20) NOT NULL DEFAULT '',
  `fullname` varchar(200) DEFAULT '',
  `shortname` varchar(50) DEFAULT NULL,
  `type` varchar(20) DEFAULT NULL,
  `username` varchar(50) DEFAULT NULL,
  `userpw` varchar(50) DEFAULT NULL,
  `groupidcode` bigint(20) DEFAULT NULL,
  `ingest` tinyint(1) DEFAULT NULL,
  `log_used_time` tinyint(1) DEFAULT '0',
  `programstart` timestamp NULL DEFAULT NULL,
  `programend` timestamp NULL DEFAULT NULL,
  `timecritical_idcode` varchar(20) DEFAULT NULL,
  `replacedwith` int(11) DEFAULT NULL,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `idcode` (`idcode`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table psns
# ------------------------------------------------------------

DROP TABLE IF EXISTS `psns`;

CREATE TABLE `psns` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(50) DEFAULT NULL,
  `ra` float DEFAULT NULL,
  `dec` float DEFAULT NULL,
  `discodate` datetime DEFAULT NULL,
  `discomag` float DEFAULT NULL,
  `discofilter` int(11) DEFAULT NULL,
  `offsete` float DEFAULT NULL,
  `offsetn` float DEFAULT NULL,
  `galaxy` varchar(50) DEFAULT NULL,
  `discoverer` char(11) DEFAULT NULL,
  `arc` char(11) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table reference_status
# ------------------------------------------------------------

DROP TABLE IF EXISTS `reference_status`;

CREATE TABLE `reference_status` (
  `id` bigint(11) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(11) DEFAULT NULL,
  `status` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table scheduling_run
# ------------------------------------------------------------

DROP TABLE IF EXISTS `scheduling_run`;

CREATE TABLE `scheduling_run` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `last_successful_run` datetime DEFAULT NULL,
  `failure_alert_sent` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table schedulinglog
# ------------------------------------------------------------

DROP TABLE IF EXISTS `schedulinglog`;

CREATE TABLE `schedulinglog` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `runstart` timestamp NULL DEFAULT NULL,
  `sequenceid` bigint(20) NOT NULL,
  `action` text NOT NULL,
  `datecreated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `notvisible` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `datecreated_schedlog_index` (`datecreated`),
  KEY `notvisible` (`notvisible`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Observation Scheduling Log';



# Dump of table spec
# ------------------------------------------------------------

DROP TABLE IF EXISTS `spec`;

CREATE TABLE `spec` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `objname` varchar(50) NOT NULL,
  `dayobs` text COMMENT 'Date used to identify the night',
  `dateobs` date DEFAULT NULL COMMENT 'Date of observation',
  `ut` time DEFAULT NULL COMMENT 'UT time of observation',
  `mjd` double DEFAULT NULL COMMENT 'MJD of observation',
  `exptime` float DEFAULT NULL,
  `filter` varchar(10) DEFAULT NULL,
  `grism` varchar(20) DEFAULT NULL,
  `telescopeid` bigint(20) DEFAULT NULL,
  `instrumentid` bigint(20) DEFAULT NULL,
  `telescope` text,
  `instrument` text,
  `airmass` float DEFAULT NULL,
  `slit` varchar(20) DEFAULT NULL,
  `original` varchar(100) DEFAULT NULL COMMENT 'Raw data filename',
  `ra0` double DEFAULT '9999' COMMENT 'RA in the image header',
  `dec0` double DEFAULT '9999' COMMENT 'Dec in the image header',
  `observer` text,
  `reducer` text,
  `filename` varchar(100) DEFAULT NULL COMMENT 'Spectrum fits file',
  `filepath` text COMMENT 'Path to image file',
  `groupidcode` bigint(20) DEFAULT NULL COMMENT 'Which groups can see this',
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `uploadedby` int(11) DEFAULT '0',
  `deleted` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  UNIQUE KEY `filename` (`filename`),
  KEY `objname` (`objname`),
  KEY `original` (`original`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='All reduced spectroscopy (including external sources)';



# Dump of table speclcoguider
# ------------------------------------------------------------

DROP TABLE IF EXISTS `speclcoguider`;

CREATE TABLE `speclcoguider` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `tracknumber` int(11) DEFAULT NULL,
  `blockid` int(11) DEFAULT NULL,
  `link` text,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table speclcoraw
# ------------------------------------------------------------

DROP TABLE IF EXISTS `speclcoraw`;

CREATE TABLE `speclcoraw` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `objname` varchar(50) DEFAULT NULL,
  `dayobs` text COMMENT 'Date used to identify the night',
  `dateobs` date DEFAULT NULL COMMENT 'Date of observation',
  `ut` time DEFAULT NULL COMMENT 'UT time of observation',
  `mjd` double DEFAULT NULL COMMENT 'MJD of observation',
  `exptime` float DEFAULT NULL,
  `filter` varchar(10) DEFAULT NULL,
  `grism` varchar(20) DEFAULT NULL,
  `telescopeid` bigint(20) DEFAULT NULL,
  `instrumentid` bigint(20) DEFAULT NULL,
  `telescope` text,
  `instrument` text,
  `type` varchar(20) DEFAULT NULL,
  `airmass` float DEFAULT NULL,
  `slit` varchar(20) DEFAULT NULL,
  `lamp` varchar(20) DEFAULT NULL,
  `ra0` double DEFAULT '9999' COMMENT 'RA in the image header',
  `dec0` double DEFAULT '9999' COMMENT 'Dec in the image header',
  `cat_ra` double DEFAULT NULL,
  `cat_dec` double DEFAULT NULL,
  `temperature` float DEFAULT NULL,
  `propid` varchar(30) DEFAULT NULL,
  `obid` int(11) DEFAULT NULL,
  `rotskypa` float DEFAULT NULL,
  `observer` varchar(30) DEFAULT NULL,
  `userid` varchar(30) DEFAULT NULL,
  `fwhm` float DEFAULT '9999',
  `tracknumber` int(11) DEFAULT NULL,
  `filename` varchar(100) DEFAULT NULL COMMENT 'Spectrum fits file',
  `filepath` text COMMENT 'Path to image file',
  `groupidcode` bigint(20) DEFAULT NULL COMMENT 'Which groups can see this',
  `markedasbad` tinyint(1) DEFAULT '0',
  `markedasbad_by` varchar(50) DEFAULT NULL,
  `markedasbad_because` text,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `filename` (`filename`),
  KEY `objname` (`objname`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Raw LCOGT spectroscopy';



# Dump of table tags
# ------------------------------------------------------------

DROP TABLE IF EXISTS `tags`;

CREATE TABLE `tags` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `tag` varchar(50) NOT NULL,
  `userid` bigint(20) NOT NULL COMMENT 'User who created this tag',
  `datecreated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Tags of science projects to be used when scheduling observat';



# Dump of table targetnames
# ------------------------------------------------------------

DROP TABLE IF EXISTS `targetnames`;

CREATE TABLE `targetnames` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `targetid` bigint(20) NOT NULL,
  `name` varchar(50) NOT NULL,
  `groupidcode` bigint(20) DEFAULT NULL COMMENT 'Which groups can see this',
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Target names';



# Dump of table targets
# ------------------------------------------------------------

DROP TABLE IF EXISTS `targets`;

CREATE TABLE `targets` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `ra0` double DEFAULT NULL,
  `dec0` double DEFAULT NULL,
  `redshift` double DEFAULT NULL,
  `classification` text,
  `classificationid` int(11) DEFAULT NULL,
  `sloan_cat` text,
  `landolt_cat` text,
  `apass_cat` text,
  `pm_ra` double NOT NULL DEFAULT '0' COMMENT 'proper motion RA times cos(dec)',
  `pm_dec` double NOT NULL DEFAULT '0' COMMENT 'proper motion Dec',
  `groupidcode` bigint(20) NOT NULL DEFAULT '32769' COMMENT 'Which groups can see this',
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Info about targets';



# Dump of table telescopes
# ------------------------------------------------------------

DROP TABLE IF EXISTS `telescopes`;

CREATE TABLE `telescopes` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `name` text NOT NULL,
  `shortname` text,
  `lat` double DEFAULT NULL,
  `lon` double DEFAULT NULL,
  `elevation` double DEFAULT NULL,
  `diameter` double DEFAULT NULL,
  `webpage` text,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table timecharged
# ------------------------------------------------------------

DROP TABLE IF EXISTS `timecharged`;

CREATE TABLE `timecharged` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `datetime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `programidcode` varchar(20) DEFAULT NULL,
  `1m_std_used` float DEFAULT NULL,
  `1m_std_allocated` float DEFAULT NULL,
  `1m_too_used` float DEFAULT NULL,
  `1m_too_allocated` float DEFAULT NULL,
  `2m_std_used` float DEFAULT NULL,
  `2m_std_allocated` float DEFAULT NULL,
  `2m_too_used` float DEFAULT NULL,
  `2m_too_allocated` float DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table useractionlog
# ------------------------------------------------------------

DROP TABLE IF EXISTS `useractionlog`;

CREATE TABLE `useractionlog` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `userid` bigint(20) NOT NULL,
  `tablemodified` varchar(50) DEFAULT NULL,
  `idmodified` bigint(20) DEFAULT NULL,
  `columnmodified` varchar(50) DEFAULT NULL,
  `prevvalue` text,
  `newvalue` text,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Log of user actions';



# Dump of table userrequests
# ------------------------------------------------------------

DROP TABLE IF EXISTS `userrequests`;

CREATE TABLE `userrequests` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `firstname` text NOT NULL,
  `lastname` text NOT NULL,
  `email` text NOT NULL,
  `groupidcode` bigint(20) DEFAULT NULL,
  `othergroups` text,
  `reason` text,
  `datecreated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `granteduserid` bigint(20) DEFAULT NULL COMMENT '0 = rejected',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Pending Users';



# Dump of table users
# ------------------------------------------------------------

DROP TABLE IF EXISTS `users`;

CREATE TABLE `users` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(20) NOT NULL,
  `pw` text,
  `groupidcode` bigint(20) DEFAULT NULL,
  `firstname` text,
  `lastname` text,
  `email` text,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='User info';



# Dump of table voevent_amon
# ------------------------------------------------------------

DROP TABLE IF EXISTS `voevent_amon`;

CREATE TABLE `voevent_amon` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `ivorn` text NOT NULL,
  `role` varchar(45) NOT NULL,
  `version` text NOT NULL,
  `xmlns_voe` text NOT NULL,
  `xmlns_xsi` text NOT NULL,
  `xsi_schemalocation` text NOT NULL,
  `author_ivorn` text NOT NULL,
  `shortname` text,
  `contactname` text,
  `contactemail` text,
  `who_description` text,
  `date` text NOT NULL,
  `packet_type` text NOT NULL,
  `pkt_ser_num` text NOT NULL,
  `trig_id` text NOT NULL,
  `event_tjd` text NOT NULL,
  `event_sod` text NOT NULL,
  `nevents` text NOT NULL,
  `stream` text NOT NULL,
  `rev` varchar(45) NOT NULL,
  `false_pos` varchar(45) NOT NULL,
  `pvalue` varchar(45) NOT NULL,
  `deltat` text NOT NULL,
  `sigmat` text NOT NULL,
  `charge` text NOT NULL,
  `signalness` text NOT NULL,
  `hesetypeindex` text NOT NULL,
  `trigger_id` varchar(45) NOT NULL,
  `misc_flags` varchar(45) DEFAULT NULL,
  `subtype` text,
  `test` varchar(45) DEFAULT NULL,
  `radec_valid` varchar(45) DEFAULT NULL,
  `retraction` varchar(45) DEFAULT NULL,
  `internal_test` text NOT NULL,
  `observatorylocation_id` text NOT NULL,
  `astrocoordsystem_id` text NOT NULL,
  `timeunit` varchar(45) NOT NULL,
  `isotime` text NOT NULL,
  `ra0` text,
  `dec0` text,
  `error2radius` text,
  `how_description` text,
  `reference_uri` text,
  `importance` text,
  `inference_probability` text,
  `concept` text,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table voevent_lvc
# ------------------------------------------------------------

DROP TABLE IF EXISTS `voevent_lvc`;

CREATE TABLE `voevent_lvc` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `ivorn` text,
  `role` varchar(45) DEFAULT NULL,
  `version` text,
  `xmlns_voe` text,
  `xmlns_xsi` text,
  `xsi_schemalocation` text,
  `author_ivorn` text,
  `shortname` text,
  `contactname` text,
  `contactemail` text,
  `who_description` text,
  `date` text,
  `packet_type` text,
  `pkt_ser_num` text,
  `alert_type` text,
  `graceid` text,
  `id_letter` text,
  `trig_id` text,
  `trigger_tjd` text,
  `trigger_sod` text,
  `eventpage` text,
  `_group` text,
  `search` varchar(45) DEFAULT NULL,
  `pipeline` varchar(45) DEFAULT NULL,
  `internal` varchar(45) DEFAULT NULL,
  `far` text,
  `chirpmass` text,
  `eta` text,
  `maxdistance` text,
  `probhasns` text,
  `probhasremnant` text,
  `hardwareinj` text,
  `vetted` text,
  `openalert` text,
  `temporalcoinc` text,
  `trigger_id` varchar(45) DEFAULT NULL,
  `misc_flags` varchar(45) DEFAULT NULL,
  `lvc_internal` text,
  `test` varchar(45) DEFAULT NULL,
  `retraction` varchar(45) DEFAULT NULL,
  `internal_test` varchar(45) DEFAULT NULL,
  `num_det_participated` text,
  `lho_participated` varchar(45) DEFAULT NULL,
  `llo_participated` varchar(45) DEFAULT NULL,
  `virgo_participated` varchar(45) DEFAULT NULL,
  `geo600_participated` varchar(45) DEFAULT NULL,
  `kagra_participated` varchar(45) DEFAULT NULL,
  `lio_participated` varchar(45) DEFAULT NULL,
  `sequence_number` text,
  `skymap_url_fits_basic` text,
  `observatorylocation_id` text,
  `astrocoordsystem_id` text,
  `timeunit` varchar(45) DEFAULT NULL,
  `isotime` text,
  `how_description` text,
  `reference_uri` text,
  `importance` text,
  `inference_probability` text,
  `concept` text,
  `datecreated` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;



# Dump of table ztf_alerts
# ------------------------------------------------------------

DROP TABLE IF EXISTS `ztf_alerts`;

CREATE TABLE `ztf_alerts` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `ztf_targetid` int(11) DEFAULT NULL,
  `jd` double DEFAULT NULL,
  `avro` text,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



# Dump of table ztf_targets
# ------------------------------------------------------------

DROP TABLE IF EXISTS `ztf_targets`;

CREATE TABLE `ztf_targets` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(50) DEFAULT NULL,
  `ra0` double DEFAULT NULL,
  `dec0` double DEFAULT NULL,
  `interest` varchar(50) DEFAULT NULL,
  `snoozed` int(11) NOT NULL DEFAULT '0',
  `targetid` int(11) DEFAULT NULL,
  `lastmodified` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;




/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;
/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
