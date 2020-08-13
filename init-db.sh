#!/bin/sh
LCOSNDBNAME=${LCOSNDBNAME:-supernova};
LCODBSNUSER=${LCOSNDBUSER:-supernova}
LCOSNDBPASS=${LCOSNDBPASS:-supernova}

envsubst < ${LCOSNPIPE:-/lcogtsnpipe}/supernova.sql > tmp.sql
mysql -h ${LCOSNDBHOST:-supernovadb} -u root -p${LCOSNDBROOTPASS-password} < tmp.sql
rm -rf tmp.sql
