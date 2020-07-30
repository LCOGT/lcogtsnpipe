#!/bin/sh
LCOSNDBNAME=${LCOSNDBNAME:-supernova};
LCODBSNUSER=${LCOSNDBUSER:-supernova}
LCOSNDBPASS=${LCOSNDBPASS:-supernova}

envsubst < /supernova/github/lcogtsnpipe/supernova.sql > test.sql
mysql -h ${LCOSNDBHOST:-supernovadb} -u root -p${LCOSNDBROOTPASS-password} < test.sql
