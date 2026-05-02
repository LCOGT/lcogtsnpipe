#!/bin/sh
LCOSNDBNAME=${LCOSNDBNAME:-supernova};
LCODBSNUSER=${LCOSNDBUSER:-supernova}
LCOSNDBPASS=${LCOSNDBPASS:-supernova}

envsubst < ${LCOSNPIPE:-/lcogtsnpipe}/supernova.sql > tmp.sql
# MySQL client behavior differs across distros/architectures; force non-SSL for local Docker DB.
if mysql --help 2>/dev/null | grep -q -- '--ssl-mode'; then
	MYSQL_SSL_OPTS='--ssl-mode=DISABLED'
else
	MYSQL_SSL_OPTS='--ssl=0'
fi

mysql ${MYSQL_SSL_OPTS} -h ${LCOSNDBHOST:-supernovadb} -u root -p${LCOSNDBROOTPASS:-password} < tmp.sql
rm -rf tmp.sql
