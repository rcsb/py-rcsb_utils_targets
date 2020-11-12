#!/bin/bash
# File: Pharos_dump_compound_tables.sh
#
# Dump tables from the Pharos schema (TCRDv6.7) containing ChEMBL references...
#
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "drop database if exists tcrd6; create database tcrd6;"
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW tcrd6 < latest.sql.gz
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "use tcrd6; select * from drug_activity;" >  drug_activity.tdd
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "use tcrd6; select * from cmpd_activity;" >  cmpd_activity.tdd
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "use tcrd6; select * from target;" >  target.tdd
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "use tcrd6; select * from protein;" >  protein.tdd
