#!/bin/bash
# File: ChEMBL_dump_compound_tables.sh
#
# Dump tables from the ChEMBL schema (chembl_27) ...
#
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "drop database if exists chembl_27;"
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "create database chembl_27 DEFAULT CHARACTER SET utf8 DEFAULT COLLATE utf8_general_ci;"
mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW chembl_27 < chembl_27_mysql.dmp
#mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "use chembl_27; select * from drug_activity;" >  drug_activity.tdd
#mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "use chembl_27; select * from cmpd_activity;" >  cmpd_activity.tdd
#mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "use chembl_27; select * from target;" >  target.tdd
#mysql --user=$MYSQL_DB_USER --password=$MYSQL_DB_PW  -e "use chembl_27; select * from protein;" >  protein.tdd
