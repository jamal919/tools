#! /bin/bash

#
# tg_update_gloss_db.sh - get a brand new copy of the GLOSS DB
#
# Author: Sara Fleury  CTOH/LEGOS
# Date: Thu Apr 18 22:37:12 CEST 2013
#

SYNTAX="
SYNTAX: tg_update_gloss_db.sh <out_path>

        To be called from a LEGOS machin

EXEMPLE: 
        tg_update_gloss_db.sh  fles@nuwa.aero.obs-mip.fr:/gsa/fles/gloss/

"

# GLOSS DB at LEGOS (eg, from 'malo')
GLOSS_PATH="/data/soa/niveau_mer/gloss/data/hourly/"

#------------------------------------------------------------
# get params

if [ $# -ne 1 ] ; then
    printf "$SYNTAX"
    exit
fi

OUT_PATH=$1
DATE=`date +%Y%m%d`

cmd="rsync -pr ${GLOSS_PATH} ${OUT_PATH}/hourly-${DATE}/origin/"

echo
echo "CMD TO EXECUTE FROM LEGOS:" $cmd
echo "let's try ..."
echo
$cmd
