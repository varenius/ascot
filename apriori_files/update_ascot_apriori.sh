EMAIL_ADDR=eskil.varenius@chalmers.se
ASCOT_APRIORI=/opt/ascot/apriori_files/
#ECCDAT.ecc not updated automatically due to manual modifications (added dome numbers for OTT)
curl https://raw.githubusercontent.com/anothnagel/antenna-info/master/antenna-info.txt > $ASCOT_APRIORI/antenna-info.txt
curl ftp://hpiers.obspm.fr/iers/series/opa/eopc04_IAU2000 > $ASCOT_APRIORI/eopc04_IAU2000
curl ftp://hpiers.obspm.fr/iers/series/opa/eopc04 > $ASCOT_APRIORI/eopc04
curl https://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now > $ASCOT_APRIORI/eopc04_IAU2000.62-now
curl ftp://ivs.bkg.bund.de/pub/vlbi/ivscontrol/ns-codes.txt > $ASCOT_APRIORI/ns-codes.txt
# ns-codes.txt need to be in two places, should be fixed in ASCOT...
curl ftp://ivs.bkg.bund.de/pub/vlbi/ivscontrol/ns-codes.txt > $ASCOT_APRIORI/masterfiles/ns-codes.txt
curl ftp://ftp.iers.org/products/eop/rapid/standard/finals2000A.all > $ASCOT_APRIORI/finals2000A.all
curl -u anonymous:$EMAIL_ADDR --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/IVS_SrcNamesTable.txt > $ASCOT_APRIORI/IVS_SrcNamesTable.txt
curl -u anonymous:$EMAIL_ADDR --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp > $ASCOT_APRIORI/usno_finals.erp
curl -u anonymous:$EMAIL_ADDR --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/blokq.c11.dat > $ASCOT_APRIORI/blokq.c11.dat
# Extract ocean loading data from blokq file
/usr/bin/python3 blokq.c11_to_ol.py $ASCOT_APRIORI/blokq.c11.dat  $ASCOT_APRIORI/blokq.c11.ol
# Update VMF1 and VMF3 files
/opt/ascot/bin/get_ext_data -d $ASCOT_APRIORI

# Get all masterfiles
startyear=1979
endyear=$(date +%Y)
for i in $(seq $startyear $endyear)
  do
    yy=${i:2:2}
    echo $yy
    # Get masterfile for this year
    curl ftp://ivs.bkg.bund.de/pub/vlbi/ivscontrol/master$yy.txt > $ASCOT_APRIORI/masterfiles/master$yy.txt
    # Get INT-masterfile for this year, if it exists
    inturl="ftp://ivs.bkg.bund.de/pub/vlbi/ivscontrol/master$yy-int.txt"
    if curl --output /dev/null --silent --head --fail "$inturl"; then
      curl $inturl > $ASCOT_APRIORI/masterfiles/master$yy-int.txt
    else
      echo "URL does not exist: $inturl"
    fi
done
