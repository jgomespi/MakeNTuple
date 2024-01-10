#!/bin/bash

# Call CMSDriver in an even more intelligent way

if [[ $1 == *"chain"* ]]; 
then
  echo "This is a chained request"
  wget --no-check-certificate https://cms-pdmv.cern.ch/mcm/public/restapi/chained_requests/get_test/$1 
else
  wget --no-check-certificate  https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_test/$1
fi
mv $1 $1.sh

### change number of events for testing
sed /^cmsDriv/s/"||"/"-n $2 ||"/g $1.sh > $1_new.sh
mv $1_new.sh $1.sh
### redirect all outputs not to clog AFS area
sed s/"file:"/"file:\/tmp\/$USER\/"/g $1.sh > $1_new.sh
mv $1_new.sh $1.sh
### remove all DIGI-RECO, DQM crap etc. if it is a chain
if [[ `more $1.sh | grep nothing` ]];
then
  head -n `awk '/nothing/{ print NR-1; exit }' $1.sh` $1.sh > $1_new.sh
  mv $1_new.sh $1.sh 
fi
 
if [[ $3 == 1 ]]; 
then
   sed s/GEN,SIM/GEN/g $1.sh > $1_new.sh
   mv $1_new.sh $1.sh
fi

source $1.sh


