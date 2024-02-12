#!/bin/bash

setup(){

    export MODE="NUMI"
    export BOOSTROOT=${BOOST_DIR}/source/boost_1_66_0
    #DK2NU:
    export DK2NU_INC=${DK2NU}/include/dk2nu/tree
    export DK2NU_LIB=${DK2NU}/lib

    # setup for jobsub client
    # according to the prescription in Mike Kirby's talk
    # minerva doc-10551, Dec 2014 (same doc can be found for other experiments)
    export NUMIANA_DIR=${PWD}
    export NUMIANA_INC=${NUMIANA_DIR}/include
    export LD_LIBRARY_PATH=${NUMIANA_DIR}/lib:${NUMIANA_DIR}/dict:${PPFX_DIR}/lib:${PPFX_DIR}/slf7.x86_64.e17.prof/lib:$LD_LIBRARY_PATH
    export LIBRARY_PATH=$LIBRARY_PATH:$LD_LIBRARY_PATH
}
HOST=$(hostname -f)
echo $HOST
if echo "$HOST" | grep 'uboone';then
    echo "executing for the $HOST"
    setup
else
    echo "This is not a uboone machine. Sorry!"
fi

