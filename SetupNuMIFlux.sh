source /nusoft/app/externals/setup
if(env | grep -q ^ROOT)


then 
  echo "ROOT has already been setup."
  echo "If the version is not correct, this may cause errors."
  echo "If this occurs, try starting a fresh terminal session."
else
  setup root   v5_34_18a   -q e5:debug:nu
fi


if(env | grep -q ^GENIEXSECPATH)
then
  echo "GENIEXSECPATH has already been setup."
  echo "If the version is not correct, this may cause errors."
  echo "If this occurs, try starting a fresh terminal session."
else
  setup genie_xsec R-2_8_0   -q default
fi 

