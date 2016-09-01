#!/bin/bash


if [ -z $NUSKYBGD ]
then
	echo "Initializing nuskybgd environment."

	if [[ $SHELL == *"bash"* ]]
	then	
		echo
		echo "Adding setup lines to .bash_profile"
		echo "" >> ~/.bash_profile
		export NUSKYBGD=`pwd`
		echo 
		echo

		echo "## Lines added by nuskbygd init script: ##" >> ~/.bash_profile
		echo "export NUSKYBGD=${NUSKYBGD}" >> ~/.bash_profile
		echo "export NUSKYBGD_AUXIL=${NUSKYBGD}/auxil" >> ~/.bash_profile
		echo "## ##" >> ~/.bash_profile


		echo "Setting environment variaibles: "
		echo "export NUSKYBGD=${NUSKYBGD}"
		echo "export NUSKYBGD_AUXIL=${NUSKYBGD}/auxil"
		echo
	
		export NUSKYBGD_AUXIL=${NUSKYBGD}/auxil		
	else
		echo "Please use bash for now or set your own environment as in this script."
	fi
else 
	echo "\$NUSKYBGD is already set to: $NUSKYBGD"
fi





