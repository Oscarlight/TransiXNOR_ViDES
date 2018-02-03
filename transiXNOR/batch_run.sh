version=D9
for vds in 0.4 
	do
		# python device.py $version --vdsmin $vds --vdsmax $vds --vdsN 1 &
		# let vds-=0.01
		# python device.py $version --vdsmin $vds --vdsmax $vds --vdsN 1 &
		# let vds-=0.01
		# python device.py $version --vdsmin $vds --vdsmax $vds --vdsN 1 &
		let vds-=0.01
		python device.py $version --vdsmin $vds --vdsmax $vds --vdsN 1
		let vds-=0.01
		wait
	done
