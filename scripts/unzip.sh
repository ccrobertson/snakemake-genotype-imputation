password=$1
for i in {1..22}; do
	echo chr$i
	unzip -P ${password} chr_$i.zip
done
