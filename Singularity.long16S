Bootstrap: docker
From: ubuntu:16.04

%post
	apt update -y
	apt upgrade -y
	apt install -y wget bzip2 zip git python3 python3-pip r-base libz-dev bc xvfb software-properties-common
	apt clean -y
	pip3 install --upgrade pip

	# Install OpenJDK 12
	wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public |  apt-key add -
	add-apt-repository --yes https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/
	apt install -y apt-transport-https
	apt update -y
	apt install -y adoptopenjdk-12-hotspot openjfx
	apt autoremove -y

	## Install NanoPlot nanofilt and (deprecated -> qcat)
	pip install NanoPlot==1.28.4 nanofilt # qcat

	## Install porechop
	pip install git+https://github.com/rrwick/Porechop.git

	## Install LAST
	cd /opt
	wget http://last.cbrc.jp/last-1047.zip
	unzip last-1047.zip
	cd last-1047
	make
	make install
	cd ..

	## Install megan6
	wget https://software-ab.informatik.uni-tuebingen.de/download/megan6/MEGAN_Community_unix_6_18_0.sh
	bash MEGAN_Community_unix_6_18_0.sh -q 
	MEGAN="/opt/megan"
	find "$MEGAN"/tools -type f | while read -r file; do b=$(basename "$file"); ln -fs "$file" /usr/local/bin/"$b"; sed -i -e "s@\"\$0\"@\"$file\"@" "$file"; done
	ln -fs "$MEGAN"/MEGAN /usr/local/bin/MEGAN
	sed -i -e "s@\"\$0\"@\"$MEGAN/MEGAN\"@" "$MEGAN"/MEGAN

	## Download DAA_Converter (jar) and make a executable script to call it
	wget https://github.com/BenjaminAlbrecht84/DAA_Converter/releases/download/v0.9.0/DAA_Converter_v0.9.0.jar

	## Clean
	rm -rf MEGAN_Community_unix_6_15_0.sh last-1047 last-1047.zip

%labels
	Authors: Cecilia Salazar (csalazar@pasteur.edu.uy) & Ignacio Ferres (iferres@pasteur.edu.uy)
	Maintainer: Ignacio Ferres (iferres@pasteur.edu.uy)
	
