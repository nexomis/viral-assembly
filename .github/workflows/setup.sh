#!/bin/bash

mkdir -p "$NXF_HOME/plugins"

for plugin in jfouret/nf-schema@2.4.2-dev nexomis/nf-nexomis@0.1.0; do
    repo=$(echo $plugin | cut -d'@' -f1)
    org=$(echo $repo | cut -d'/' -f1)
    name=$(echo $repo | cut -d'/' -f2)
    ver=$(echo $plugin | cut -d'@' -f2)
    wget https://github.com/${org}/${name}/releases/download/${ver}/${name}-${ver}.zip -O $NXF_HOME/plugins/${name}-${ver}.zip
    unzip $NXF_HOME/plugins/${name}-${ver}.zip -d $NXF_HOME/plugins/${name}-${ver}
    rm $NXF_HOME/plugins/${name}-${ver}.zip
done

bash data/test/inputs/get_data_test.sh
