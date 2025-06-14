#!/bin/bash

mkdir -p "$NXF_HOME/plugins"

cd "$NXF_HOME/plugins"

for plugin in jfouret/nf-schema@2.4.2-dev nexomis/nf-nexomis@0.1.0; do
    repo=$(echo $plugin | cut -d'@' -f1)
    org=$(echo $repo | cut -d'/' -f1)
    name=$(echo $repo | cut -d'/' -f2)
    ver=$(echo $plugin | cut -d'@' -f2)
    rm -rf $plugin
    wget https://github.com/${org}/${name}/releases/download/${ver}/${name}-${ver}.zip
    unzip ${name}-${ver}.zip -d ${name}-${ver}
    rm ${name}-${ver}.zip
    rm -rf $plugin
done
