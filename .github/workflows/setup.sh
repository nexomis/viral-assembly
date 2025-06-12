#!/bin/bash

NXF_PLUGINS_TEST_REPOSITORY=https://github.com/jfouret/nf-schema/releases/download/2.4.2-dev/nf-schema-2.4.2-dev-meta.json nextflow plugin install nf-schema@2.4.2-dev
NXF_PLUGINS_TEST_REPOSITORY=https://github.com/nexomis/nf-nexomis/releases/download/0.1.0/nf-schema-0.1.0-meta.json nextflow plugin install nf-nexomis@0.1.0

