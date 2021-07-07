#!/bin/bash
conda activate artic
if ! command -v artic &> /dev/null
then
  cd /var/lib/minotour/apps/fieldbioinformatics/
  echo "artic could not be found, installing" >&2
  python setup.py install
fi
artic -v;
conda activate pangolin
if ! command -v pangolin &> /dev/null
then
  echo "pangolin command not found, installing" >&2
  cd /var/lib/minotour/apps/pangolin/
  pip install .
fi
pangolin --version
conda deactivate;
cd /var/lib/minotour/apps/minotourapp