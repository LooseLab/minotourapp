cd /var/lib/minotour/apps/fieldbioinformatics/
conda activate artic
python setup.py install
artic -v;
cd /var/lib/minotour/apps/pangolin/
conda activate pangolin
pip install .
pangolin --version
conda deactivate;
