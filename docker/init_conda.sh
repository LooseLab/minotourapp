conda activate artic
python setup.py install
artic -v
conda activate pangolin
pip install .
pangolin --version
conda deactivate