python ../manage.py dumpdata --indent 2 \
    reads.fastqreadtype \
    reads.jobtype \
    reads.minioneventtype \
    auth.user \
    authtoken \
    > ../fixtures/auxiliary_data.json
