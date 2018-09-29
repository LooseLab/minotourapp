python ../manage.py dumpdata --indent 2 \
    reads.fastqreadtype \
    jobs.jobtype \
    reads.minioneventtype \
    reference.referenceinfo \
    reference.referenceline \
    > ../fixtures/auxiliary_data.json
