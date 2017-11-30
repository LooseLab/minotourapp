python ../manage.py dumpdata --indent 2 \
    reads.fastqreadtype \
    reads.groupruntype \
    reads.jobtype \
    reads.minioneventtype \
    > ../fixtures/auxiliary_data.json
