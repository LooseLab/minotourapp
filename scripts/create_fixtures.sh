python ../manage.py dumpdata --indent 2 \
    reads.fastqreadtype \
    jobs.jobtype \
    reads.minioneventtype \
    > ../fixtures/auxiliary_data.json
