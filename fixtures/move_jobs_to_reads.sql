CREATE TABLE `reads_jobtype` (`id` integer AUTO_INCREMENT NOT NULL PRIMARY KEY, `name` varchar(256) NOT NULL, `description` longtext NULL, `long_description` longtext NULL, `reference` bool NOT NULL, `transcriptome` bool NOT NULL, `readcount` bool NOT NULL, `private` bool NOT NULL);

CREATE TABLE `reads_jobmaster` (`id` integer AUTO_INCREMENT NOT NULL PRIMARY KEY, `last_read` bigint NOT NULL, `tempfile_name` varchar(256) NULL, `read_count` bigint NOT NULL, `complete` bool NOT NULL, `running` bool NOT NULL, `paused` bool NOT NULL, `iteration_count` integer NULL, `target_set` varchar(100) NULL, `flowcell_id` integer NULL, `job_type_id` integer NOT NULL, `reference_id` integer NULL, `run_id` integer NULL);

ALTER TABLE `reads_jobmaster` ADD CONSTRAINT `reads_jobmaster_flowcell_id_49c1aab0_fk_reads_flowcell_id` FOREIGN KEY (`flowcell_id`) REFERENCES `reads_flowcell` (`id`);

ALTER TABLE `reads_jobmaster` ADD CONSTRAINT `reads_jobmaster_job_type_id_8b9e2415_fk_reads_jobtype_id` FOREIGN KEY (`job_type_id`) REFERENCES `reads_jobtype` (`id`);

ALTER TABLE `reads_jobmaster` ADD CONSTRAINT `reads_jobmaster_reference_id_42760c83_fk_reference` FOREIGN KEY (`reference_id`) REFERENCES `reference_referenceinfo` (`id`);

ALTER TABLE `reads_jobmaster` ADD CONSTRAINT `reads_jobmaster_run_id_d84205cf_fk_reads_run_id` FOREIGN KEY (`run_id`) REFERENCES `reads_run` (`id`);

insert into reads_jobtype select * from jobs_jobtype;

insert into reads_jobmaster select * from jobs_jobmaster;

