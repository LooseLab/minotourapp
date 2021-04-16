--
-- Create model UncertaintyProbability
--
CREATE TABLE `metagenomics_uncertaintyprobability` (`id` integer AUTO_INCREMENT NOT NULL PRIMARY KEY, `upper_ci_value` double precision NOT NULL, `lower_ci_value` double precision NOT NULL, `classification_number` integer NOT NULL, `tax_ids` varchar(256) NULL, `task_id` integer NOT NULL);
ALTER TABLE `metagenomics_uncertaintyprobability` ADD CONSTRAINT `metagenomics_uncerta_task_id_545a3579_fk_reads_job` FOREIGN KEY (`task_id`) REFERENCES `reads_jobmaster` (`id`);
