--
-- Create model EstimatedAbundance
--
CREATE TABLE `metagenomics_estimatedabundance` (`id` integer AUTO_INCREMENT NOT NULL PRIMARY KEY, `tax_id` integer NOT NULL, `abundance` double precision NOT NULL, `name` varchar(256) NULL, `task_id` integer NOT NULL);
ALTER TABLE `metagenomics_estimatedabundance` ADD CONSTRAINT `metagenomics_estimat_task_id_7655d94a_fk_reads_job` FOREIGN KEY (`task_id`) REFERENCES `reads_jobmaster` (`id`);
