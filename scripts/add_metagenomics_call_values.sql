--
-- Create model CAllValues
--
CREATE TABLE `metagenomics_callvalues` (`id` integer AUTO_INCREMENT NOT NULL PRIMARY KEY, `classification_count` integer NOT NULL, `classification_number` integer NOT NULL, `task_id` integer NOT NULL);
ALTER TABLE `metagenomics_callvalues` ADD CONSTRAINT `metagenomics_callvalues_task_id_bb805c9e_fk_reads_jobmaster_id` FOREIGN KEY (`task_id`) REFERENCES `reads_jobmaster` (`id`);
