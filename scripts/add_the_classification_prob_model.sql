--
-- Create model ClassificationProbabilites
--
CREATE TABLE `metagenomics_classificationprobabilites` (`id` integer AUTO_INCREMENT NOT NULL PRIMARY KEY, `tax_ids` integer NOT NULL, `probability` double precision NOT NULL, `classification_number` integer NOT NULL, `target_set_name` varchar(256) NOT NULL, `task_id` integer NOT NULL);
ALTER TABLE `metagenomics_classificationprobabilites` ADD CONSTRAINT `metagenomics_classif_task_id_7dd21256_fk_reads_job` FOREIGN KEY (`task_id`) REFERENCES `reads_jobmaster` (`id`);

--
-- Alter field tax_ids on classificationprobabilites
--
ALTER TABLE `metagenomics_classificationprobabilites` MODIFY `tax_ids` varchar(256) NOT NULL;

--
-- Remove field task from classificationprobabilites
--
ALTER TABLE `metagenomics_classificationprobabilites` DROP FOREIGN KEY `metagenomics_classif_task_id_7dd21256_fk_reads_job`;
ALTER TABLE `metagenomics_classificationprobabilites` DROP COLUMN `task_id`;
