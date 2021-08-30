--
-- Add field max_read_len to primerscheme
--
ALTER TABLE `reads_primerscheme` ADD COLUMN `max_read_len` integer DEFAULT 800 NOT NULL;
ALTER TABLE `reads_primerscheme` ALTER COLUMN `max_read_len` DROP DEFAULT;
--
-- Add field min_read_len to primerscheme
--
ALTER TABLE `reads_primerscheme` ADD COLUMN `min_read_len` integer DEFAULT 400 NOT NULL;
ALTER TABLE `reads_primerscheme` ALTER COLUMN `min_read_len` DROP DEFAULT;
--
-- Alter field primer_scheme on jobmaster
--
ALTER TABLE `reads_jobmaster` DROP FOREIGN KEY `reads_jobmaster_primer_scheme_id_64932f56_fk_reads_pri`;
ALTER TABLE `reads_jobmaster` ADD CONSTRAINT `reads_jobmaster_primer_scheme_id_64932f56_fk_reads_pri` FOREIGN KEY (`primer_scheme_id`) REFERENCES `reads_primerscheme` (`id`);
