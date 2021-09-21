--
-- Create model CompletedBarcodes
--
CREATE TABLE `readfish_completedbarcodes` (`id` bigint AUTO_INCREMENT NOT NULL PRIMARY KEY, `barcode` varchar(10) NOT NULL, `run_id` varchar(64) NOT NULL);
--
-- Create model TomlFile
--
CREATE TABLE `readfish_tomlfile` (`id` bigint AUTO_INCREMENT NOT NULL PRIMARY KEY, `toml_json` longtext NOT NULL, `sha_hash` varchar(64) NOT NULL, `run_id` integer NOT NULL);
--
-- Create constraint unique completed barcode on model completedbarcodes
--
ALTER TABLE `readfish_completedbarcodes` ADD CONSTRAINT `unique completed barcode` UNIQUE (`barcode`, `run_id`);
ALTER TABLE `readfish_tomlfile` ADD CONSTRAINT `readfish_tomlfile_run_id_b404cd77_fk_minknow_data_run_id` FOREIGN KEY (`run_id`) REFERENCES `minknow_data_run` (`id`);