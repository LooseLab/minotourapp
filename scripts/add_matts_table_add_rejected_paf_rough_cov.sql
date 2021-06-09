--
-- Add field rejected_barcode to pafroughcovintermediate
--
ALTER TABLE `alignment_pafroughcovintermediate` ADD COLUMN `rejected_barcode_id` integer NULL , ADD CONSTRAINT `alignment_pafroughco_rejected_barcode_id_9fc074ec_fk_reads_bar` FOREIGN KEY (`rejected_barcode_id`) REFERENCES `reads_barcode`(`id`);
--
-- Create model MattsAmazingAlignmentSum
--
CREATE TABLE `alignment_mattsamazingalignmentsum` (`id` integer AUTO_INCREMENT NOT NULL PRIMARY KEY, `is_pass` bool NOT NULL, `bin_position_start_str` longtext NOT NULL, `bin_coverage_str` longtext NOT NULL, `bin_window` integer NOT NULL, `reference_pk` integer NOT NULL, `reference_name` varchar(256) NULL, `chromosome_length` integer NOT NULL, `chromosome_pk` integer NOT NULL, `chromosome_name` varchar(256) NULL, `barcode_id` integer NULL, `chromosome_id` integer NULL, `flowcell_id` integer NULL, `job_master_id` integer NOT NULL, `read_type_id` integer NULL, `reference_id` integer NULL, `rejected_barcode_id` integer NULL, `run_id` integer NULL);
ALTER TABLE `alignment_mattsamazingalignmentsum` ADD CONSTRAINT `alignment_mattsamazi_barcode_id_4cd5fd91_fk_reads_bar` FOREIGN KEY (`barcode_id`) REFERENCES `reads_barcode` (`id`);
ALTER TABLE `alignment_mattsamazingalignmentsum` ADD CONSTRAINT `alignment_mattsamazi_chromosome_id_4346ffe6_fk_reference` FOREIGN KEY (`chromosome_id`) REFERENCES `reference_referenceline` (`id`);
ALTER TABLE `alignment_mattsamazingalignmentsum` ADD CONSTRAINT `alignment_mattsamazi_flowcell_id_fde2caee_fk_minknow_d` FOREIGN KEY (`flowcell_id`) REFERENCES `minknow_data_flowcell` (`id`);
ALTER TABLE `alignment_mattsamazingalignmentsum` ADD CONSTRAINT `alignment_mattsamazi_job_master_id_b7c9dae6_fk_reads_job` FOREIGN KEY (`job_master_id`) REFERENCES `reads_jobmaster` (`id`);
ALTER TABLE `alignment_mattsamazingalignmentsum` ADD CONSTRAINT `alignment_mattsamazi_read_type_id_a7839ebb_fk_reads_fas` FOREIGN KEY (`read_type_id`) REFERENCES `reads_fastqreadtype` (`id`);
ALTER TABLE `alignment_mattsamazingalignmentsum` ADD CONSTRAINT `alignment_mattsamazi_reference_id_c0e19a16_fk_reference` FOREIGN KEY (`reference_id`) REFERENCES `reference_referenceinfo` (`id`);
ALTER TABLE `alignment_mattsamazingalignmentsum` ADD CONSTRAINT `alignment_mattsamazi_rejected_barcode_id_99b26e9b_fk_reads_bar` FOREIGN KEY (`rejected_barcode_id`) REFERENCES `reads_barcode` (`id`);
ALTER TABLE `alignment_mattsamazingalignmentsum` ADD CONSTRAINT `alignment_mattsamazi_run_id_719cc994_fk_minknow_d` FOREIGN KEY (`run_id`) REFERENCES `minknow_data_run` (`id`);
