--
-- Add field alt_name to barcode
--
ALTER TABLE `reads_barcode` ADD COLUMN `alt_name` varchar(128) NULL;