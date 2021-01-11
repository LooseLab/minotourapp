--
-- Create model NotificationConditionsBarcode
--
CREATE TABLE `communication_notificationconditionsbarcode` (`id` integer AUTO_INCREMENT NOT NULL PRIMARY KEY, `barcode_name` varchar(256) NULL, `barcode_id` integer NOT NULL, `condition_id` integer NOT NULL);
ALTER TABLE `communication_notificationconditionsbarcode` ADD CONSTRAINT `communication_notifi_barcode_id_c78a8425_fk_reads_bar` FOREIGN KEY (`barcode_id`) REFERENCES `reads_barcode` (`id`);
ALTER TABLE `communication_notificationconditionsbarcode` ADD CONSTRAINT `communication_notifi_condition_id_3130a6a7_fk_communica` FOREIGN KEY (`condition_id`) REFERENCES `communication_notificationconditions` (`id`);
