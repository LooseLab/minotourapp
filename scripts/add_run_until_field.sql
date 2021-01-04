--
-- Add field run_until to notificationconditions
--
ALTER TABLE `communication_notificationconditions` ADD COLUMN `run_until` bool DEFAULT b'0' NOT NULL;
ALTER TABLE `communication_notificationconditions` ALTER COLUMN `run_until` DROP DEFAULT;
