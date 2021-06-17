--
-- Add field final_completion to articbarcodemetadata
--
ALTER TABLE `artic_articbarcodemetadata` ADD COLUMN `final_completion` bool DEFAULT b'0' NOT NULL;
ALTER TABLE `artic_articbarcodemetadata` ALTER COLUMN `final_completion` DROP DEFAULT;
