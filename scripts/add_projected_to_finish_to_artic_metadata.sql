--
-- Add field projected_to_finish to articbarcodemetadata
--
ALTER TABLE `artic_articbarcodemetadata` ADD COLUMN `projected_to_finish` bool DEFAULT b'0' NOT NULL;
ALTER TABLE `artic_articbarcodemetadata` ALTER COLUMN `projected_to_finish` DROP DEFAULT;
--
-- Alter field x_coverage on articfireconditions
--
