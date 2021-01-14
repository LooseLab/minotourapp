--
-- Rename field ninety_nine_percent_bases_at on articfireconditions to x_coverage
--
ALTER TABLE `artic_articfireconditions` CHANGE `ninety_nine_percent_bases_at` `x_coverage` integer NOT NULL;
--
-- Remove field ninety_five_percent_bases_at from articfireconditions
--
ALTER TABLE `artic_articfireconditions` DROP COLUMN `ninety_five_percent_bases_at`;
--
-- Remove field ninety_percent_bases_at from articfireconditions
--
ALTER TABLE `artic_articfireconditions` DROP COLUMN `ninety_percent_bases_at`;
--
-- Add field percent_of_amplicons to articfireconditions
--
ALTER TABLE `artic_articfireconditions` ADD COLUMN `percent_of_amplicons` integer DEFAULT 90 NOT NULL;
ALTER TABLE `artic_articfireconditions` ALTER COLUMN `percent_of_amplicons` DROP DEFAULT;
