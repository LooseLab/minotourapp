--
-- Add field artic_fire_conditions to articbarcodemetadata
--
CREATE TABLE `artic_articbarcodemetadata_artic_fire_conditions` (`id` integer AUTO_INCREMENT NOT NULL PRIMARY KEY, `articbarcodemetadata_id` integer NOT NULL, `articfireconditions_id` integer NOT NULL);
ALTER TABLE `artic_articbarcodemetadata_artic_fire_conditions` ADD CONSTRAINT `artic_articbarcodemetada_articbarcodemetadata_id__f2e7a042_uniq` UNIQUE (`articbarcodemetadata_id`, `articfireconditions_id`);
ALTER TABLE `artic_articbarcodemetadata_artic_fire_conditions` ADD CONSTRAINT `artic_articbarcodeme_articbarcodemetadata_7bf1d4a0_fk_artic_art` FOREIGN KEY (`articbarcodemetadata_id`) REFERENCES `artic_articbarcodemetadata` (`id`);
ALTER TABLE `artic_articbarcodemetadata_artic_fire_conditions` ADD CONSTRAINT `artic_articbarcodeme_articfireconditions__61a2c562_fk_artic_art` FOREIGN KEY (`articfireconditions_id`) REFERENCES `artic_articfireconditions` (`id`);
