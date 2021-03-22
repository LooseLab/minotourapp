--
-- Alter field species on mappingresult
--
ALTER TABLE `metagenomics_mappingresult` MODIFY `species` varchar(256) NOT NULL;
--
-- Alter field species on mappingtarget
--
ALTER TABLE `metagenomics_mappingtarget` MODIFY `species` varchar(100) NULL;
