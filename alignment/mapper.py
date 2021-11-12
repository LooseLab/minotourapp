"""Mapping class, holds mappy aligner objects open for use around minotour in a threaded class"""
import time

import mappy as mp
from celery.utils.log import get_task_logger

logger = get_task_logger(__name__)


class MappingServer:

    """
    This class exists to provide a stable list of references to map against.
    The class should exist in perpetuity and ensure that we only ever have one instance of a reference loaded in memory.

    Attributes
    ----------
    references : set of str
        The reference names that are loaded.
    mapping_objects : dict
        Dictionary of mp.Aligner classes, keyed to reference and the last time the aligner was used
    interval: int
        The amount of time that a reference going unused will be tolerated before it is closed.
    """
    def __init__(self):
        self.references = {""}
        self.mapping_objects = {}
        self.interval = 1800  # Check every 30 minutes to see if a reference has been used - delete it if it hasn't been used.

    def reference_monitor(self):
        """
        Check every 30 minutes (Called by django celery beat) to see if a reference has been used in the previosu 30 minutes.
        If not, unload it.
        Returns
        -------

        """
        logger.info("checking references")
        logger.info(self.references)
        loop_me_baby = dict(self.mapping_objects)
        for reference in loop_me_baby:
            logger.info(reference)
            if (
                    self.mapping_objects[reference]["last_used"]
                    < time.time() - self.interval
            ):
                logger.info(f"This reference is old {reference}")
                self.delete_reference(reference)
                self.mapping_objects.pop(reference)
            else:
                logger.info(f"This reference is OK: {reference}")
            
    def add_reference(self, reference, filepath):
        """
        Add a reference to the set of reference names, and load the mappy index.
        Parameters
        ----------
        reference: str
            The string name of the reference
        filepath: str
            The file path to the reference MMI

        Returns
        -------

        """
        if reference not in self.references:
            logger.info("loading reference")
            self.references.add(reference)
            self.load_index(reference, str(filepath))
        else:
            logger.info("Already loaded")
        logger.info(self.references)

    def delete_reference(self, reference):
        """
        Delete a reference name from the reference set.
        Parameters
        ----------
        reference: str
            Name of the reference

        Returns
        -------

        """
        self.references.remove(reference)

    def list_references(self):
        """
        List the references held in the set of reference names.
        Returns
        -------
        set
            Set of reference names.
        """
        return self.references

    def valid(self, ref_name):
        """
        Check if a reference name is held in the reference set
        Parameters
        ----------
        ref_name: str
            Name of the reference.

        Returns
        -------
        bool
            True if reference is contained in references set, else false
        """
        if ref_name in self.references:
            logger.info(self.mapping_objects)
            return True
        else:
            return False

    def load_index(self, reference, filepath):
        """
        Load an index to be held in a mp.Aligner class
        Parameters
        ----------
        reference : str
            The name of the reference that is to be mapped against
        filepath : str
            The file path pointing to the minimap2 index of the reference we wish to map against

        Returns
        -------
        None
        """
        self.mapping_objects[reference] = {}
        self.mapping_objects[reference]["reference"] = mp.Aligner(
            filepath, preset="map-ont", best_n=1
        )
        self.mapping_objects[reference]["last_used"] = time.time()

    def map_sequence(self, reference, read_tuple, include_barcode=False):
        """
        This is a fast mapper that takes a sequence and returns the mapped sequence.
        Parameters
        ----------
        reference: str
            The name of the reference to be mapped against
        read_tuple: dict
            The read split into different dictionary fields
        include_barcode: bool, default False
            Include the name of the barcode in the paf results in the last index
        Returns
        -------
        list of str
            List of Paf formatted mapping hits
        """
        # For each alignment to be mapped against, returns a PAF format line
        result_me_up_inside = []
        for hit in self.mapping_objects[reference]["reference"].map(read_tuple["sequence"]):
            results = [read_tuple['read_id'], read_tuple['sequence']]
            split_hit = str(hit).split("\t")
            if include_barcode:
                barcode = read_tuple.get("barcode", False)
                split_hit.append("no_barcode" if not barcode else barcode)
            results.extend(
                split_hit
            )
            result_me_up_inside.append(results)
        self.refresh_index(reference)
        return result_me_up_inside

    def refresh_index(self, reference):
        """
        Update the last time that this reference was used
        Parameters
        ----------
        reference: str
            Name of the reference to be updated

        Returns
        -------
        None
        """
        self.mapping_objects[reference]["last_used"] = time.time()


MAP = MappingServer()
