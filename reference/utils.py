"""
Util functions that don't belong in any one view
"""
import gzip
import hashlib
from pathlib import Path

from django.conf import settings
from django.core.files.uploadedfile import TemporaryUploadedFile
from django.db.models import Q

from reference.models import ReferenceInfo


def check_media_roots_exist():
    """
    Check the django specified media root exists, and if not creates them
    Returns
    -------

    """
    media_root = Path(settings.MEDIA_ROOT)

    if not media_root.exists():
        media_root.mkdir()
    if not (media_root / "minimap2_indexes").exists():
        (media_root / "minimap2_indexes").mkdir()
    if not (media_root / "reference_files").exists():
        (media_root / "reference_files").mkdir()


def validate_reference_checks(ref_file, user):
    """

    Parameters
    ----------
    ref_file: pathlib.PosixPath
        The reference file path
    user: django.contrib.auth.models.User
        The user uploading the Reference
    Returns
    -------

    """
    """
    Perform all the checks that we need to do to.
    1. If the reference filename exists in any of the references available to the user, reject it
    2. Check if the md5 of the gzipped file matches any references, if it does, silently backref to that rather
     than add a new file.
    """
    print("Generating SHA Hash...")
    name = ref_file
    # If uploaded from web, could be a InMemoryFile, TempUploadedFile, so need a PosixPath of file name
    if not isinstance(ref_file, Path):
        name = Path(ref_file.name)
        # If TempUploadedFile, we need the file path to the temp on disk file
        if isinstance(ref_file, TemporaryUploadedFile):
            ref_file = Path(ref_file.temporary_file_path())
    read_handle = gzip.open if set(name.suffixes).intersection({".gz"}) else open
    with read_handle(
        ref_file,
        "rb",
    ) as f:
        # Read and update hash string value in blocks of 4K
        sha256_hash = hashlib.sha256()
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
        hash_string = sha256_hash.hexdigest()
        queryset = ReferenceInfo.objects.filter(sha256_checksum=hash_string).filter(
            Q(private=False) | Q(uploader=user)
        )
        if queryset:
            print(
                f"Exact reference already exists, please use {queryset.values_list('name', flat=True)[0]}"
            )
            return True, queryset.values_list('name')[0]
        return False, hash_string
