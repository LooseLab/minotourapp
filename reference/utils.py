"""
Util functions that don't belong in any one view
"""
import gzip
import hashlib
import subprocess
from pathlib import Path

from django.core.files.uploadedfile import TemporaryUploadedFile, InMemoryUploadedFile
from django.db.models import Q

from minotourapp.settings import MINIMAP2, MEDIA_ROOT
from minotourapp.utils import get_env_variable
from reference.models import ReferenceInfo


def create_minimap2_index(ref_info, file_name):
    """
    Create the minimap2 index for the reference file that we are uploading
    Parameters
    ----------
    ref_info: reference.models.ReferenceInfo
        The django ORM object
    file_name: pathlib.PosixPath
        The filename we are uploading

    Returns
    -------
    str
        File path to the newly created minimap2 index file
    """
    index_dir_path = (
        MEDIA_ROOT
        if get_env_variable("MT_MINIMAP2_INDEX_DIR").isdigit()
        else get_env_variable("MT_MINIMAP2_INDEX_DIR")
    )
    minimap2_index_file_location = (
        f"{index_dir_path}/minimap2_indexes/{file_name.stem}.mmi"
    )
    out, err = subprocess.Popen(
        f"{MINIMAP2} -d {minimap2_index_file_location}"
        f" {ref_info.file_location.path}".split()
    ).communicate()
    print(out)
    print(err)
    return minimap2_index_file_location


def check_for_duplicate_references(hash_string, user):
    """
    Check for duplicate files using the hash string and user
    Parameters
    ----------
    hash_string: str
        Sha256 hash of the file
    user: django.contrib.auth.models.User
        The user uploading the file

    Returns
    -------
    boolean, tuple

    """
    queryset = ReferenceInfo.objects.filter(sha256_checksum=hash_string).filter(
        Q(private=False) | Q(uploader=user)
    )
    if queryset:
        print(
            f"Exact reference already exists, please use {queryset.values_list('name', flat=True)[0]}"
        )
        return True, queryset.values_list("name")[0]
    else:
        return False, []


def generate_sha256_hash(ref_file, user):
    """
    Generate the hash

    Parameters
    ----------
    ref_file: InMemoryUploadedFile or TemporaryUploadedFile
        The reference file being uploaded
    user: django.contrib.auth.models.User
        The user uploading the reference
    Returns
    -------
    boolean, str
    """
    sha256_hash = hashlib.sha256()
    for byte_block in iter(lambda: ref_file.read(4096), b""):
        sha256_hash.update(byte_block)
    hash_string = sha256_hash.hexdigest()
    duplicated, duplicated_file = check_for_duplicate_references(hash_string, user=user)

    return duplicated, hash_string


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
    print(f"Generating SHA Hash for {ref_file}...")
    name = ref_file
    # If uploaded from web, could be a InMemoryFile, TempUploadedFile, so need a PosixPath of file name
    if not isinstance(ref_file, Path):
        name = Path(ref_file.name)
        # If TempUploadedFile, we need the file path to the temp on disk file
        if isinstance(ref_file, TemporaryUploadedFile):
            ref_file = Path(ref_file.temporary_file_path())
    if set(name.suffixes).intersection({".gz"}):
        read_handle = gzip.open
    elif isinstance(ref_file, InMemoryUploadedFile):
        duplicated, hash_string = generate_sha256_hash(ref_file, user)
        return duplicated, hash_string
    else:
        read_handle = open
    with read_handle(ref_file, "rb",) as f:
        duplicated, hash_string = generate_sha256_hash(f, user)
        return duplicated, hash_string
