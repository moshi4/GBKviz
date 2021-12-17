from pathlib import Path
from typing import List

from streamlit.uploaded_file_manager import UploadedFile, UploadedFileRec


def load_files(files: List[Path]) -> List[UploadedFile]:
    """Load files as uploaded files

    Args:
        files (List[Path]): files

    Returns:
        List[UploadedFile]: Uploaded files
    """
    uploaded_files: List[UploadedFile] = []
    for id, file in enumerate(files, 1):
        filename = Path(file).name
        type = "application/octet-stream"
        with open(file, "rb") as f:
            filebytes = f.read()
        record = UploadedFileRec(id, filename, type, filebytes)
        uploaded_files.append(UploadedFile(record))

    return uploaded_files
