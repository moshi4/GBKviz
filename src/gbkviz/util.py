import os
import shutil
import time
from pathlib import Path
from typing import List

from streamlit.script_run_context import get_script_run_ctx
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


def get_session_id() -> str:
    """Get session id

    Returns:
        str: Session id
    """
    ctx = get_script_run_ctx()
    if ctx:
        return ctx.session_id
    else:
        raise ValueError("Failed to get session id.")


def make_session_dir(outdir: Path) -> Path:
    """Make session directory

    Args:
        outdir (Path): Output directory

    Returns:
        Path: Session directory path
    """
    session_dir = outdir / get_session_id()
    os.makedirs(session_dir, exist_ok=True)
    return session_dir


def remove_olddir(target_dir: Path, ttl: int = 600) -> None:
    """Remove old last updated directory

    Args:
        target_dir (Path): Target directory
        ttl (int, optional): Time to live. Defaults to 600[s].
    """
    elapsed_time = time.time() - target_dir.stat().st_mtime
    if elapsed_time > ttl:
        shutil.rmtree(target_dir, ignore_errors=True)
