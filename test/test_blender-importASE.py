import os
import shutil
import subprocess

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))


def _find_blender():
    for cand in (os.environ.get("BLENDER"), "blender",
                 "/Applications/Blender.app/Contents/MacOS/Blender"):
        if cand and (shutil.which(cand) or os.path.exists(cand)):
            return shutil.which(cand) or cand
    return None


@pytest.mark.skipif(_find_blender() is None, reason="blender not on PATH")
def test_import_filepath_creates_objects():
    blender = _find_blender()
    script = os.path.join(HERE, "run_in_blender.py")
    out = subprocess.run(
        [blender, "--background", "--python", script],
        capture_output=True, text=True, timeout=120,
    )
    assert out.returncode == 0, f"stdout:\n{out.stdout}\nstderr:\n{out.stderr}"
    assert "OK" in out.stdout, out.stdout
