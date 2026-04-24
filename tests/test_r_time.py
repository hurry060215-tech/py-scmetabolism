import shutil
import subprocess
import pytest


@pytest.mark.skipif(not shutil.which("Rscript"), reason="Rscript not found on PATH")
def test_rscript_available():
    result = subprocess.run(
        ["Rscript", "--version"],
        capture_output=True,
        text=True,
        timeout=10,
    )
    assert result.returncode == 0
