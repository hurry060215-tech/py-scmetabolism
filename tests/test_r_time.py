import subprocess
import time
import os

rscript_path = "C:/Program Files/R/R-4.5.2/bin/Rscript.exe"
script_path = "tests/R_scripts/run_aucell_real_genes.R"

if os.path.exists(rscript_path) and os.path.exists(script_path):
    start = time.perf_counter()
    try:
        result = subprocess.run(
            [rscript_path, script_path],
            capture_output=True,
            text=True,
            cwd=".",  # root of project
            timeout=30
        )
        elapsed = time.perf_counter() - start
        print(f"R AUCell script completed in {elapsed:.2f} seconds")
        print(f"Return code: {result.returncode}")
        if result.stdout:
            print("Output (first 200 chars):", result.stdout[:200])
    except subprocess.TimeoutExpired:
        print("R script timed out after 30 seconds")
    except Exception as e:
        print(f"Error: {e}")
else:
    print("Rscript or script not found")